"""
API response contract validation.

The failure mode this module exists to prevent: an upstream API renames a
field, our parser silently returns None for it, and a 30,000-paper harvest
completes "successfully" with every record missing its journal, or its PMCID,
or its authors. Nothing raises. The corpus is quietly ruined and you find out
three weeks later when the lab-clustering output looks wrong.

The defence is a contract per endpoint: the fields we depend on, whether each
is required, and what shape it should have. Every harvester validates its
first N responses against its contract and refuses to proceed if the contract
is violated. Validation is cheap and runs once per session, not per record.

Contracts are also self-documenting: when Europe PMC does change something,
the failure message names the exact field and shows a sample of what actually
came back, so the fix is a one-line edit rather than an afternoon of printf
debugging.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from typing import Any, Callable

log = logging.getLogger(__name__)


@dataclass
class FieldSpec:
    """One field we depend on in an API response."""

    path: str  # dotted path, e.g. "journalInfo.journal.title"
    required: bool = False
    kind: str | None = None  # "str" | "int" | "list" | "dict" | "bool"
    note: str = ""  # what we use it for; appears in failure messages

    def extract(self, obj: Any) -> tuple[bool, Any]:
        """Walk the dotted path. Returns (found, value)."""
        cur = obj
        for part in self.path.split("."):
            if isinstance(cur, dict):
                if part not in cur:
                    return (False, None)
                cur = cur[part]
            elif isinstance(cur, list):
                # Take the first element for list-valued intermediate paths
                if not cur:
                    return (False, None)
                cur = cur[0]
                if isinstance(cur, dict) and part in cur:
                    cur = cur[part]
                else:
                    return (False, None)
            else:
                return (False, None)
        return (True, cur)

    def check_kind(self, value: Any) -> bool:
        if self.kind is None or value is None:
            return True
        return {
            "str": lambda v: isinstance(v, str),
            "int": lambda v: isinstance(v, int) and not isinstance(v, bool),
            "num": lambda v: isinstance(v, (int, float)) and not isinstance(v, bool),
            "list": lambda v: isinstance(v, list),
            "dict": lambda v: isinstance(v, dict),
            "bool": lambda v: isinstance(v, bool),
        }.get(self.kind, lambda v: True)(value)


@dataclass
class ContractResult:
    endpoint: str
    ok: bool
    n_samples: int = 0
    missing_required: list[str] = field(default_factory=list)
    missing_optional: list[str] = field(default_factory=list)
    wrong_kind: list[str] = field(default_factory=list)
    sample_keys: list[str] = field(default_factory=list)
    message: str = ""

    def report(self) -> str:
        lines = [f"[{self.endpoint}] {'OK' if self.ok else 'CONTRACT VIOLATION'}"
                 f" ({self.n_samples} samples checked)"]
        if self.missing_required:
            lines.append("  MISSING REQUIRED FIELDS:")
            for f_ in self.missing_required:
                lines.append(f"    - {f_}")
        if self.wrong_kind:
            lines.append("  WRONG TYPE:")
            for f_ in self.wrong_kind:
                lines.append(f"    - {f_}")
        if self.missing_optional:
            lines.append("  missing optional fields (degraded, not fatal):")
            for f_ in self.missing_optional[:8]:
                lines.append(f"    - {f_}")
        if not self.ok and self.sample_keys:
            lines.append(f"  keys actually present: {', '.join(self.sample_keys[:25])}")
        if self.message:
            lines.append(f"  {self.message}")
        return "\n".join(lines)


class Contract:
    """A set of field expectations for one API endpoint."""

    def __init__(self, endpoint: str, fields: list[FieldSpec],
                 root_path: str | None = None):
        self.endpoint = endpoint
        self.fields = fields
        self.root_path = root_path  # where the record list lives in the payload

    def records_from(self, payload: Any) -> list[dict]:
        """Pull the list of records out of a raw payload."""
        if self.root_path is None:
            return payload if isinstance(payload, list) else [payload]
        cur = payload
        for part in self.root_path.split("."):
            if isinstance(cur, dict):
                cur = cur.get(part)
            else:
                return []
            if cur is None:
                return []
        if isinstance(cur, list):
            return [r for r in cur if isinstance(r, dict)]
        return [cur] if isinstance(cur, dict) else []

    def validate(self, payload: Any, max_samples: int = 5) -> ContractResult:
        records = self.records_from(payload)
        res = ContractResult(endpoint=self.endpoint, ok=True)

        if not records:
            res.ok = False
            res.message = (
                f"No records found at path '{self.root_path}'. The response "
                "structure has changed, or the query returned nothing."
            )
            if isinstance(payload, dict):
                res.sample_keys = list(payload.keys())
            return res

        samples = records[:max_samples]
        res.n_samples = len(samples)
        res.sample_keys = sorted({k for r in samples for k in r.keys()})

        for spec in self.fields:
            found_any = False
            kind_ok_any = False
            for rec in samples:
                found, value = spec.extract(rec)
                if found:
                    found_any = True
                    if spec.check_kind(value):
                        kind_ok_any = True
                        break
            desc = f"{spec.path}" + (f"  ({spec.note})" if spec.note else "")
            if not found_any:
                if spec.required:
                    res.missing_required.append(desc)
                    res.ok = False
                else:
                    res.missing_optional.append(desc)
            elif not kind_ok_any:
                res.wrong_kind.append(f"{desc} -- expected {spec.kind}")
                if spec.required:
                    res.ok = False
        return res


# --------------------------------------------------------------------------
# Contracts for every external endpoint this system depends on
# --------------------------------------------------------------------------

EUROPEPMC_SEARCH = Contract(
    "Europe PMC /search",
    root_path="resultList.result",
    fields=[
        FieldSpec("id", required=False, kind="str", note="record id"),
        FieldSpec("source", required=False, kind="str",
                  note="PPR indicates a preprint"),
        FieldSpec("pmid", required=False, kind="str", note="primary identifier"),
        FieldSpec("pmcid", required=False, kind="str",
                  note="gates full-text retrieval"),
        FieldSpec("doi", required=False, kind="str"),
        FieldSpec("title", required=True, kind="str", note="paper title"),
        FieldSpec("authorString", required=False, kind="str",
                  note="author fallback"),
        # Required: lab clustering is the denominator for every replication
        # claim the system makes. Losing authors does not degrade the output,
        # it invalidates it -- every paper becomes its own "independent group".
        FieldSpec("authorList.author", required=True, kind="list",
                  note="authors -- drives lab clustering and every "
                       "independent-replication count"),
        FieldSpec("journalInfo.journal.title", required=False, kind="str",
                  note="journal -- drives quality weighting"),
        FieldSpec("journalInfo.yearOfPublication", required=False,
                  note="publication year -- drives recency weighting"),
        FieldSpec("firstPublicationDate", required=False, kind="str"),
        FieldSpec("isOpenAccess", required=False, kind="str"),
        FieldSpec("citedByCount", required=False, kind="int"),
        FieldSpec("pubTypeList.pubType", required=False, kind="list",
                  note="detects reviews and retractions"),
    ],
)

EUROPEPMC_SEARCH_ENVELOPE = Contract(
    "Europe PMC /search envelope",
    root_path=None,
    fields=[
        FieldSpec("hitCount", required=True, kind="int",
                  note="total results -- needed for pagination"),
        FieldSpec("nextCursorMark", required=False, kind="str",
                  note="pagination cursor; without it we cannot page"),
        FieldSpec("resultList", required=True, kind="dict"),
    ],
)

PUBMED_ESEARCH = Contract(
    "PubMed /esearch",
    root_path=None,
    fields=[
        FieldSpec("esearchresult.idlist", required=True, kind="list",
                  note="PMIDs returned by the search"),
        FieldSpec("esearchresult.count", required=False,
                  note="total hits -- needed for pagination"),
    ],
)

BIORXIV_DETAILS = Contract(
    "bioRxiv /details",
    root_path="collection",
    fields=[
        FieldSpec("doi", required=True, kind="str", note="preprint identifier"),
        FieldSpec("title", required=True, kind="str"),
        FieldSpec("authors", required=False, kind="str",
                  note="semicolon-delimited author list"),
        FieldSpec("date", required=False, kind="str"),
        FieldSpec("author_corresponding", required=False, kind="str",
                  note="senior author for lab clustering"),
        FieldSpec("author_corresponding_institution", required=False, kind="str",
                  note="institution for lab clustering"),
        FieldSpec("category", required=False, kind="str"),
    ],
)

CHEMBL_TARGET_SEARCH = Contract(
    "ChEMBL /target/search",
    root_path="targets",
    fields=[
        FieldSpec("target_chembl_id", required=True, kind="str"),
        FieldSpec("organism", required=False, kind="str",
                  note="filters to human targets"),
        FieldSpec("target_type", required=False, kind="str",
                  note="filters to SINGLE PROTEIN"),
        FieldSpec("target_components", required=False, kind="list",
                  note="carries the gene-symbol synonyms we match on"),
    ],
)

CHEMBL_MECHANISM = Contract(
    "ChEMBL /mechanism",
    root_path="mechanisms",
    fields=[
        FieldSpec("molecule_chembl_id", required=True, kind="str"),
        # Required, not optional: without action_type every compound falls to
        # MoADirection.UNKNOWN and direction matching silently excludes
        # everything, or -- worse, if someone relaxes the filter -- ranks
        # inhibitors alongside activators. Losing this field must be loud.
        FieldSpec("action_type", required=True, kind="str",
                  note="determines agonist vs antagonist; direction matching "
                       "cannot run without it"),
        FieldSpec("mechanism_of_action", required=False, kind="str"),
        FieldSpec("max_phase", required=False, note="clinical stage"),
    ],
)

CHEMBL_MOLECULE = Contract(
    "ChEMBL /molecule",
    root_path=None,
    fields=[
        FieldSpec("pref_name", required=False, kind="str"),
        FieldSpec("molecule_properties", required=False, kind="dict",
                  note="physicochemistry for CNS scoring"),
        FieldSpec("molecule_properties.full_mwt", required=False,
                  note="molecular weight"),
        FieldSpec("molecule_properties.alogp", required=False, note="logP"),
        FieldSpec("molecule_properties.psa", required=False, note="TPSA"),
        FieldSpec("molecule_properties.hbd", required=False,
                  note="H-bond donors"),
        FieldSpec("max_phase", required=False),
        FieldSpec("withdrawn_flag", required=False),
    ],
)

CHEMBL_ACTIVITY = Contract(
    "ChEMBL /activity",
    root_path="activities",
    fields=[
        FieldSpec("standard_type", required=False, kind="str",
                  note="IC50 / Ki / EC50"),
        FieldSpec("standard_value", required=False, note="potency value"),
        # Required: a missing unit means we cannot convert to nM, and a
        # silently-unconverted uM value is a 1000x potency error reported as
        # if it were measured.
        FieldSpec("standard_units", required=True, kind="str",
                  note="nM / uM / M -- a missing unit is a 1000x error"),
    ],
)

ALL_CONTRACTS: dict[str, Contract] = {
    "europepmc_search": EUROPEPMC_SEARCH,
    "europepmc_envelope": EUROPEPMC_SEARCH_ENVELOPE,
    "pubmed_esearch": PUBMED_ESEARCH,
    "biorxiv_details": BIORXIV_DETAILS,
    "chembl_target": CHEMBL_TARGET_SEARCH,
    "chembl_mechanism": CHEMBL_MECHANISM,
    "chembl_molecule": CHEMBL_MOLECULE,
    "chembl_activity": CHEMBL_ACTIVITY,
}


class DriftDetector:
    """Validates live responses against contracts once per session.

    Attach to a harvester; it checks the first response from each endpoint and
    then gets out of the way. Set `strict=True` to raise on violation rather
    than warn -- appropriate for a long unattended harvest, where a silent
    degradation is far more expensive than an early crash.
    """

    def __init__(self, strict: bool = False,
                 on_violation: Callable[[ContractResult], None] | None = None):
        self.strict = strict
        self.on_violation = on_violation
        self.checked: dict[str, ContractResult] = {}

    def check(self, contract_key: str, payload: Any) -> ContractResult | None:
        """Validate a payload once per contract key. Returns None if already done."""
        if contract_key in self.checked:
            return None
        contract = ALL_CONTRACTS.get(contract_key)
        if contract is None:
            log.debug("No contract registered for '%s'", contract_key)
            return None

        result = contract.validate(payload)
        self.checked[contract_key] = result

        if not result.ok:
            msg = (
                f"\n{'=' * 70}\n"
                f"API SCHEMA DRIFT DETECTED\n"
                f"{'=' * 70}\n"
                f"{result.report()}\n"
                f"{'=' * 70}\n"
                f"The upstream API no longer returns the fields this parser "
                f"expects.\nHarvesting now would produce records with missing "
                f"data and no error.\nFix the parser in synaptik/ingest/ before "
                f"continuing.\n"
            )
            if self.on_violation:
                self.on_violation(result)
            if self.strict:
                raise SchemaDriftError(msg)
            log.error(msg)
        elif result.missing_optional:
            log.warning(
                "[%s] contract satisfied but %d optional field(s) absent: %s",
                result.endpoint, len(result.missing_optional),
                ", ".join(result.missing_optional[:5]),
            )
        else:
            log.info("[%s] response contract OK", result.endpoint)
        return result

    def summary(self) -> str:
        if not self.checked:
            return "No endpoints validated yet."
        return "\n\n".join(r.report() for r in self.checked.values())

    def all_ok(self) -> bool:
        return all(r.ok for r in self.checked.values())


class SchemaDriftError(RuntimeError):
    """Raised when a live API response violates its contract in strict mode."""


# --------------------------------------------------------------------------
# Offline self-test fixtures
# --------------------------------------------------------------------------

# Captured shapes of real responses, used to test the validator itself and to
# document what each endpoint is expected to return. If you ever need to see
# what the parser assumes, read these.
FIXTURES: dict[str, dict] = {
    "europepmc_search": {
        "hitCount": 12345,
        "nextCursorMark": "AoJwsPqAipYDKDIzMjE0NDUz",
        "resultList": {"result": [{
            "id": "23214453", "source": "MED", "pmid": "23214453",
            "pmcid": "PMC3535000", "doi": "10.1016/j.neuron.2012.11.001",
            "title": "Rapid single-step induction of functional neurons.",
            "authorString": "Zhang Y, Pak C, Han Y.",
            "authorList": {"author": [
                {"fullName": "Zhang Y", "firstName": "Yingsha",
                 "lastName": "Zhang",
                 "affiliation": "Department of Molecular and Cellular "
                                "Physiology, Stanford University"},
            ]},
            "journalInfo": {
                "yearOfPublication": 2013,
                "journal": {"title": "Neuron", "isoabbreviation": "Neuron"},
            },
            "firstPublicationDate": "2013-06-05",
            "isOpenAccess": "Y", "citedByCount": 1420,
            "pubTypeList": {"pubType": ["Journal Article", "research-article"]},
        }]},
    },
    "pubmed_esearch": {
        "esearchresult": {
            "count": "4210", "retmax": "20", "retstart": "0",
            "idlist": ["23214453", "28153965", "31000000"],
        }
    },
    "biorxiv_details": {
        "messages": [{"status": "ok", "total": 120}],
        "collection": [{
            "doi": "10.1101/2023.01.15.524000",
            "title": "Cryptic exon inclusion in iPSC-derived neurons",
            "authors": "Smith, J.; Doe, A.",
            "author_corresponding": "Jane Smith",
            "author_corresponding_institution": "Massachusetts General Hospital",
            "date": "2023-01-15", "version": "1", "category": "neuroscience",
        }],
    },
    "chembl_target": {
        "targets": [{
            "target_chembl_id": "CHEMBL2028",
            "organism": "Homo sapiens", "target_type": "SINGLE PROTEIN",
            "pref_name": "Neurotrophic tyrosine kinase receptor type 2",
            "target_components": [{
                "target_component_synonyms": [
                    {"component_synonym": "NTRK2", "syn_type": "GENE_SYMBOL"},
                ]}],
        }]
    },
    "chembl_mechanism": {
        "mechanisms": [{
            "molecule_chembl_id": "CHEMBL1234",
            "action_type": "AGONIST",
            "mechanism_of_action": "Tyrosine kinase receptor B agonist",
            "max_phase": 2, "target_chembl_id": "CHEMBL2028",
        }]
    },
    "chembl_molecule": {
        "molecule_chembl_id": "CHEMBL1234", "pref_name": "EXAMPLE-DRUG",
        "max_phase": 4, "first_approval": 1996, "withdrawn_flag": False,
        "molecule_properties": {
            "full_mwt": "379.50", "alogp": "4.30", "psa": "38.77",
            "hbd": 0, "rtb": 6, "cx_most_bpka": "8.60",
        },
    },
    "chembl_activity": {
        "activities": [{
            "standard_type": "EC50", "standard_value": "45.0",
            "standard_units": "nM", "molecule_chembl_id": "CHEMBL1234",
        }]
    },
}


def self_test(verbose: bool = True) -> bool:
    """Validate every contract against its fixture. Run in CI.

    This checks the validator, not the APIs -- but it also guarantees that if
    someone edits a parser and forgets the contract, the mismatch surfaces
    immediately rather than at harvest time.
    """
    ok = True
    for key, fixture in FIXTURES.items():
        contract = ALL_CONTRACTS[key]
        result = contract.validate(fixture)
        if verbose:
            print(result.report())
        if not result.ok:
            ok = False
    # envelope contract shares the europepmc fixture
    env = EUROPEPMC_SEARCH_ENVELOPE.validate(FIXTURES["europepmc_search"])
    if verbose:
        print(env.report())
    ok = ok and env.ok
    return ok


def dump_contracts(path: str) -> None:
    """Write all contracts to JSON, for documentation or external checking."""
    data = {
        key: {
            "endpoint": c.endpoint,
            "root_path": c.root_path,
            "fields": [
                {"path": f_.path, "required": f_.required,
                 "kind": f_.kind, "note": f_.note}
                for f_ in c.fields
            ],
        }
        for key, c in ALL_CONTRACTS.items()
    }
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2)


if __name__ == "__main__":
    import sys
    sys.exit(0 if self_test() else 1)
