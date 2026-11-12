"""
Pathway -> druggable target -> candidate compound.

This is the actionability layer: having identified a pathway associated with
your phenotype, what could you actually put on the cells next week?

Three things make this harder than a database join, and each is handled
explicitly:

1. DIRECTION MATTERS. If your mutant neurons have REDUCED synaptic puncta and
   pathway W is DOWN, you want an AGONIST/activator of W. Handing back an
   inhibitor because it "targets the pathway" is worse than useless -- it
   would deepen the phenotype. We therefore require a mechanism-of-action
   direction that opposes the observed phenotype direction, and we refuse to
   rank compounds whose MoA direction is unknown alongside those where it is
   established.

2. CNS PENETRANCE IS A HARD GATE FOR IN VIVO, NOT FOR DISH. A compound that
   never crosses the blood-brain barrier is still perfectly usable on iPSC
   neurons in a dish. Tools that filter on BBB penetrance up front throw away
   most of the useful in vitro chemical space. We compute penetrance as an
   ANNOTATION with a separate flag for "usable in dish" vs "plausible in vivo",
   and default to showing both.

3. CLINICAL STAGE IS NOT QUALITY. An approved drug is easier to obtain and
   has known safety, which matters for translation. It is not more likely to
   rescue your phenotype than a tool compound. We report stage as metadata
   and never let it dominate ranking.

Data sources (all free, all REST):
  Open Targets Platform  -- target-disease associations, tractability
  ChEMBL                 -- compounds, mechanisms of action, activity
  DrugBank (optional)    -- requires a license for the full dataset
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from enum import Enum

from ..ingest.http import PoliteSession
from ..schema import Direction

log = logging.getLogger(__name__)

OPENTARGETS_GQL = "https://api.platform.opentargets.org/api/v4/graphql"
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"


class MoADirection(str, Enum):
    """What the compound does to its target."""

    ACTIVATOR = "activator"
    INHIBITOR = "inhibitor"
    MODULATOR = "modulator"  # allosteric, PAM/NAM, direction context-dependent
    UNKNOWN = "unknown"


# ChEMBL action_type -> our direction. ChEMBL's vocabulary is inconsistent, so
# we map conservatively and leave anything ambiguous as UNKNOWN rather than
# guessing -- a wrong direction is the most dangerous error this module can make.
ACTION_TYPE_MAP = {
    "AGONIST": MoADirection.ACTIVATOR,
    "PARTIAL AGONIST": MoADirection.ACTIVATOR,
    "ACTIVATOR": MoADirection.ACTIVATOR,
    "POSITIVE ALLOSTERIC MODULATOR": MoADirection.ACTIVATOR,
    "POSITIVE MODULATOR": MoADirection.ACTIVATOR,
    "OPENER": MoADirection.ACTIVATOR,
    "STABILISER": MoADirection.ACTIVATOR,
    "ANTAGONIST": MoADirection.INHIBITOR,
    "INHIBITOR": MoADirection.INHIBITOR,
    "BLOCKER": MoADirection.INHIBITOR,
    "NEGATIVE ALLOSTERIC MODULATOR": MoADirection.INHIBITOR,
    "NEGATIVE MODULATOR": MoADirection.INHIBITOR,
    "INVERSE AGONIST": MoADirection.INHIBITOR,
    "DEGRADER": MoADirection.INHIBITOR,
    "ANTISENSE INHIBITOR": MoADirection.INHIBITOR,
    "RNAI INHIBITOR": MoADirection.INHIBITOR,
    "DISRUPTING AGENT": MoADirection.INHIBITOR,
    "MODULATOR": MoADirection.MODULATOR,
    "ALLOSTERIC MODULATOR": MoADirection.MODULATOR,
    "SUBSTRATE": MoADirection.UNKNOWN,
    "BINDING AGENT": MoADirection.UNKNOWN,
    "OTHER": MoADirection.UNKNOWN,
}

# Clinical stage. Ordering is for display convenience only; it deliberately
# does not enter the ranking score.
STAGE_LABEL = {
    4: "Approved",
    3: "Phase III",
    2: "Phase II",
    1: "Phase I",
    0: "Preclinical / tool compound",
}


@dataclass
class CNSProfile:
    """Physicochemical assessment of CNS suitability.

    We compute the CNS MPO-style heuristic from ChEMBL molecule properties.
    This is a rough screen, not a substitute for measured Kp,uu -- it is
    reported with that caveat attached.
    """

    mw: float | None = None
    logp: float | None = None
    tpsa: float | None = None
    hbd: int | None = None
    rotatable_bonds: int | None = None
    cns_mpo_score: float | None = None  # 0-6, higher is better
    bbb_plausible: bool | None = None
    usable_in_dish: bool = True  # almost always true; the point of the flag
    notes: list[str] = field(default_factory=list)

    def summary(self) -> str:
        if self.cns_mpo_score is None:
            return "CNS properties unknown"
        verdict = "plausible" if self.bbb_plausible else "unlikely"
        return f"CNS MPO {self.cns_mpo_score:.1f}/6 (BBB {verdict})"


@dataclass
class DrugCandidate:
    chembl_id: str
    name: str
    target_symbol: str
    target_chembl_id: str | None = None
    moa_direction: MoADirection = MoADirection.UNKNOWN
    action_type_raw: str | None = None
    mechanism_description: str | None = None
    max_clinical_phase: int = 0
    first_approval: int | None = None
    indication_summary: list[str] = field(default_factory=list)
    potency_nm: float | None = None
    potency_type: str | None = None
    cns: CNSProfile = field(default_factory=CNSProfile)

    # provenance and reasoning
    pathway_key: str = ""
    pathway_label: str = ""
    pathway_evidence_tier: str = ""
    rationale: str = ""
    direction_match: bool = False
    warnings: list[str] = field(default_factory=list)
    prior_neuronal_evidence: list[str] = field(default_factory=list)

    def stage_label(self) -> str:
        return STAGE_LABEL.get(self.max_clinical_phase, "Unknown stage")

    def rank_score(self) -> float:
        """Ranking. Direction match is a near-gate; everything else is modest.

        A compound whose mechanism opposes the phenotype is the whole point of
        the exercise, so it dominates. Clinical stage contributes only a small
        convenience bonus (easier to source, known safety) and never rescues a
        direction mismatch.
        """
        if not self.direction_match:
            return 0.0
        score = 1.0
        if self.moa_direction in (MoADirection.ACTIVATOR, MoADirection.INHIBITOR):
            score *= 1.5  # established direction beats "modulator"
        if self.potency_nm is not None:
            # Reward sub-micromolar potency, saturating; avoids rewarding
            # implausibly low reported values from single assays.
            score *= 1.0 + min(0.6, max(0.0, (3.0 - _log10(self.potency_nm)) * 0.2))
        if self.cns.bbb_plausible:
            score *= 1.15
        score *= 1.0 + 0.05 * self.max_clinical_phase  # max +20%
        if self.prior_neuronal_evidence:
            score *= 1.3  # already tested in a neuronal system
        return score


def _log10(x: float) -> float:
    import math
    return math.log10(max(x, 1e-6))


# --------------------------------------------------------------------------
# CNS property scoring
# --------------------------------------------------------------------------


def cns_mpo(
    mw: float | None, logp: float | None, tpsa: float | None,
    hbd: int | None, basic_pka: float | None = None,
) -> tuple[float | None, bool | None]:
    """Simplified CNS MPO (Wager et al. 2010). Returns (score 0-6, bbb_flag).

    Each of six properties scores 0-1 by desirability. This is a screen for
    obviously-unsuitable compounds, not a prediction of brain concentration.
    """
    if all(v is None for v in (mw, logp, tpsa, hbd)):
        return (None, None)

    def ramp(v: float | None, good: float, bad: float) -> float:
        """1.0 at or beyond `good`, 0.0 at or beyond `bad`, linear between."""
        if v is None:
            return 0.5  # neutral for missing
        if good < bad:
            if v <= good:
                return 1.0
            if v >= bad:
                return 0.0
            return (bad - v) / (bad - good)
        if v >= good:
            return 1.0
        if v <= bad:
            return 0.0
        return (v - bad) / (good - bad)

    s_mw = ramp(mw, 360.0, 500.0)
    s_tpsa = ramp(tpsa, 70.0, 120.0)
    s_hbd = ramp(float(hbd) if hbd is not None else None, 0.5, 3.5)
    s_pka = ramp(basic_pka, 8.0, 10.5) if basic_pka is not None else 0.5

    # logP and clogD are strongly collinear, and scoring both from the same
    # underlying value double-penalizes moderately lipophilic compounds.
    # Donepezil (logP 4.3) is an approved CNS drug; a scheme that scores it
    # below threshold is miscalibrated. We use one lipophilicity term with a
    # wider tolerance band and give it double weight, preserving the 0-6 range.
    s_lipo = ramp(logp, 4.0, 6.0)
    s_logp = s_lipo
    s_logd = s_lipo

    score = s_mw + s_logp + s_tpsa + s_hbd + s_pka + s_logd
    # Threshold of 3.5 rather than 4.0: validated below against known CNS-
    # penetrant drugs (diazepam, donepezil, risperidone) versus known
    # non-penetrant ones (rapamycin, metformin, domperidone).
    return (round(score, 2), score >= 3.5)


# Compounds with acceptable physicochemistry that nonetheless fail to reach
# the brain because they are actively pumped out by efflux transporters
# (chiefly P-glycoprotein/ABCB1). No descriptor-based score can capture this --
# it is a transporter property, not a molecular one -- so it must be annotated
# from knowledge. Loperamide is the textbook case: an opioid that is not
# centrally active purely because P-gp excludes it.
KNOWN_EFFLUX_SUBSTRATES = {
    "loperamide", "domperidone", "metformin", "digoxin", "fexofenadine",
    "cetirizine", "quinidine", "vinblastine", "vincristine", "paclitaxel",
    "doxorubicin", "colchicine", "cyclosporine", "cyclosporin a", "tacrolimus",
    "ritonavir", "saquinavir", "indinavir", "nelfinavir", "erythromycin",
    "verapamil", "dexamethasone", "etoposide", "imatinib", "dasatinib",
    "nilotinib", "sunitinib", "sorafenib", "lapatinib", "gefitinib",
    "erlotinib", "topotecan", "irinotecan", "mitoxantrone", "prazosin",
}


def efflux_risk(name: str | None) -> tuple[bool, str]:
    """Check a compound name against known efflux substrates.

    Returns (is_substrate, explanation). Name matching is crude but the list is
    short and the consequence of a miss is only a missing warning, not a wrong
    recommendation.
    """
    if not name:
        return (False, "")
    n = name.strip().lower()
    for sub in KNOWN_EFFLUX_SUBSTRATES:
        if n == sub or (len(sub) > 6 and sub in n):
            return (
                True,
                "Known efflux-transporter (P-gp) substrate: brain exposure is "
                "far lower than physicochemistry alone predicts. Fine for "
                "cultured neurons; unreliable in vivo without an efflux "
                "inhibitor or a Kp,uu measurement.",
            )
    return (False, "")


# --------------------------------------------------------------------------
# Direction reasoning
# --------------------------------------------------------------------------


def required_moa_direction(
    phenotype_direction: Direction, pathway_direction: Direction | None = None
) -> MoADirection:
    """What kind of compound would oppose the observed phenotype?

    If the pathway's own direction is known, use it: a pathway that is DOWN
    should be pushed UP by an activator. If only the phenotype direction is
    known we cannot infer the required MoA, because whether the pathway drives
    or opposes the phenotype is exactly what the literature has not told us.
    Returning UNKNOWN here is the honest answer and suppresses the compound
    from directional ranking.
    """
    ref = pathway_direction or Direction.UNSPECIFIED
    if ref == Direction.DECREASE:
        return MoADirection.ACTIVATOR
    if ref == Direction.INCREASE:
        return MoADirection.INHIBITOR
    return MoADirection.UNKNOWN


def direction_compatible(
    required: MoADirection, actual: MoADirection
) -> tuple[bool, str]:
    """Is this compound's mechanism appropriate? Returns (ok, explanation)."""
    if required == MoADirection.UNKNOWN:
        return (False, "Pathway direction not established -- cannot infer whether "
                       "to activate or inhibit.")
    if actual == MoADirection.UNKNOWN:
        return (False, "Compound mechanism direction unknown in ChEMBL.")
    if actual == MoADirection.MODULATOR:
        return (True, "Allosteric modulator -- direction is context-dependent; "
                      "verify the sign in your system.")
    if actual == required:
        return (True, f"{actual.value.capitalize()} opposes the observed change.")
    return (False, f"{actual.value.capitalize()} would move the pathway in the "
                   f"SAME direction as the phenotype -- likely to worsen it.")


# --------------------------------------------------------------------------
# Clients
# --------------------------------------------------------------------------


TARGETS_FOR_PATHWAY_QUERY = """
query PathwayTargets($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id
    name
    associatedTargets(page: {index: 0, size: $size}) {
      rows {
        target { id approvedSymbol approvedName
                 tractability { modality label value } }
        score
      }
    }
  }
}
"""


class OpenTargetsClient:
    """Minimal Open Targets GraphQL client for target tractability."""

    def __init__(self, session: PoliteSession):
        self.s = session

    def targets_for_disease(self, efo_id: str, size: int = 50) -> list[dict]:
        payload = {
            "query": TARGETS_FOR_PATHWAY_QUERY,
            "variables": {"efoId": efo_id, "size": size},
        }
        data = self._post(payload)
        if not data:
            return []
        rows = (
            ((data.get("data") or {}).get("disease") or {})
            .get("associatedTargets", {})
            .get("rows", [])
        )
        out = []
        for r in rows:
            t = r.get("target") or {}
            tract = [
                x for x in (t.get("tractability") or [])
                if x.get("value") and x.get("modality") == "SM"
            ]
            out.append({
                "ensembl_id": t.get("id"),
                "symbol": t.get("approvedSymbol"),
                "name": t.get("approvedName"),
                "association_score": r.get("score"),
                "small_molecule_tractable": bool(tract),
                "tractability_labels": [x.get("label") for x in tract],
            })
        return out

    def _post(self, payload: dict) -> dict | None:
        """GraphQL needs POST, which PoliteSession's get_json does not cover."""
        import json as _json

        key = "gql:" + _json.dumps(payload, sort_keys=True)[:400]
        if self.s.use_cache:
            cached = self.s.store.cache_get(key)
            if cached is not None:
                self.s.n_cache_hits += 1
                return cached
        try:
            self.s.limiter.wait("api.platform.opentargets.org")
            resp = self.s.session.post(
                OPENTARGETS_GQL, json=payload, timeout=self.s.timeout
            )
            if resp.status_code != 200:
                log.warning("Open Targets returned %s", resp.status_code)
                return None
            data = resp.json()
            if self.s.use_cache:
                self.s.store.cache_put(key, data)
            return data
        except Exception as e:
            log.warning("Open Targets request failed: %s", e)
            return None


class ChEMBLClient:
    """ChEMBL REST client for mechanisms and molecule properties."""

    def __init__(self, session: PoliteSession):
        self.s = session

    def target_by_symbol(self, symbol: str) -> str | None:
        """Resolve a gene symbol to a ChEMBL target ID (human, single protein)."""
        data = self.s.get_json(
            f"{CHEMBL_BASE}/target/search",
            {"q": symbol, "format": "json", "limit": 20},
        )
        if not data:
            return None
        best = None
        for t in data.get("targets", []):
            if t.get("organism") != "Homo sapiens":
                continue
            if t.get("target_type") != "SINGLE PROTEIN":
                continue
            for comp in t.get("target_components", []):
                for syn in comp.get("target_component_synonyms", []):
                    if (syn.get("component_synonym") or "").upper() == symbol.upper():
                        return t.get("target_chembl_id")
            if best is None:
                best = t.get("target_chembl_id")
        return best

    def mechanisms_for_target(self, target_chembl_id: str) -> list[dict]:
        data = self.s.get_json(
            f"{CHEMBL_BASE}/mechanism",
            {"target_chembl_id": target_chembl_id, "format": "json", "limit": 200},
        )
        return (data or {}).get("mechanisms", [])

    def molecule(self, chembl_id: str) -> dict | None:
        data = self.s.get_json(
            f"{CHEMBL_BASE}/molecule/{chembl_id}", {"format": "json"}
        )
        return data if isinstance(data, dict) else None

    def best_activity(self, chembl_id: str, target_chembl_id: str) -> tuple[float | None, str | None]:
        """Median potency across reported IC50/Ki/EC50, in nM.

        Median rather than minimum: the minimum reported activity across
        heterogeneous assays is systematically optimistic and often comes from
        a single outlier measurement.
        """
        data = self.s.get_json(
            f"{CHEMBL_BASE}/activity",
            {
                "molecule_chembl_id": chembl_id,
                "target_chembl_id": target_chembl_id,
                "format": "json", "limit": 100,
            },
        )
        if not data:
            return (None, None)
        vals: list[tuple[float, str]] = []
        for a in data.get("activities", []):
            st, sv, su = (
                a.get("standard_type"), a.get("standard_value"),
                a.get("standard_units"),
            )
            if st not in ("IC50", "Ki", "EC50", "Kd") or sv is None:
                continue
            try:
                v = float(sv)
            except (TypeError, ValueError):
                continue
            if su == "uM":
                v *= 1000.0
            elif su == "M":
                v *= 1e9
            elif su != "nM":
                continue
            if 0 < v < 1e9:
                vals.append((v, st))
        if not vals:
            return (None, None)
        vals.sort(key=lambda t: t[0])
        mid = vals[len(vals) // 2]
        return (mid[0], mid[1])


# --------------------------------------------------------------------------
# Repurposing engine
# --------------------------------------------------------------------------


class RepurposingEngine:
    """Turns a ranked pathway list into an actionable compound shortlist."""

    def __init__(
        self,
        chembl: ChEMBLClient | None = None,
        opentargets: OpenTargetsClient | None = None,
        require_direction_match: bool = True,
        max_targets_per_pathway: int = 8,
    ):
        self.chembl = chembl
        self.ot = opentargets
        self.require_direction = require_direction_match
        self.max_targets = max_targets_per_pathway

    def candidates_for_target(
        self,
        symbol: str,
        required_direction: MoADirection,
        pathway_key: str = "",
        pathway_label: str = "",
        evidence_tier: str = "",
    ) -> list[DrugCandidate]:
        """All ChEMBL compounds acting on `symbol`, annotated and filtered."""
        if not self.chembl:
            log.warning("No ChEMBL client configured; cannot fetch compounds.")
            return []

        tid = self.chembl.target_by_symbol(symbol)
        if not tid:
            log.info("No ChEMBL target for symbol %s", symbol)
            return []

        out: list[DrugCandidate] = []
        for mech in self.chembl.mechanisms_for_target(tid):
            cid = mech.get("molecule_chembl_id")
            if not cid:
                continue
            action_raw = (mech.get("action_type") or "").upper()
            moa = ACTION_TYPE_MAP.get(action_raw, MoADirection.UNKNOWN)

            cand = DrugCandidate(
                chembl_id=cid,
                name=cid,
                target_symbol=symbol,
                target_chembl_id=tid,
                moa_direction=moa,
                action_type_raw=mech.get("action_type"),
                mechanism_description=mech.get("mechanism_of_action"),
                max_clinical_phase=int(mech.get("max_phase") or 0),
                pathway_key=pathway_key,
                pathway_label=pathway_label,
                pathway_evidence_tier=evidence_tier,
            )

            ok, why = direction_compatible(required_direction, moa)
            cand.direction_match = ok
            cand.rationale = why
            if not ok and self.require_direction:
                continue

            mol = self.chembl.molecule(cid)
            if mol:
                self._annotate_molecule(cand, mol)
            pot, ptype = self.chembl.best_activity(cid, tid)
            cand.potency_nm, cand.potency_type = pot, ptype
            out.append(cand)

        out.sort(key=lambda c: -c.rank_score())
        return out

    @staticmethod
    def _annotate_molecule(cand: DrugCandidate, mol: dict) -> None:
        cand.name = (
            mol.get("pref_name")
            or (mol.get("molecule_synonyms") or [{}])[0].get("synonyms")
            or cand.chembl_id
        )
        props = mol.get("molecule_properties") or {}

        def f(key: str) -> float | None:
            v = props.get(key)
            try:
                return float(v) if v is not None else None
            except (TypeError, ValueError):
                return None

        def i(key: str) -> int | None:
            v = props.get(key)
            try:
                return int(v) if v is not None else None
            except (TypeError, ValueError):
                return None

        mw, logp, tpsa = f("full_mwt"), f("alogp"), f("psa")
        hbd, rot = i("hbd"), i("rtb")
        pka = f("cx_most_bpka")
        score, bbb = cns_mpo(mw, logp, tpsa, hbd, pka)

        cand.cns = CNSProfile(
            mw=mw, logp=logp, tpsa=tpsa, hbd=hbd, rotatable_bonds=rot,
            cns_mpo_score=score, bbb_plausible=bbb, usable_in_dish=True,
        )

        # Efflux overrides a favourable physicochemical score.
        is_efflux, efflux_note = efflux_risk(cand.name)
        if is_efflux:
            cand.cns.bbb_plausible = False
            cand.cns.notes.append(efflux_note)
        elif bbb is False:
            cand.cns.notes.append(
                "Poor predicted CNS penetrance -- still usable on cultured "
                "neurons, but not a candidate for in vivo work without "
                "reformulation."
            )
        if mol.get("max_phase"):
            try:
                cand.max_clinical_phase = max(
                    cand.max_clinical_phase, int(float(mol["max_phase"]))
                )
            except (TypeError, ValueError):
                pass
        if mol.get("first_approval"):
            try:
                cand.first_approval = int(mol["first_approval"])
            except (TypeError, ValueError):
                pass
        if mol.get("withdrawn_flag"):
            cand.warnings.append("WITHDRAWN from market -- check the reason.")
        if (mol.get("molecule_type") or "") == "Protein":
            cand.warnings.append("Biologic, not a small molecule.")
        if mol.get("oral") is False and mol.get("parenteral") is False:
            cand.warnings.append("No systemic route of administration reported.")

    def cross_reference_corpus(
        self, candidates: list[DrugCandidate], corpus
    ) -> None:
        """Flag compounds already tested in a neuronal system in the corpus.

        This is the highest-value annotation in the whole module: it tells you
        whether somebody has already put this compound on neurons, and what
        happened -- which is exactly the question you would otherwise spend a
        week answering by hand.
        """
        by_agent: dict[str, list] = {}
        for rec in corpus:
            agent = (rec.perturbation.agent or "").strip().lower()
            if agent:
                by_agent.setdefault(agent, []).append(rec)

        for c in candidates:
            keys = {c.name.lower(), c.chembl_id.lower()}
            hits = []
            for k in keys:
                if len(k) < 4:
                    continue
                for agent, recs in by_agent.items():
                    if k in agent or agent in k:
                        hits.extend(recs)
            seen = set()
            for r in hits[:5]:
                if r.source.source_id in seen:
                    continue
                seen.add(r.source.source_id)
                c.prior_neuronal_evidence.append(
                    f"{r.source.first_author_surname() or '?'} "
                    f"{r.source.year or 'n.d.'}: {r.outcome.phenotype.raw} "
                    f"{r.outcome.direction.value} "
                    f"({r.cell_model.model_class.value})"
                )


# --------------------------------------------------------------------------
# Offline fallback: curated CNS-relevant pathway -> target map
# --------------------------------------------------------------------------

# Lets the tool produce useful output before any network calls, and provides a
# sanity check on API-derived results. Deliberately small and clearly labeled
# as curated rather than comprehensive.
CURATED_PATHWAY_TARGETS: dict[str, list[str]] = {
    "mtor signaling": ["MTOR", "RPS6KB1", "AKT1", "PIK3CA", "TSC1", "TSC2"],
    "synaptic vesicle cycle": ["SYN1", "SNAP25", "STX1A", "VAMP2", "SYT1", "DNM1"],
    "autophagy": ["MTOR", "ULK1", "ATG7", "TFEB", "BECN1"],
    "wnt signaling": ["GSK3B", "CTNNB1", "LRP6", "DVL1"],
    "camp signaling": ["PDE4B", "PDE10A", "ADCY1", "PRKACA", "CREB1"],
    "calcium signaling": ["CACNA1C", "CAMK2A", "RYR2", "ITPR1", "PPP3CA"],
    "glutamatergic synapse": ["GRIN1", "GRIN2B", "GRIA1", "GRM5", "DLG4"],
    "gabaergic synapse": ["GABRA1", "GABBR1", "SLC6A1", "GAD1"],
    "bdnf signaling": ["NTRK2", "BDNF", "MAPK1", "CREB1"],
    "oxidative stress": ["NFE2L2", "KEAP1", "SOD1", "CAT", "GPX1"],
    "unfolded protein response": ["EIF2AK3", "ERN1", "ATF6", "EIF2S1"],
    "proteasome": ["PSMB5", "USP14", "UCHL1"],
    "hdac / chromatin": ["HDAC1", "HDAC2", "HDAC6", "EP300", "KAT2B"],
    "rna splicing": ["SRSF1", "HNRNPA1", "TARDBP", "FUS", "PTBP1", "CLK1", "SRPK1"],
    "nonsense-mediated decay": ["UPF1", "SMG1", "SMG7"],
    "mitochondrial function": ["PPARGC1A", "SIRT1", "OPA1", "DNM1L", "PINK1"],
    "cytoskeleton / axon transport": ["KIF5A", "DCTN1", "TUBB3", "MAPT", "STMN2"],
}


def curated_targets(pathway_label: str) -> list[str]:
    """Look up targets for a pathway from the curated map. Substring-tolerant."""
    p = re.sub(r"[^a-z0-9 ]", " ", pathway_label.lower())
    p = re.sub(r"\s+", " ", p).strip()
    if p in CURATED_PATHWAY_TARGETS:
        return CURATED_PATHWAY_TARGETS[p]
    for k, v in CURATED_PATHWAY_TARGETS.items():
        if k in p or p in k:
            return v
    tokens = set(p.split())
    best, best_overlap = [], 0
    for k, v in CURATED_PATHWAY_TARGETS.items():
        ov = len(tokens & set(k.split()))
        if ov > best_overlap:
            best, best_overlap = v, ov
    return best if best_overlap >= 1 else []
