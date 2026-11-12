"""
Persistence layer. SQLite because it needs zero administration on a lab
server, handles tens of millions of rows fine for this access pattern, and
the whole knowledge base stays a single file you can scp to a collaborator.

Records are stored denormalized-as-JSON plus extracted columns for the fields
we filter on. This keeps schema migration cheap: adding a field to
ExperimentRecord does not require an ALTER TABLE unless we want to index it.
"""

from __future__ import annotations

import json
import sqlite3
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator

from .schema import ExperimentRecord, SourceMetadata

SCHEMA_VERSION = 3

DDL = """
PRAGMA journal_mode=WAL;
PRAGMA synchronous=NORMAL;
PRAGMA foreign_keys=ON;

CREATE TABLE IF NOT EXISTS meta (
    key TEXT PRIMARY KEY,
    value TEXT
);

CREATE TABLE IF NOT EXISTS sources (
    source_id       TEXT PRIMARY KEY,
    pmid            TEXT,
    pmcid           TEXT,
    doi             TEXT,
    title           TEXT,
    journal         TEXT,
    year            INTEGER,
    is_preprint     INTEGER DEFAULT 0,
    is_review       INTEGER DEFAULT 0,
    is_retracted    INTEGER DEFAULT 0,
    lab_id          TEXT,
    citation_count  INTEGER,
    sjr_quartile    TEXT,
    sjr_score       REAL,
    open_access     INTEGER DEFAULT 0,
    full_text_available INTEGER DEFAULT 0,
    payload         TEXT NOT NULL,
    fetched_at      TEXT
);
CREATE INDEX IF NOT EXISTS idx_sources_lab ON sources(lab_id);
CREATE INDEX IF NOT EXISTS idx_sources_year ON sources(year);
CREATE INDEX IF NOT EXISTS idx_sources_journal ON sources(journal);

CREATE TABLE IF NOT EXISTS documents (
    source_id   TEXT PRIMARY KEY REFERENCES sources(source_id) ON DELETE CASCADE,
    abstract    TEXT,
    full_text   TEXT,
    methods     TEXT,
    results     TEXT,
    sections_json TEXT
);

CREATE TABLE IF NOT EXISTS records (
    record_id       TEXT PRIMARY KEY,
    source_id       TEXT NOT NULL REFERENCES sources(source_id) ON DELETE CASCADE,
    model_class     TEXT,
    maturity_bin    TEXT,
    days_in_vitro   INTEGER,
    perturbation_type TEXT,
    perturbation_target_key TEXT,
    variant_hgvs    TEXT,
    assay_class     TEXT,
    ortho_group     TEXT,
    phenotype_key   TEXT,
    phenotype_raw   TEXT,
    direction       TEXT,
    effect_size     REAL,
    effect_unit     TEXT,
    p_value         REAL,
    evidence_type   TEXT,
    is_negative     INTEGER DEFAULT 0,
    is_quantified   INTEGER DEFAULT 0,
    lab_id          TEXT,
    year            INTEGER,
    payload         TEXT NOT NULL,
    ingested_at     TEXT
);
CREATE INDEX IF NOT EXISTS idx_rec_phen ON records(phenotype_key);
CREATE INDEX IF NOT EXISTS idx_rec_src ON records(source_id);
CREATE INDEX IF NOT EXISTS idx_rec_assay ON records(assay_class);
CREATE INDEX IF NOT EXISTS idx_rec_model ON records(model_class, maturity_bin);
CREATE INDEX IF NOT EXISTS idx_rec_target ON records(perturbation_target_key);
CREATE INDEX IF NOT EXISTS idx_rec_lab ON records(lab_id);
CREATE INDEX IF NOT EXISTS idx_rec_variant ON records(variant_hgvs);

-- record -> pathway / gene edges, exploded for join performance
CREATE TABLE IF NOT EXISTS record_pathways (
    record_id   TEXT NOT NULL REFERENCES records(record_id) ON DELETE CASCADE,
    pathway_key TEXT NOT NULL,
    pathway_label TEXT,
    PRIMARY KEY (record_id, pathway_key)
);
CREATE INDEX IF NOT EXISTS idx_rp_pathway ON record_pathways(pathway_key);

CREATE TABLE IF NOT EXISTS record_genes (
    record_id TEXT NOT NULL REFERENCES records(record_id) ON DELETE CASCADE,
    gene_key  TEXT NOT NULL,
    gene_label TEXT,
    PRIMARY KEY (record_id, gene_key)
);
CREATE INDEX IF NOT EXISTS idx_rg_gene ON record_genes(gene_key);

CREATE TABLE IF NOT EXISTS splicing_events (
    record_id   TEXT PRIMARY KEY REFERENCES records(record_id) ON DELETE CASCADE,
    event_type  TEXT,
    gene_key    TEXT,
    coordinates TEXT,
    genome_build TEXT,
    delta_psi   REAL,
    detection_method TEXT,
    regulator_key TEXT
);
CREATE INDEX IF NOT EXISTS idx_spl_gene ON splicing_events(gene_key);
CREATE INDEX IF NOT EXISTS idx_spl_reg ON splicing_events(regulator_key);

-- Full text search over source titles/abstracts for candidate retrieval
CREATE VIRTUAL TABLE IF NOT EXISTS source_fts USING fts5(
    source_id UNINDEXED,
    title,
    abstract,
    tokenize='porter unicode61'
);

-- Cache of external API calls so re-runs are cheap and reproducible
CREATE TABLE IF NOT EXISTS api_cache (
    cache_key   TEXT PRIMARY KEY,
    payload     TEXT,
    fetched_at  TEXT
);
"""


class Store:
    """Thin wrapper over SQLite. Not thread-safe; open one per worker."""

    def __init__(self, path: str | Path = "data/synaptik.db"):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(str(self.path), timeout=60.0)
        self.conn.row_factory = sqlite3.Row
        self._init_schema()

    def _init_schema(self) -> None:
        self.conn.executescript(DDL)
        cur = self.conn.execute("SELECT value FROM meta WHERE key='schema_version'")
        row = cur.fetchone()
        if row is None:
            self.conn.execute(
                "INSERT INTO meta(key,value) VALUES('schema_version',?)",
                (str(SCHEMA_VERSION),),
            )
        self.conn.commit()

    @contextmanager
    def tx(self) -> Iterator[sqlite3.Connection]:
        try:
            yield self.conn
            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise

    # ---------------- sources ----------------

    def upsert_source(self, src: SourceMetadata, fetched_at: str | None = None) -> None:
        self.conn.execute(
            """INSERT INTO sources
               (source_id,pmid,pmcid,doi,title,journal,year,is_preprint,is_review,
                is_retracted,lab_id,citation_count,sjr_quartile,sjr_score,
                open_access,full_text_available,payload,fetched_at)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
               ON CONFLICT(source_id) DO UPDATE SET
                 title=excluded.title, journal=excluded.journal, year=excluded.year,
                 is_retracted=excluded.is_retracted, lab_id=excluded.lab_id,
                 citation_count=excluded.citation_count,
                 sjr_quartile=excluded.sjr_quartile, sjr_score=excluded.sjr_score,
                 open_access=excluded.open_access,
                 full_text_available=excluded.full_text_available,
                 payload=excluded.payload, fetched_at=excluded.fetched_at""",
            (
                src.source_id, src.pmid, src.pmcid, src.doi, src.title, src.journal,
                src.year, int(src.is_preprint), int(src.is_review),
                int(src.is_retracted), src.lab_id, src.citation_count,
                src.sjr_quartile, src.sjr_score, int(src.open_access),
                int(src.full_text_available), src.model_dump_json(),
                fetched_at or datetime.now(timezone.utc).isoformat(),
            ),
        )

    def index_source_text(self, source_id: str, title: str, abstract: str) -> None:
        self.conn.execute("DELETE FROM source_fts WHERE source_id=?", (source_id,))
        self.conn.execute(
            "INSERT INTO source_fts(source_id,title,abstract) VALUES (?,?,?)",
            (source_id, title or "", abstract or ""),
        )

    def store_document(
        self,
        source_id: str,
        abstract: str = "",
        full_text: str = "",
        methods: str = "",
        results: str = "",
        sections: dict | None = None,
    ) -> None:
        self.conn.execute(
            """INSERT INTO documents(source_id,abstract,full_text,methods,results,sections_json)
               VALUES (?,?,?,?,?,?)
               ON CONFLICT(source_id) DO UPDATE SET
                 abstract=excluded.abstract, full_text=excluded.full_text,
                 methods=excluded.methods, results=excluded.results,
                 sections_json=excluded.sections_json""",
            (source_id, abstract, full_text, methods, results,
             json.dumps(sections or {})),
        )

    def get_source(self, source_id: str) -> SourceMetadata | None:
        row = self.conn.execute(
            "SELECT payload FROM sources WHERE source_id=?", (source_id,)
        ).fetchone()
        return SourceMetadata.model_validate_json(row["payload"]) if row else None

    def get_document(self, source_id: str) -> dict | None:
        row = self.conn.execute(
            "SELECT * FROM documents WHERE source_id=?", (source_id,)
        ).fetchone()
        return dict(row) if row else None

    def has_source(self, source_id: str) -> bool:
        return (
            self.conn.execute(
                "SELECT 1 FROM sources WHERE source_id=?", (source_id,)
            ).fetchone()
            is not None
        )

    # ---------------- records ----------------

    def insert_record(self, rec: ExperimentRecord) -> None:
        self.conn.execute(
            """INSERT INTO records
               (record_id,source_id,model_class,maturity_bin,days_in_vitro,
                perturbation_type,perturbation_target_key,variant_hgvs,assay_class,
                ortho_group,phenotype_key,phenotype_raw,direction,effect_size,
                effect_unit,p_value,evidence_type,is_negative,is_quantified,
                lab_id,year,payload,ingested_at)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
               ON CONFLICT(record_id) DO UPDATE SET payload=excluded.payload""",
            (
                rec.record_id, rec.source.source_id, rec.cell_model.model_class.value,
                rec.cell_model.maturity_bin(), rec.cell_model.days_in_vitro,
                rec.perturbation.ptype.value,
                rec.perturbation.target.key() if rec.perturbation.target else None,
                rec.perturbation.variant_hgvs, rec.assay_class.value,
                rec.orthogonality_group(), rec.outcome.phenotype.key(),
                rec.outcome.phenotype.raw, rec.outcome.direction.value,
                rec.outcome.effect_size, rec.outcome.effect_size_unit,
                rec.outcome.p_value, rec.evidence_type.value,
                int(rec.is_negative_result()), int(rec.outcome.is_quantified),
                rec.source.lab_id, rec.source.year, rec.model_dump_json(),
                rec.ingested_at or datetime.now(timezone.utc).isoformat(),
            ),
        )
        for p in rec.pathways:
            self.conn.execute(
                """INSERT OR IGNORE INTO record_pathways(record_id,pathway_key,pathway_label)
                   VALUES (?,?,?)""",
                (rec.record_id, p.key(), p.ontology_label or p.raw),
            )
        for g in rec.genes:
            self.conn.execute(
                """INSERT OR IGNORE INTO record_genes(record_id,gene_key,gene_label)
                   VALUES (?,?,?)""",
                (rec.record_id, g.key(), g.ontology_label or g.raw),
            )
        if rec.splicing:
            s = rec.splicing
            self.conn.execute(
                """INSERT OR REPLACE INTO splicing_events
                   (record_id,event_type,gene_key,coordinates,genome_build,
                    delta_psi,detection_method,regulator_key)
                   VALUES (?,?,?,?,?,?,?,?)""",
                (
                    rec.record_id, s.event_type.value,
                    s.gene.key() if s.gene else None, s.coordinates, s.genome_build,
                    s.delta_psi, s.detection_method,
                    s.regulator.key() if s.regulator else None,
                ),
            )

    def insert_records(self, recs: Iterable[ExperimentRecord]) -> int:
        n = 0
        with self.tx():
            for r in recs:
                self.insert_record(r)
                n += 1
        return n

    def get_record(self, record_id: str) -> ExperimentRecord | None:
        row = self.conn.execute(
            "SELECT payload FROM records WHERE record_id=?", (record_id,)
        ).fetchone()
        return ExperimentRecord.model_validate_json(row["payload"]) if row else None

    def iter_records(self, where: str = "", params: tuple = ()) -> Iterator[ExperimentRecord]:
        sql = "SELECT payload FROM records"
        if where:
            sql += " WHERE " + where
        for row in self.conn.execute(sql, params):
            yield ExperimentRecord.model_validate_json(row["payload"])

    def records_for_phenotype(self, phenotype_keys: list[str]) -> list[ExperimentRecord]:
        if not phenotype_keys:
            return []
        q = ",".join("?" * len(phenotype_keys))
        rows = self.conn.execute(
            f"SELECT payload FROM records WHERE phenotype_key IN ({q})",
            tuple(phenotype_keys),
        ).fetchall()
        return [ExperimentRecord.model_validate_json(r["payload"]) for r in rows]

    def search_sources(self, query: str, limit: int = 100) -> list[str]:
        """FTS5 candidate retrieval. Returns source_ids ranked by bm25."""
        try:
            rows = self.conn.execute(
                """SELECT source_id FROM source_fts WHERE source_fts MATCH ?
                   ORDER BY bm25(source_fts) LIMIT ?""",
                (query, limit),
            ).fetchall()
        except sqlite3.OperationalError:
            return []
        return [r["source_id"] for r in rows]

    # ---------------- cache ----------------

    def cache_get(self, key: str) -> dict | None:
        row = self.conn.execute(
            "SELECT payload FROM api_cache WHERE cache_key=?", (key,)
        ).fetchone()
        return json.loads(row["payload"]) if row else None

    def cache_put(self, key: str, payload: dict) -> None:
        self.conn.execute(
            """INSERT OR REPLACE INTO api_cache(cache_key,payload,fetched_at)
               VALUES (?,?,?)""",
            (key, json.dumps(payload), datetime.now(timezone.utc).isoformat()),
        )
        self.conn.commit()

    # ---------------- stats ----------------

    def stats(self) -> dict:
        def one(sql: str) -> int:
            return self.conn.execute(sql).fetchone()[0]

        return {
            "sources": one("SELECT COUNT(*) FROM sources"),
            "documents": one("SELECT COUNT(*) FROM documents"),
            "records": one("SELECT COUNT(*) FROM records"),
            "negative_records": one("SELECT COUNT(*) FROM records WHERE is_negative=1"),
            "quantified_records": one("SELECT COUNT(*) FROM records WHERE is_quantified=1"),
            "splicing_events": one("SELECT COUNT(*) FROM splicing_events"),
            "distinct_labs": one("SELECT COUNT(DISTINCT lab_id) FROM sources WHERE lab_id IS NOT NULL"),
            "distinct_phenotypes": one("SELECT COUNT(DISTINCT phenotype_key) FROM records"),
            "retracted": one("SELECT COUNT(*) FROM sources WHERE is_retracted=1"),
        }

    def close(self) -> None:
        self.conn.close()
