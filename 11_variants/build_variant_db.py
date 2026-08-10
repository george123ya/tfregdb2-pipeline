"""Load the parsed ClinVar TSV (and later gnomAD / COSMIC TSVs) into a single
SQLite DB indexed on (chrom, pos) so per-TF projection can slice the
genomic ranges of each canonical ENSP's CDS exons in O(log n).

Schema is intentionally wide-and-flat — one row per variant per source. The
projector joins across sources by (chrom, pos, ref, alt) in code; we don't
need server-side joins because the read patterns are local range queries.
"""

from __future__ import annotations

import argparse
import csv
import sqlite3
import sys
from pathlib import Path


SCHEMA = """
CREATE TABLE IF NOT EXISTS variants (
    chrom         TEXT NOT NULL,
    pos           INTEGER NOT NULL,
    ref           TEXT NOT NULL,
    alt           TEXT NOT NULL,
    source        TEXT NOT NULL,        -- 'clinvar' | 'gnomad' | 'cosmic'
    clinvar_id    TEXT,
    cosmic_id     TEXT,
    rsid          TEXT,
    gene          TEXT,
    name          TEXT,                  -- HGVS or display name
    type          TEXT,                  -- SNV / Insertion / Deletion / etc.
    consequence   TEXT,                  -- missense / nonsense / synonymous / ...
    significance  TEXT,                  -- ClinVar ClinicalSignificance
    review_status TEXT,
    origin        TEXT,
    af            REAL,                  -- gnomAD allele frequency
    ac            INTEGER,
    an            INTEGER,
    popmax_af     REAL,                  -- gnomAD popmax AF across ancestries
    pmids         TEXT,
    -- COSMIC per-transcript fields. COSMIC v103 ships one row per
    -- (variant × transcript), so the same genomic SNV produces multiple
    -- rows here, each carrying its own transcript-specific HGVS.p and
    -- SO-term consequence. The projector matches `transcript_id`
    -- against the TF's per-iso ENSTs to attach iso-correct annotations.
    -- Empty / NULL for ClinVar + gnomAD rows where this granularity
    -- isn't available.
    transcript_id TEXT,                  -- Ensembl ENST (unversioned)
    hgvs_p        TEXT,
    so_term       TEXT,
    is_canonical  INTEGER                -- 1 if COSMIC marked this row canonical
);
CREATE INDEX IF NOT EXISTS idx_variants_pos ON variants(chrom, pos);
CREATE INDEX IF NOT EXISTS idx_variants_gene ON variants(gene);
CREATE INDEX IF NOT EXISTS idx_variants_enst ON variants(transcript_id);
"""


def _bulk_insert_clinvar(con: sqlite3.Connection, tsv_path: Path) -> int:
    cur = con.cursor()
    n = 0
    with tsv_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        batch: list[tuple] = []
        for row in reader:
            batch.append((
                row["chrom"],
                int(row["pos"]),
                row["ref"],
                row["alt"],
                "clinvar",
                row.get("clinvar_id") or None,
                None,
                row.get("rsid") or None,
                row.get("gene") or None,
                row.get("name") or None,
                row.get("type") or None,
                row.get("consequence") or None,
                row.get("significance") or None,
                row.get("review_status") or None,
                row.get("origin") or None,
                None, None, None, None,
                row.get("pmids") or None,
                # transcript_id, hgvs_p, so_term, is_canonical (NULL — ClinVar
                # variant_summary has no per-transcript granularity beyond
                # one preferred-transcript HGVS in `name`).
                None, None, None, None,
            ))
            if len(batch) >= 10000:
                cur.executemany(
                    "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    batch,
                )
                n += len(batch)
                batch.clear()
                if n % 200000 == 0:
                    print(f"  inserted {n:,}")
        if batch:
            cur.executemany(
                "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                batch,
            )
            n += len(batch)
    con.commit()
    return n


def _bulk_insert_gnomad(con: sqlite3.Connection, tsv_path: Path) -> int:
    cur = con.cursor()
    n = 0
    with tsv_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        batch: list[tuple] = []
        for row in reader:
            try:
                af = float(row["af"]) if row.get("af") else None
                ac = int(row["ac"]) if row.get("ac") else None
                an = int(row["an"]) if row.get("an") else None
                popmax = float(row["popmax_af"]) if row.get("popmax_af") else None
            except ValueError:
                af = ac = an = popmax = None
            batch.append((
                row["chrom"],
                int(row["pos"]),
                row["ref"],
                row["alt"],
                "gnomad",
                None, None, None, None, None, None, None, None, None, None,
                af, ac, an, popmax, None,
                # transcript_id, hgvs_p, so_term, is_canonical (NULL — we
                # don't currently keep gnomAD CSQ on ingest. Step 3 of the
                # VEP-CSQ ramp would fill these in.)
                None, None, None, None,
            ))
            if len(batch) >= 10000:
                cur.executemany(
                    "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    batch,
                )
                n += len(batch)
                batch.clear()
                if n % 500000 == 0:
                    print(f"  inserted {n:,}")
        if batch:
            cur.executemany(
                "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                batch,
            )
            n += len(batch)
    con.commit()
    return n


def _bulk_insert_cosmic(con: sqlite3.Connection, tsv_path: Path) -> int:
    cur = con.cursor()
    n = 0
    with tsv_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        batch: list[tuple] = []
        for row in reader:
            try:
                count = int(row["count"]) if row.get("count") else 1
            except ValueError:
                count = 1
            try:
                is_canon = int(row.get("is_canonical") or "0")
            except ValueError:
                is_canon = 0
            batch.append((
                row["chrom"],
                int(row["pos"]),
                row["ref"],
                row["alt"],
                "cosmic",
                None,
                row.get("cosmic_id") or None,
                None,
                row.get("gene") or None,
                row.get("hgvs_p") or row.get("consequence") or None,  # name = HGVS.p when available
                row.get("consequence") or None,
                row.get("so_term") or row.get("consequence") or None,
                None, None, None,
                None,
                count,  # store COSMIC sample count in `ac` so the canvas log-scale gets a real signal
                None, None, None,
                row.get("transcript") or None,
                row.get("hgvs_p") or None,
                row.get("so_term") or None,
                is_canon,
            ))
            if len(batch) >= 10000:
                cur.executemany(
                    "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    batch,
                )
                n += len(batch)
                batch.clear()
        if batch:
            cur.executemany(
                "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                batch,
            )
            n += len(batch)
    con.commit()
    return n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clinvar", type=Path, default=Path("data/variants/clinvar.tsv"))
    ap.add_argument("--gnomad", type=Path, default=Path("data/variants/gnomad.tsv"))
    ap.add_argument("--cosmic", type=Path, default=Path("data/variants/cosmic.tsv"))
    ap.add_argument("--out", type=Path, default=Path("data/variants/variants.sqlite"))
    ap.add_argument("--reset", action="store_true",
                    help="Drop and re-create the variants table before ingest.")
    args = ap.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(args.out)
    cur = con.cursor()

    if args.reset:
        print("[build_variant_db] resetting table")
        cur.executescript("DROP TABLE IF EXISTS variants;")

    cur.executescript(SCHEMA)
    con.commit()

    for src_name, src_path, inserter in (
        ("clinvar", args.clinvar, _bulk_insert_clinvar),
        ("gnomad", args.gnomad, _bulk_insert_gnomad),
        ("cosmic", args.cosmic, _bulk_insert_cosmic),
    ):
        if not src_path.exists():
            print(f"[build_variant_db] no {src_name} TSV at {src_path} — skip")
            continue
        print(f"[build_variant_db] ingesting {src_path}")
        cur.execute("DELETE FROM variants WHERE source = ?", (src_name,))
        con.commit()
        n = inserter(con, src_path)
        print(f"  inserted {n:,} {src_name} rows")

    cur.execute("SELECT COUNT(*) FROM variants")
    total = cur.fetchone()[0]
    print(f"\n  total rows in DB: {total:,}")
    print(f"  → {args.out} ({args.out.stat().st_size / 1e6:.1f} MB)")
    con.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
