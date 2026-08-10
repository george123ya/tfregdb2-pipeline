"""Apply ELM database short-linear-motif patterns to every curated TF sequence.

Replaces the ~7 hand-crafted SLiM regexes in `biophys.ts` (`predictSLiMs`)
with ELM's curated set of ~350 functional-site motifs. ELM's column
`Probability` is per-amino-acid match probability under a null model —
high-probability patterns produce noise, low-probability ones are
specific. We use `1 - log10(probability)` as a confidence weight and
keep only the top-scoring hits per TF.

Hits inside ordered Pfam domains (DBD / Oligo) are downweighted because
SLiMs functionally live in intrinsically-disordered regions; the
existing in-app heuristic already does this and we preserve the rule.

Output: `public/mock/slims/{uniprot}.json` shaped:
    [{ id: str, label: str, start: int (1-based), end: int, score: float }, ...]

The [symbol].astro page reads it and falls back to the in-app
`predictSLiMs` when missing — same pattern as pLDDT / secondary.
"""

from __future__ import annotations

import csv
import json
import math
import re
import sys
from pathlib import Path


ELM_TSV = Path("data/elm_classes.tsv")
TFS_DIR = Path("public/mock/tfs")
OUT_DIR = Path("public/mock/slims")
# Maximum hits per TF — keeps the per-TF JSON small + the SLIMTrack canvas
# readable. Sorted by score desc, so the most-specific motifs survive.
MAX_HITS_PER_TF = 80
# Minimum score (0..1 range) for a hit to be emitted.
MIN_SCORE = 0.30


def _load_elm() -> list[tuple[str, str, str, float]]:
    """Returns [(elm_id, label, regex_str, base_score), …].

    `base_score` is derived from ELM's Probability column: less probable =
    more specific = higher base score. We use 1 - log10(prob * 100) clipped
    to [0.3, 1.0] so the score range is intuitive in the UI.
    """
    rows: list[tuple[str, str, str, float]] = []
    with ELM_TSV.open() as fh:
        # Skip ELM's comment header (`#ELM_Classes_Download_Version: …`).
        lines = [ln for ln in fh if not ln.startswith("#")]
    reader = csv.DictReader(lines, dialect="excel-tab")
    for r in reader:
        elm_id = (r.get("ELMIdentifier") or "").strip().strip('"')
        regex_s = (r.get("Regex") or "").strip().strip('"')
        try:
            prob = float((r.get("Probability") or "").strip().strip('"'))
        except ValueError:
            prob = 1.0
        if not elm_id or not regex_s:
            continue
        # Probability log-score: ELM probabilities range ~1e-7 to ~1e-2.
        # log10(1e-7)=-7, log10(1e-2)=-2 → score 0.7..0.2.
        base = max(0.3, min(1.0, 1.0 + math.log10(max(prob, 1e-8) * 100)))
        # ELM IDs encode motif class as the prefix (CLV_, LIG_, MOD_, DEG_,
        # DOC_, TRG_). Use the first underscore-segment as the label so
        # the SLIMTrack legend stays readable.
        label = elm_id.split("_", 1)[0] if "_" in elm_id else elm_id
        rows.append((elm_id, label, regex_s, base))
    return rows


def _compile_elm(rows: list[tuple[str, str, str, float]]) -> list[tuple[str, str, re.Pattern, float]]:
    compiled: list[tuple[str, str, re.Pattern, float]] = []
    for elm_id, label, regex_s, base in rows:
        try:
            pat = re.compile(regex_s)
        except re.error:
            # A handful of ELM patterns use Perl-specific extensions
            # (lookarounds, named groups) Python's `re` may reject. Skip
            # rather than abort — the rest still apply.
            continue
        compiled.append((elm_id, label, pat, base))
    return compiled


def _build_ordered_mask(length: int, domains: list[dict]) -> bytearray:
    mask = bytearray(length)
    for d in domains:
        t = d.get("type", "")
        if t not in ("DNA", "Oligo"):
            continue
        s = max(0, int(d.get("start", 0)) - 1)
        e = min(length, int(d.get("end", 0)))
        for i in range(s, e):
            mask[i] = 1
    return mask


def extract_hits(sequence: str, domains: list[dict],
                 patterns: list[tuple[str, str, re.Pattern, float]]) -> list[dict]:
    if not sequence:
        return []
    ordered = _build_ordered_mask(len(sequence), domains)
    hits: list[dict] = []
    for elm_id, label, pat, base in patterns:
        for m in pat.finditer(sequence):
            start, end = m.start(), m.end()
            if start == end:
                continue
            span = end - start
            in_ord = sum(ordered[start:end])
            frac_ord = in_ord / span
            # SLiMs ARE biased toward IDRs (Tompa 2014). Penalize hits
            # buried in ordered DBD / Oligo regions by 60%; drop if the
            # adjusted score falls below MIN_SCORE.
            score = base * (1.0 - 0.6 * frac_ord)
            if score < MIN_SCORE:
                continue
            hits.append({
                "id": elm_id,
                "label": label,
                "start": start + 1,  # 1-based for the UI
                "end": end,
                "score": round(score, 3),
            })
    # Same-id overlaps collapse (e.g. greedy PEST runs).
    hits.sort(key=lambda h: (h["start"], h["end"]))
    merged: list[dict] = []
    for h in hits:
        prev = merged[-1] if merged else None
        if prev and prev["id"] == h["id"] and h["start"] <= prev["end"] + 1:
            prev["end"] = max(prev["end"], h["end"])
            prev["score"] = max(prev["score"], h["score"])
        else:
            merged.append(dict(h))
    # Keep the top-N by score so a 2,000-aa TF doesn't explode the JSON.
    merged.sort(key=lambda h: -h["score"])
    return merged[:MAX_HITS_PER_TF]


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    elm_rows = _load_elm()
    patterns = _compile_elm(elm_rows)
    print(f"[extract_slims_elm] loaded {len(elm_rows)} ELM classes · "
          f"{len(patterns)} compiled regex patterns")

    tf_files = sorted(TFS_DIR.glob("*.json"))
    n_ok = n_skip = total_hits = 0
    for i, fp in enumerate(tf_files, 1):
        tf = json.loads(fp.read_text())
        seq = tf.get("sequence") or ""
        upro = (tf.get("uniprot") or "").strip()
        if not seq or not upro:
            n_skip += 1
            continue
        hits = extract_hits(seq, tf.get("domains") or [], patterns)
        (OUT_DIR / f"{upro}.json").write_text(json.dumps(hits, separators=(",", ":")))
        n_ok += 1
        total_hits += len(hits)
        if i % 200 == 0:
            print(f"  [{i}/{len(tf_files)}] {n_ok} TFs · {total_hits:,} hits")
    print(f"  wrote {n_ok:,} SLiM files · {total_hits:,} total hits · "
          f"{total_hits/max(1, n_ok):.1f} avg/TF · → {OUT_DIR}")

    # ── Per-isoform pass: re-run ELM over each isoform's own sequence, keyed by
    # iso NAME (matches conservation / conway / adpred). The ordered-region
    # mask uses the isoform's own Pfam DBDs. Pure regex — fast, no extra data.
    n_iso = iso_hits = 0
    for fp in tf_files:
        tf = json.loads(fp.read_text())
        for iso in tf.get("isoforms") or []:
            name = iso.get("name")
            iseq = iso.get("sequence") or ""
            if not name or not iseq:
                continue
            iso_domains = [
                {"start": d.get("start"), "end": d.get("end"), "type": "DNA"}
                for d in (iso.get("pfamDBDs") or [])
            ]
            hits = extract_hits(iseq, iso_domains, patterns)
            (OUT_DIR / f"{name}.json").write_text(json.dumps(hits, separators=(",", ":")))
            n_iso += 1
            iso_hits += len(hits)
    print(f"  wrote {n_iso:,} isoform SLiM files · {iso_hits:,} hits")
    return 0


if __name__ == "__main__":
    sys.exit(main())
