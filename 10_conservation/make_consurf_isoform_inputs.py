"""Generate ConSurf input FASTAs for every TF isoform that needs scoring.

ConSurf grades depend only on the protein sequence, so the unit of work is a
*distinct* sequence. This script:

  - reads every isoform sequence from public/mock/tfs/*.json
  - skips isoforms whose sequence == their TF's canonical (they reuse the
    existing canonical /conservation/{uniprot}.json)
  - dedups the remaining (non-canonical) sequences globally
  - writes one FASTA per distinct sequence, named + headed exactly like the
    pipeline that made the canonical data:
        consurf_isoforms/inputs/{symbol}_{isoformName}_{uniprot}.fasta
        >{symbol}|{isoformName}|{uniprot}
  - writes consurf_isoforms/manifest.json mapping each FASTA basename to ALL
    isoform names sharing that sequence, plus the canonical-equal isoforms.

The FASTAs + the adapted array job (consurf_array_isoforms.sh) go to the SCC;
parse_consurf_isoforms.py turns the returned grades into per-isoform JSON.
"""

from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
TF_DIR = REPO / "public" / "mock" / "tfs"
OUT = REPO / "consurf_isoforms"
INPUTS = OUT / "inputs"


def main() -> None:
    INPUTS.mkdir(parents=True, exist_ok=True)
    for old in INPUTS.glob("*.fasta"):
        old.unlink()

    by_seq: dict[str, dict] = {}          # sequence -> run record
    canonical_equal: dict[str, str] = {}  # isoform name -> canonical uniprot

    for f in sorted(TF_DIR.glob("*.json")):
        d = json.loads(f.read_text())
        symbol = d["symbol"]
        canon_uni = d["uniprot"]
        isos = d.get("isoforms") or []
        canon = next((i for i in isos if i.get("canonical")), None)
        canon_seq = canon.get("sequence") if canon else None

        for iso in isos:
            seq = iso.get("sequence")
            name = iso.get("name")
            if not seq or not name:
                continue
            if seq == canon_seq:
                # Identical to canonical -> reuse the canonical ConSurf file.
                canonical_equal[name] = canon_uni
                continue
            rec = by_seq.get(seq)
            if rec is None:
                uni = iso.get("uniprotIso") or "NA"
                base = f"{symbol}_{name}_{uni}"
                by_seq[seq] = {"base": base, "symbol": symbol, "name": name,
                               "uniprot": uni, "names": [name]}
            else:
                rec["names"].append(name)

    # Write one FASTA per distinct non-canonical sequence.
    for seq, rec in by_seq.items():
        b = rec["base"]
        (INPUTS / f"{b}.fasta").write_text(
            f">{rec['symbol']}|{rec['name']}|{rec['uniprot']}\n{seq}\n")

    manifest = {
        "runs": {rec["base"]: rec["names"] for rec in by_seq.values()},
        "canonical_equal": canonical_equal,
        "n_runs": len(by_seq),
        "n_isoforms_covered": sum(len(r["names"]) for r in by_seq.values())
                              + len(canonical_equal),
        "n_canonical_equal": len(canonical_equal),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"[consurf-in] distinct sequences needing a run: {manifest['n_runs']}")
    print(f"[consurf-in] isoforms reusing canonical: {manifest['n_canonical_equal']}")
    print(f"[consurf-in] total isoforms covered: {manifest['n_isoforms_covered']}")
    print(f"[consurf-in] FASTAs -> {INPUTS}")
    print(f"[consurf-in] set the array job to:  #$ -t 1-{manifest['n_runs']}")


if __name__ == "__main__":
    main()
