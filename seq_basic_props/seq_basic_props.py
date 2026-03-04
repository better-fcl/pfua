#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Dict, Any, Tuple

AA_CANONICAL = set("ACDEFGHIKLMNPQRSTVWY")
HYDROPHOBIC = set("AILMFWV")

def clean_sequence(raw: str) -> Tuple[str, Dict[str, Any]]:
    s = re.sub(r"\s+", "", raw.strip().upper())
    cleaned = []
    removed = 0
    bad = Counter()
    for ch in s:
        if ch in AA_CANONICAL:
            cleaned.append(ch)
        else:
            removed += 1
            bad[ch] += 1
    notes = {
        "removed_noncanonical": removed,
        "noncanonical_chars": [{"char": k, "count": v} for k, v in bad.most_common()],
    }
    return "".join(cleaned), notes

def longest_hydrophobic_run(seq: str) -> int:
    best = 0
    cur = 0
    for ch in seq:
        if ch in HYDROPHOBIC:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best

def shannon_entropy(seq: str) -> float:
    n = len(seq)
    if n == 0:
        return 0.0
    c = Counter(seq)
    ent = 0.0
    for v in c.values():
        p = v / n
        ent -= p * math.log2(p)
    return ent

def low_complexity_index(seq: str) -> float:
    ent = shannon_entropy(seq)
    ent_norm = min(max(ent / 4.32, 0.0), 1.0)
    return round(1.0 - ent_norm, 4)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True, help="Per-call output directory (created if not exists).")
    ap.add_argument("--aaseq", required=True, help="Amino acid sequence.")
    ap.add_argument("--membrane-run-thr", type=int, default=18, help="Hydrophobic-run threshold for membrane-like heuristic.")
    ap.add_argument("--low-complexity-thr", type=float, default=0.25, help="Low-complexity index threshold heuristic.")
    ap.add_argument("--include-cleaning-notes", action="store_true", help="Include cleaning notes in evidence.")
    args = ap.parse_args()

    run_dir = Path(args.run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    seq, clean_notes = clean_sequence(args.aaseq)
    length = len(seq)
    hydro_run = longest_hydrophobic_run(seq)
    lcx = low_complexity_index(seq)

    evidence: Dict[str, Any] = {
        "length": length,
        "hydrophobic_run_max": hydro_run,
        "low_complexity_index_0to1": lcx,
        "heuristics": {
            "looks_membrane_like": hydro_run >= int(args.membrane_run_thr),
            "looks_low_complexity_like": lcx >= float(args.low_complexity_thr),
        },
    }
    if args.include_cleaning_notes:
        evidence["cleaning_notes"] = clean_notes

    (run_dir / "evidence.json").write_text(json.dumps(evidence, ensure_ascii=False, indent=2), encoding="utf-8")

if __name__ == "__main__":
    main()
