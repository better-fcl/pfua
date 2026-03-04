#!/usr/bin/env python3
import argparse
import json
import os
import re
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

AA_RE = re.compile(r"[^ACDEFGHIKLMNPQRSTVWY]")

def sanitize_aaseq(aaseq: str) -> Tuple[str, Dict[str, Any]]:
    s = aaseq.strip().upper()
    s = re.sub(r"\s+", "", s)
    bad = AA_RE.search(s)
    notes = {"removed_noncanonical": 0, "noncanonical_chars": []}
    if bad:
        raise ValueError(f"Invalid amino acid character: '{bad.group(0)}' in sequence.")
    if len(s) < 20:
        raise ValueError(f"Sequence too short (len={len(s)}).")
    return s, notes

def write_fasta(seq: str, out_fa: Path, header: str = "query") -> None:
    with out_fa.open("w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")

def _safe_float(x: str) -> Optional[float]:
    try:
        return float(x.replace("E", "e"))
    except Exception:
        return None

def _safe_int(x: str) -> Optional[int]:
    try:
        return int(float(x))
    except Exception:
        return None

def parse_domtblout(domtbl_path: Path, query_len: int, query_id: str) -> List[Dict[str, Any]]:
    hits: List[Dict[str, Any]] = []
    with domtbl_path.open("r", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 23:
                continue

            target_name = parts[0]
            target_acc = parts[1]
            qname = parts[3]
            i_evalue = _safe_float(parts[12])
            d_score = _safe_float(parts[13])
            hmm_from = _safe_int(parts[15])
            hmm_to = _safe_int(parts[16])
            ali_from = _safe_int(parts[17])
            ali_to = _safe_int(parts[18])
            desc = " ".join(parts[22:]).strip()

            if qname != query_id:
                continue
            if i_evalue is None or d_score is None or ali_from is None or ali_to is None:
                continue
            if ali_from <= 0 or ali_to <= 0 or ali_to < ali_from:
                continue

            aln_len = ali_to - ali_from + 1
            cov_q = aln_len / float(query_len) if query_len > 0 else 0.0

            hits.append({
                "pfam_id": target_name,
                "pfam_acc": target_acc,
                "query": qname,
                "evalue": i_evalue,
                "score": d_score,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "ali_from": ali_from,
                "ali_to": ali_to,
                "coverage_query": round(cov_q, 4),
                "desc": desc,
            })

    hits.sort(key=lambda h: (h["evalue"], -h["score"], -h["coverage_query"]))
    return hits

def select_domains_clean(hits: List[Dict[str, Any]], evalue_thr: float, min_cov: float) -> List[Dict[str, Any]]:
    strong = [h for h in hits if h["evalue"] <= evalue_thr and h["coverage_query"] >= min_cov]
    best_by_acc: Dict[str, Dict[str, Any]] = {}
    for h in strong:
        acc = h.get("pfam_acc") or h.get("pfam_id")
        if acc not in best_by_acc:
            best_by_acc[acc] = h
        else:
            cur = best_by_acc[acc]
            if (h["evalue"], -h["score"], -h["coverage_query"]) < (cur["evalue"], -cur["score"], -cur["coverage_query"]):
                best_by_acc[acc] = h
    selected = list(best_by_acc.values())
    selected.sort(key=lambda h: (h["evalue"], -h["score"], -h["coverage_query"]))
    return selected

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--aaseq", required=True)
    ap.add_argument("--query-id", default="query")

    ap.add_argument("--pfam-hmm", default=os.environ.get("PFAM_HMM", ""))
    ap.add_argument("--cpu", type=int, default=int(os.environ.get("HMMER_CPU", "8")))
    ap.add_argument("--incE", default=os.environ.get("HMMER_INCE", "1e-3"))
    ap.add_argument("--domE", default=os.environ.get("HMMER_DOME", "1e-3"))

    ap.add_argument("--select-evalue", type=float, default=float(os.environ.get("PFAM_SELECT_EVALUE", "1e-5")))
    ap.add_argument("--select-min-cov", type=float, default=float(os.environ.get("PFAM_SELECT_MIN_COV", "0.30")))
    ap.add_argument("--topk", type=int, default=10)
    ap.add_argument("--no-seq-in-evidence", action="store_true")
    args = ap.parse_args()

    run_dir = Path(args.run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    if not args.pfam_hmm or not Path(args.pfam_hmm).exists():
        raise RuntimeError("PFAM_HMM is not set or does not exist.")

    seq, _notes = sanitize_aaseq(args.aaseq)

    query_fa = run_dir / "query.fasta"
    domtbl = run_dir / "domtblout.txt"
    stdout_txt = run_dir / "hmmscan.stdout.txt"
    stderr_txt = run_dir / "hmmscan.stderr.txt"

    write_fasta(seq, query_fa, header=args.query_id)

    cmd = [
        "hmmscan",
        "--cpu", str(args.cpu),
        "--incE", str(args.incE),
        "--domE", str(args.domE),
        "--domtblout", str(domtbl),
        str(args.pfam_hmm),
        str(query_fa),
    ]

    with stdout_txt.open("w") as out_f:
        p = subprocess.run(cmd, cwd=str(run_dir), stdout=out_f, stderr=subprocess.PIPE, text=True)
    stderr_txt.write_text(p.stderr or "", encoding="utf-8", errors="ignore")

    if p.returncode != 0:
        raise RuntimeError(
            "hmmscan failed:\n"
            f"  cmd: {' '.join(cmd)}\n"
            f"  run_dir: {run_dir}\n"
            f"  stderr:\n{p.stderr}\n"
        )

    hits_all = parse_domtblout(domtbl, query_len=len(seq), query_id=args.query_id)
    hits = hits_all[: max(0, args.topk)]

    selected = select_domains_clean(hits_all, evalue_thr=args.select_evalue, min_cov=args.select_min_cov)
    selected_compact = selected[: max(0, args.topk)]

    evidence: Dict[str, Any] = {
        # "query": {"id": args.query_id, "len": len(seq), **({} if args.no_seq_in_evidence else {"sequence": seq})},
        # "hmmscan_params": {
        #     "pfam_hmm": str(Path(args.pfam_hmm)),
        #     "cpu": args.cpu,
        #     "incE": str(args.incE),
        #     "domE": str(args.domE),
        #     "select_evalue": args.select_evalue,
        #     "select_min_cov": args.select_min_cov,
        #     "topk": args.topk,
        # },
        "hits": hits,
        "selected_domains": selected_compact,
    }

    (run_dir / "evidence.json").write_text(json.dumps(evidence, ensure_ascii=False, indent=2), encoding="utf-8")

if __name__ == "__main__":
    main()
