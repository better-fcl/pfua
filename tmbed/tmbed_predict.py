#!/usr/bin/env python3
import argparse
import json
import os
import re
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple

AA_RE = re.compile(r"[^ACDEFGHIKLMNPQRSTVWY]")

def sanitize_aaseq(aaseq: str) -> str:
    s = aaseq.strip().upper()
    s = re.sub(r"\s+", "", s)
    bad = AA_RE.search(s)
    if bad:
        raise ValueError(f"Invalid amino acid character: '{bad.group(0)}' in sequence.")
    if len(s) < 10:
        raise ValueError(f"Sequence too short (len={len(s)}).")
    return s

def write_fasta(seq: str, out_fa: Path, header: str) -> None:
    with out_fa.open("w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")

def run_cmd(cmd: List[str], cwd: Optional[Path] = None, env: Optional[Dict[str, str]] = None) -> Tuple[str, str]:
    p = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
    )
    if p.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            f"  cmd: {' '.join(cmd)}\n"
            f"  cwd: {cwd}\n"
            f"  stdout:\n{p.stdout}\n"
            f"  stderr:\n{p.stderr}\n"
        )
    return p.stdout, p.stderr

def parse_tmbed_pred_minimal(pred_text: str) -> Dict[str, Any]:
    raw = pred_text.strip()
    lines = [ln.strip() for ln in raw.splitlines() if ln.strip()]
    tm_letters = set("HMB")
    tm_hits = 0
    for ln in lines:
        tm_hits += sum(1 for ch in ln if ch in tm_letters)
    return {
        "raw_pred": raw,
        "tm_signal_letter_hits": tm_hits,
        "has_tm_signal_heuristic": (tm_hits > 10),
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--aaseq", required=True)
    ap.add_argument("--query-id", default="query")

    ap.add_argument("--model-dir", default=os.environ.get("TMBED_MODEL_DIR", "").strip())
    ap.add_argument("--use-gpu", type=int, default=1 if os.environ.get("TMBED_USE_GPU", "1").strip() == "1" else 0)
    ap.add_argument("--threads", type=int, default=int(os.environ.get("TMBED_THREADS", "8")))
    ap.add_argument("--batch-size", type=int, default=int(os.environ.get("TMBED_BATCH_SIZE", "2000")))
    ap.add_argument("--cuda-visible-devices", default=os.environ.get("CUDA_VISIBLE_DEVICES", "").strip())
    ap.add_argument("--no-seq-in-evidence", action="store_true")
    args = ap.parse_args()

    seq = sanitize_aaseq(args.aaseq)
    if not args.model_dir:
        raise RuntimeError("TMBED_MODEL_DIR is not set (pass --model-dir or set env).")

    run_dir = Path(args.run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    fa = run_dir / "query.fasta"
    pred_path = run_dir / "tmbed.pred"
    stdout_path = run_dir / "tmbed.stdout.txt"
    stderr_path = run_dir / "tmbed.stderr.txt"

    write_fasta(seq, fa, header=args.query_id)

    cmd = [
        "tmbed", "predict",
        "-f", str(fa),
        "-p", str(pred_path),
        "--model-dir", str(args.model_dir),
        "--batch-size", str(args.batch_size),
    ]
    if int(args.use_gpu) == 1:
        cmd += ["--use-gpu"]
    else:
        cmd += ["--no-use-gpu", "--threads", str(args.threads)]

    env = os.environ.copy()
    if args.cuda_visible_devices:
        env["CUDA_VISIBLE_DEVICES"] = args.cuda_visible_devices

    out, err = run_cmd(cmd, cwd=run_dir, env=env)
    stdout_path.write_text(out or "", encoding="utf-8", errors="ignore")
    stderr_path.write_text(err or "", encoding="utf-8", errors="ignore")

    pred_text = pred_path.read_text(encoding="utf-8", errors="ignore")
    parsed = parse_tmbed_pred_minimal(pred_text)

    evidence: Dict[str, Any] = {
        # "query": {"id": args.query_id, "len": len(seq), **({} if args.no_seq_in_evidence else {"sequence": seq})},
        # "tmbed_params": {
        #     "model_dir": str(args.model_dir),
        #     "use_gpu": int(args.use_gpu),
        #     "threads": int(args.threads),
        #     "batch_size": int(args.batch_size),
        #     "cuda_visible_devices": args.cuda_visible_devices,
        # },
        "prediction": parsed,
    }

    (run_dir / "evidence.json").write_text(json.dumps(evidence, ensure_ascii=False, indent=2), encoding="utf-8")

if __name__ == "__main__":
    main()
