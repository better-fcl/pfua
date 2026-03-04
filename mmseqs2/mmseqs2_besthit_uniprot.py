#!/usr/bin/env python3
import argparse
import hashlib
import json
import os
import re
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Any, Optional

AA_RE = re.compile(r"[^ACDEFGHIKLMNPQRSTVWY]")

def run(cmd: List[str], cwd: Optional[Path] = None) -> None:
    p = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if p.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            f"  cmd: {' '.join(cmd)}\n"
            f"  cwd: {cwd}\n"
            f"  stdout:\n{p.stdout}\n"
            f"  stderr:\n{p.stderr}\n"
        )

def sanitize_aaseq(aaseq: str) -> str:
    s = aaseq.strip().upper()
    s = re.sub(r"\s+", "", s)
    bad = AA_RE.search(s)
    if bad:
        raise ValueError(f"Invalid amino acid character: '{bad.group(0)}' in sequence.")
    # if len(s) < 20:
    if len(s) < 15:
        raise ValueError(f"Sequence too short (len={len(s)}).")
    return s

def write_fasta(seq: str, out_fa: Path, header: str = "query") -> None:
    with out_fa.open("w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")

def parse_convertalis_tsv(tsv_path: Path) -> List[Dict[str, Any]]:
    hits: List[Dict[str, Any]] = []
    with tsv_path.open() as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            q, t, pident, alnlen, evalue, bits = parts[:6]
            try:
                hits.append({
                    "query": q,
                    "target": t,
                    "pident": float(pident),
                    "alnlen": int(float(alnlen)),
                    "evalue": float(evalue.replace("E", "e")),
                    "bits": float(bits),
                })
            except Exception:
                continue
    return hits

def select_best_hit(hits: List[Dict[str, Any]]) -> Dict[str, Any]:
    hits_sorted = sorted(hits, key=lambda h: (h["evalue"], -h["bits"], -h["pident"], -h["alnlen"]))
    return hits_sorted[0]

def parse_uniprot_dat(dat_path: Path, wanted_acc: str) -> Dict[str, Any]:
    wanted = wanted_acc.strip()
    entry: Dict[str, Any] = {
        "accessions": [wanted],
        "protein_name": None,
        "function": [],
        "catalytic_activity": [],
        "ec": [],
        "cofactor": [],
        "subcellular_location": [],
        "go": [],
    }

    buf: List[str] = []
    in_entry = False

    def parse_entry(lines: List[str]) -> Dict[str, Any]:
        e = {k: ([] if isinstance(entry[k], list) else None) for k in entry.keys()}
        e["accessions"] = []
        current_cc_topic = None

        for ln in lines:
            if ln.startswith("AC   "):
                accs = [a.strip() for a in ln[5:].split(";") if a.strip()]
                e["accessions"].extend(accs)

            elif ln.startswith("DE   RecName: Full=") and e["protein_name"] is None:
                e["protein_name"] = ln.split("Full=", 1)[1].rstrip(";\n")

            elif ln.startswith("CC   -!- "):
                m = re.match(r"CC   -!- ([A-Z \-]+):\s*(.*)$", ln)
                if m:
                    current_cc_topic = m.group(1).strip()
                    payload = m.group(2).strip()
                else:
                    current_cc_topic = None
                    payload = ""

                if current_cc_topic == "FUNCTION" and payload:
                    e["function"].append(payload)
                elif current_cc_topic == "CATALYTIC ACTIVITY" and payload:
                    e["catalytic_activity"].append(payload)
                    for ec in re.findall(r"EC=([0-9]+\.[0-9]+\.[0-9]+\.[0-9]+)", payload):
                        e["ec"].append(ec)
                elif current_cc_topic == "COFACTOR" and payload:
                    e["cofactor"].append(payload)
                elif current_cc_topic == "SUBCELLULAR LOCATION" and payload:
                    e["subcellular_location"].append(payload)

            elif ln.startswith("CC       ") and current_cc_topic:
                payload = ln[9:].strip()
                if current_cc_topic == "FUNCTION" and payload:
                    e["function"].append(payload)
                elif current_cc_topic == "CATALYTIC ACTIVITY" and payload:
                    e["catalytic_activity"].append(payload)
                    for ec in re.findall(r"EC=([0-9]+\.[0-9]+\.[0-9]+\.[0-9]+)", payload):
                        e["ec"].append(ec)
                elif current_cc_topic == "COFACTOR" and payload:
                    e["cofactor"].append(payload)
                elif current_cc_topic == "SUBCELLULAR LOCATION" and payload:
                    e["subcellular_location"].append(payload)

            elif ln.startswith("DR   GO;"):
                e["go"].append(ln.strip())

        for k in ["accessions", "function", "catalytic_activity", "ec", "cofactor", "subcellular_location", "go"]:
            e[k] = list(dict.fromkeys(e[k]))
        return e

    def flush() -> Optional[Dict[str, Any]]:
        nonlocal buf
        if not buf:
            return None
        e = parse_entry(buf)
        if wanted in set(e.get("accessions", [])):
            return e
        return None

    with dat_path.open("r", errors="ignore") as f:
        for ln in f:
            if ln.startswith("ID   "):
                in_entry = True
                buf = [ln]
            elif in_entry:
                buf.append(ln)
                if ln.startswith("//"):
                    hit = flush()
                    if hit is not None:
                        return hit
                    buf = []
                    in_entry = False

    return entry

def sha256_hex(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--aaseq", required=True)
    ap.add_argument("--query-id", default="query")

    # fully exposed params
    ap.add_argument("--target-db", default=os.environ.get("MMSEQS_TARGET_DB", ""), help="mmseqs target DB prefix")
    ap.add_argument("--uniprot-dat", default=os.environ.get("UNIPROT_SPROT_DAT", ""), help="UniProt Swiss-Prot dat")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("MMSEQS_THREADS", "16")))
    ap.add_argument("--evalue", default=os.environ.get("MMSEQS_EVALUE", "1e-5"))
    ap.add_argument("--min-seq-id", default=os.environ.get("MMSEQS_MIN_SEQ_ID", "0.3"))
    ap.add_argument("--coverage", default=os.environ.get("MMSEQS_COVERAGE", "0.5"))
    ap.add_argument("--no-seq-in-evidence", action="store_true")
    args = ap.parse_args()

    run_dir = Path(args.run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    if not args.target_db or not Path(args.target_db).exists():
        raise RuntimeError("MMSEQS_TARGET_DB is not set or does not exist.")
    if not args.uniprot_dat or not Path(args.uniprot_dat).exists():
        raise RuntimeError("UNIPROT_SPROT_DAT is not set or does not exist.")

    seq = sanitize_aaseq(args.aaseq)

    fa = run_dir / "query.fasta"
    qdb = run_dir / "queryDB"
    result_db = run_dir / "resultDB"
    tmp_dir = run_dir / "tmp"
    tsv = run_dir / "result.tsv"
    evidence_path = run_dir / "evidence.json"

    write_fasta(seq, fa, header=args.query_id)
    tmp_dir.mkdir(exist_ok=True)

    run(["mmseqs", "createdb", str(fa), str(qdb)], cwd=run_dir)

    run([
        "mmseqs", "search",
        str(qdb), str(args.target_db), str(result_db), str(tmp_dir),
        "--threads", str(args.threads),
        "-e", str(args.evalue),
        "--min-seq-id", str(args.min_seq_id),
        "-c", str(args.coverage),
        "--cov-mode", "1",
        "-a"
    ], cwd=run_dir)

    run([
        "mmseqs", "convertalis",
        str(qdb), str(args.target_db), str(result_db), str(tsv),
        "--format-output", "query,target,pident,alnlen,evalue,bits"
    ], cwd=run_dir)

    hits = parse_convertalis_tsv(tsv)

    query_obj: Dict[str, Any] = {"id": args.query_id, "len": len(seq)}
    if args.no_seq_in_evidence:
        query_obj["sha256"] = sha256_hex(seq)
    else:
        query_obj["sequence"] = seq

    evidence: Dict[str, Any] = {
        # "query": query_obj,
        # "mmseqs_params": {
        #     "threads": args.threads,
        #     "evalue": str(args.evalue),
        #     "min_seq_id": str(args.min_seq_id),
        #     "coverage": str(args.coverage),
        # },
        "best_hit": None,
        "uniprot_annotation": None,
    }

    if hits:
        best = select_best_hit(hits)
        evidence["best_hit"] = best
        uniprot_id = best["target"]
        evidence["uniprot_annotation"] = parse_uniprot_dat(Path(args.uniprot_dat), uniprot_id)

    evidence_path.write_text(json.dumps(evidence, ensure_ascii=False, indent=2), encoding="utf-8")

if __name__ == "__main__":
    main()
