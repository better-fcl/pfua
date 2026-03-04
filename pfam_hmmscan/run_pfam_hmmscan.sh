#!/usr/bin/env bash
set -euo pipefail

source /your/path/to/conda/env
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Defaults (recommended)
PFAM_HMM="/your/path/to/mmseqs_pipeline/02_target_db/pfam/Pfam-A.hmm"
CPU="8"
INCE="1e-3"
DOME="1e-3"
SEL_E="1e-5"
SEL_COV="0.30"
TOPK="10"
NOSEQ="0"

AASEQ=""
QUERY_ID="query"
OUT_ROOT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --aaseq) AASEQ="$2"; shift 2 ;;
    --query-id) QUERY_ID="$2"; shift 2 ;;
    --out-root) OUT_ROOT="$2"; shift 2 ;;
    --pfam-hmm) PFAM_HMM="$2"; shift 2 ;;
    --cpu) CPU="$2"; shift 2 ;;
    --incE) INCE="$2"; shift 2 ;;
    --domE) DOME="$2"; shift 2 ;;
    --select-evalue) SEL_E="$2"; shift 2 ;;
    --select-min-cov) SEL_COV="$2"; shift 2 ;;
    --topk) TOPK="$2"; shift 2 ;;
    --no-seq-in-evidence) NOSEQ="1"; shift 1 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "${AASEQ}" ]]; then echo "ERROR: --aaseq is required" >&2; exit 2; fi
if [[ -z "${OUT_ROOT}" ]]; then echo "ERROR: --out-root is required" >&2; exit 2; fi

mkdir -p "${OUT_ROOT}"

OUT_LOG="${OUT_ROOT}/tool_stdout.log"
ERR_LOG="${OUT_ROOT}/tool_stderr.log"

set +e
if [[ "${NOSEQ}" == "1" ]]; then
  python "${SCRIPT_DIR}/pfam_hmmscan.py" \
    --run-dir "${OUT_ROOT}" \
    --aaseq "${AASEQ}" \
    --query-id "${QUERY_ID}" \
    --pfam-hmm "${PFAM_HMM}" \
    --cpu "${CPU}" \
    --incE "${INCE}" \
    --domE "${DOME}" \
    --select-evalue "${SEL_E}" \
    --select-min-cov "${SEL_COV}" \
    --topk "${TOPK}" \
    --no-seq-in-evidence \
    > "${OUT_LOG}" 2> "${ERR_LOG}"
else
  python "${SCRIPT_DIR}/pfam_hmmscan.py" \
    --run-dir "${OUT_ROOT}" \
    --aaseq "${AASEQ}" \
    --query-id "${QUERY_ID}" \
    --pfam-hmm "${PFAM_HMM}" \
    --cpu "${CPU}" \
    --incE "${INCE}" \
    --domE "${DOME}" \
    --select-evalue "${SEL_E}" \
    --select-min-cov "${SEL_COV}" \
    --topk "${TOPK}" \
    > "${OUT_LOG}" 2> "${ERR_LOG}"
fi
RC=$?
set -e

if [[ $RC -ne 0 ]]; then
  cat "${ERR_LOG}" >&2
  exit $RC
fi

cat "${OUT_ROOT}/evidence.json"
