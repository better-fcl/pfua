#!/usr/bin/env bash
set -euo pipefail

source /your/path/to/conda/env
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

MMSEQS_TARGET_DB="/your/path/to/mmseqs_pipeline/02_target_db/targetDB"
UNIPROT_SPROT_DAT="/your/path/to/mmseqs_pipeline/02_target_db/uniprot_sprot.dat"
THREADS="16"
EVALUE="1e-5"
MIN_SEQ_ID="0.3"
COVERAGE="0.5"
NOSEQ="0"

AASEQ=""
QUERY_ID="query"
OUT_ROOT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --aaseq) AASEQ="$2"; shift 2 ;;
    --query-id) QUERY_ID="$2"; shift 2 ;;
    --out-root) OUT_ROOT="$2"; shift 2 ;;
    --target-db) MMSEQS_TARGET_DB="$2"; shift 2 ;;
    --uniprot-dat) UNIPROT_SPROT_DAT="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --evalue) EVALUE="$2"; shift 2 ;;
    --min-seq-id) MIN_SEQ_ID="$2"; shift 2 ;;
    --coverage) COVERAGE="$2"; shift 2 ;;
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
  python "${SCRIPT_DIR}/mmseqs2_besthit_uniprot.py" \
    --run-dir "${OUT_ROOT}" \
    --aaseq "${AASEQ}" \
    --query-id "${QUERY_ID}" \
    --target-db "${MMSEQS_TARGET_DB}" \
    --uniprot-dat "${UNIPROT_SPROT_DAT}" \
    --threads "${THREADS}" \
    --evalue "${EVALUE}" \
    --min-seq-id "${MIN_SEQ_ID}" \
    --coverage "${COVERAGE}" \
    --no-seq-in-evidence \
    > "${OUT_LOG}" 2> "${ERR_LOG}"
else
  python "${SCRIPT_DIR}/mmseqs2_besthit_uniprot.py" \
    --run-dir "${OUT_ROOT}" \
    --aaseq "${AASEQ}" \
    --query-id "${QUERY_ID}" \
    --target-db "${MMSEQS_TARGET_DB}" \
    --uniprot-dat "${UNIPROT_SPROT_DAT}" \
    --threads "${THREADS}" \
    --evalue "${EVALUE}" \
    --min-seq-id "${MIN_SEQ_ID}" \
    --coverage "${COVERAGE}" \
    > "${OUT_LOG}" 2> "${ERR_LOG}"
fi
RC=$?
set -e

if [[ $RC -ne 0 ]]; then
  cat "${ERR_LOG}" >&2
  exit $RC
fi

cat "${OUT_ROOT}/evidence.json"
