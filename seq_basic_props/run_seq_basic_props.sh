#!/usr/bin/env bash
set -euo pipefail

source /your/path/to/conda/env
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

AASEQ=""
QUERY_ID="query"
OUT_ROOT=""
MEM_THR="18"
LC_THR="0.25"
INCLUDE_CLEAN="0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --aaseq) AASEQ="$2"; shift 2 ;;
    --query-id) QUERY_ID="$2"; shift 2 ;;
    --out-root) OUT_ROOT="$2"; shift 2 ;;
    --membrane-run-thr) MEM_THR="$2"; shift 2 ;;
    --low-complexity-thr) LC_THR="$2"; shift 2 ;;
    --include-cleaning-notes) INCLUDE_CLEAN="1"; shift 1 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "${AASEQ}" ]]; then echo "ERROR: --aaseq is required" >&2; exit 2; fi
if [[ -z "${OUT_ROOT}" ]]; then echo "ERROR: --out-root is required" >&2; exit 2; fi

mkdir -p "${OUT_ROOT}"

OUT_LOG="${OUT_ROOT}/tool_stdout.log"
ERR_LOG="${OUT_ROOT}/tool_stderr.log"

set +e
if [[ "${INCLUDE_CLEAN}" == "1" ]]; then
  python "${SCRIPT_DIR}/seq_basic_props.py" \
    --run-dir "${OUT_ROOT}" \
    --aaseq "${AASEQ}" \
    --membrane-run-thr "${MEM_THR}" \
    --low-complexity-thr "${LC_THR}" \
    --include-cleaning-notes \
    > "${OUT_LOG}" 2> "${ERR_LOG}"
else
  python "${SCRIPT_DIR}/seq_basic_props.py" \
    --run-dir "${OUT_ROOT}" \
    --aaseq "${AASEQ}" \
    --membrane-run-thr "${MEM_THR}" \
    --low-complexity-thr "${LC_THR}" \
    > "${OUT_LOG}" 2> "${ERR_LOG}"
fi
RC=$?
set -e

if [[ $RC -ne 0 ]]; then
  cat "${ERR_LOG}" >&2
  exit $RC
fi

cat "${OUT_ROOT}/evidence.json"
