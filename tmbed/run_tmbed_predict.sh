#!/usr/bin/env bash
set -euo pipefail

source /your/path/to/conda/env
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Defaults (recommended)
TMBED_MODEL_DIR="/your/path/to/tmbed/models"
USE_GPU="1"
THREADS="8"
BATCH_SIZE="2000"
CUDA_VISIBLE_DEVICES_DEFAULT="5"
NOSEQ="0"

AASEQ=""
QUERY_ID="query"
OUT_ROOT=""
MODEL_DIR=""
CUDA_VISIBLE_DEVICES_ARG=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --aaseq) AASEQ="$2"; shift 2 ;;
    --query-id) QUERY_ID="$2"; shift 2 ;;
    --out-root) OUT_ROOT="$2"; shift 2 ;;
    --model-dir) MODEL_DIR="$2"; shift 2 ;;
    --use-gpu) USE_GPU="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --batch-size) BATCH_SIZE="$2"; shift 2 ;;
    --cuda-visible-devices) CUDA_VISIBLE_DEVICES_ARG="$2"; shift 2 ;;
    --no-seq-in-evidence) NOSEQ="1"; shift 1 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "${AASEQ}" ]]; then echo "ERROR: --aaseq is required" >&2; exit 2; fi
if [[ -z "${OUT_ROOT}" ]]; then echo "ERROR: --out-root is required" >&2; exit 2; fi

mkdir -p "${OUT_ROOT}"
OUT_LOG="${OUT_ROOT}/tool_stdout.log"
ERR_LOG="${OUT_ROOT}/tool_stderr.log"

if [[ -n "${MODEL_DIR}" ]]; then
  TMBED_MODEL_DIR="${MODEL_DIR}"
fi

if [[ -n "${CUDA_VISIBLE_DEVICES_ARG}" ]]; then
  export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES_ARG}"
elif [[ -n "${CUDA_VISIBLE_DEVICES_DEFAULT}" ]]; then
  export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES_DEFAULT}"
fi

set +e
if [[ "${NOSEQ}" == "1" ]]; then
  python "${SCRIPT_DIR}/tmbed_predict.py" \
    --run-dir "${OUT_ROOT}" \
    --aaseq "${AASEQ}" \
    --query-id "${QUERY_ID}" \
    --model-dir "${TMBED_MODEL_DIR}" \
    --use-gpu "${USE_GPU}" \
    --threads "${THREADS}" \
    --batch-size "${BATCH_SIZE}" \
    --cuda-visible-devices "${CUDA_VISIBLE_DEVICES:-}" \
    --no-seq-in-evidence \
    > "${OUT_LOG}" 2> "${ERR_LOG}"
else
  python "${SCRIPT_DIR}/tmbed_predict.py" \
    --run-dir "${OUT_ROOT}" \
    --aaseq "${AASEQ}" \
    --query-id "${QUERY_ID}" \
    --model-dir "${TMBED_MODEL_DIR}" \
    --use-gpu "${USE_GPU}" \
    --threads "${THREADS}" \
    --batch-size "${BATCH_SIZE}" \
    --cuda-visible-devices "${CUDA_VISIBLE_DEVICES:-}" \
    > "${OUT_LOG}" 2> "${ERR_LOG}"
fi
RC=$?
set -e

if [[ $RC -ne 0 ]]; then
  cat "${ERR_LOG}" >&2
  exit $RC
fi

cat "${OUT_ROOT}/evidence.json"
