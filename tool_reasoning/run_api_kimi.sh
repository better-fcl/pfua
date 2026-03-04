#!/usr/bin/env bash
set -euo pipefail

source /your/path/to/conda/env

# Kimi API
export API_KEY="sk-your-api-key"
export BASE_URL="https://api.moonshot.cn/v1"
export MODEL_NAME="kimi-k2-thinking"

python /your/path/to/api_kimi.py \
  --input-json /your/path/to/molinstruction_input.json \
  --out-sft-jsonl /your/path/to/infer_mol_inst_kimi_tool_reasoning.jsonl \
  --max-items 120 \
  --tool-out-root /your/path/to/runs_from_orchestrator \
  --max-rounds 12 \
  --temperature 0.0
