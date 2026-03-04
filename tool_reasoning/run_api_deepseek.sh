#!/usr/bin/env bash
set -euo pipefail

source /your/path/to/conda/env

# DeepSeek API
export DEEPSEEK_API_KEY="sk-your-api-key"
export DEEPSEEK_BASE_URL="https://api.deepseek.com"
export DEEPSEEK_MODEL_NAME="deepseek-reasoner"

python /your/path/to/api_deepseek.py \
  --input-json /your/path/to/molinstruction_input.json \
  --out-sft-jsonl /your/path/to/infer_mol_inst_deepseek_tool_reasoning.jsonl \
  --max-items 120 \
  --tool-out-root /your/path/to/runs_from_orchestrator \
  --max-rounds 12 \
  --temperature 0.0
