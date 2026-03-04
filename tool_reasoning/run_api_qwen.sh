#!/usr/bin/env bash
set -euo pipefail

source /your/path/to/conda/env

export API_KEY="sk-your-api-key"
export BASE_URL="https://dashscope.aliyuncs.com/compatible-mode/v1"
export MODEL_NAME="qwen3-max-preview"
export INCLUDE_REASONING_FIELD=1

python /your/path/to/api_qwen.py \
  --input-json /your/path/to/molinstruction_input.json \
  --out-sft-jsonl /your/path/to/infer_mol_inst_qwen_tool_reasoning.jsonl \
  --max-items 120 \
  --tool-out-root /your/path/to/runs_from_orchestrator_qwen \
  --max-rounds 12 \
  --temperature 0.0
