#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import os
import re
import subprocess
import time
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from openai import OpenAI


# ==============
# Tool entrypoints
# ==============
TOOL_BASH: Dict[str, str] = {
    "seq_basic_props": "/your/path/to/seq_basic_props/run_seq_basic_props.sh",
    "mmseqs2_besthit_uniprot": "/your/path/to/mmseqs2/run_mmseqs2_besthit_uniprot.sh",
    "pfam_hmmscan": "/your/path/to/pfam_hmmscan/run_pfam_hmmscan.sh",
    "tmbed_predict": "/your/path/to/tmbed/run_tmbed_predict.sh",
}

# =========================
# Tool schemas
# =========================
TOOLS_SPEC: List[Dict[str, Any]] = [
    {
        "type": "function",
        "function": {
            "name": "seq_basic_props",
            "description": "Compute basic sequence properties such as length, composition, charge, low-complexity and hydrophobicity.",
            "parameters": {
                "type": "object",
                "properties": {"sequence_ref": {"type": "string"}},
                "additionalProperties": False,
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "mmseqs2_besthit_uniprot",
            "description": "Search UniProt to identify closest homolog and functional annotation.",
            "parameters": {
                "type": "object",
                "properties": {
                    "sequence_ref": {"type": "string"},
                    "min_seq_id": {"type": "number", "default": 0.3},
                },
                "additionalProperties": False,
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "pfam_hmmscan",
            "description": "Detect conserved domains using Pfam HMMs and summarize predicted domain architecture.",
            "parameters": {
                "type": "object",
                "properties": {"sequence_ref": {"type": "string"}},
                "additionalProperties": False,
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "tmbed_predict",
            "description": "Predict transmembrane regions and topology.",
            "parameters": {
                "type": "object",
                "properties": {"sequence_ref": {"type": "string"}},
                "additionalProperties": False,
            },
        },
    },
]

# =========================
# Prompt blocks
# =========================
COLLECT_HEADER_WITH_ROLE_SEPS = """[ROLE] You are an expert protein analysis agent.

[TASK]
Your goal is to analyze the protein sequence and produce a biologically meaningful interpretation.
You should reason step-by-step, form hypotheses, and use tools only when they help reduce uncertainty.

[REASONING REQUIREMENTS]
Before calling tools, you MUST:
- propose hypotheses about the protein
- explain which uncertainties still remain

For EVERY tool call, you MUST:
- explicitly explain WHY this tool is needed
- describe WHAT evidence you expect it to provide

After each tool result, you MUST:
- summarize what new evidence was obtained
- update or revise your hypothesis
- decide whether additional tools are needed

[TOOLS]
You may call the following tools through function calling:
- seq_basic_props: basic physicochemical properties
- pfam_hmmscan: domain and family inference
- mmseqs2_besthit_uniprot: homolog search and functional annotation
- tmbed_predict: transmembrane and topology prediction

Prefer {"sequence_ref": "query"} instead of pasting long sequences.

[OUTPUT]
When finished, wrap the final answer in <answer>...</answer>.
"""

SFT_HEADER_PLAIN = COLLECT_HEADER_WITH_ROLE_SEPS
AA_RE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]+", re.I)


def extract_aaseq(input_field: str) -> str:
    m = re.search(r"```(.*?)```", input_field, flags=re.S)
    body = m.group(1) if m else input_field
    parts = AA_RE.findall(body)
    if not parts:
        raise ValueError("Failed to extract amino-acid sequence.")
    return "".join(parts).upper()


def is_final_answer(text: str) -> bool:
    t = (text or "").lower()
    return "<answer>" in t and "</answer>" in t


def build_user_collection(instruction: str, aaseq: str) -> str:
    return (
        COLLECT_HEADER_WITH_ROLE_SEPS
        + "\n<|im_start|>user\n"
        + instruction.strip()
        + "\nProtein sequence (id=query):\n```"
        + aaseq
        + "```\n<|im_end|>\n"
    )


def build_user_sft(instruction: str, aaseq: str) -> str:
    return (
        SFT_HEADER_PLAIN
        + "\n\n"
        + instruction.strip()
        + "\nProtein sequence (id=query):\n```"
        + aaseq
        + "```"
    )


def normalize_tool_calls(tool_calls_obj: Any) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for tc in tool_calls_obj:
        arg_str = tc.function.arguments or "{}"
        try:
            args = json.loads(arg_str)
        except Exception:
            args = {"_raw_arguments": arg_str}
        out.append({"function": {"arguments": args, "name": tc.function.name}, "type": tc.type})
    return out


def run_tool_via_bash(
    tool_name: str,
    tool_args: Dict[str, Any],
    seq_map: Dict[str, str],
    out_root: Path,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    if tool_name not in TOOL_BASH:
        raise KeyError(f"Unknown tool: {tool_name}")
    bash_path = Path(TOOL_BASH[tool_name])
    if not bash_path.exists():
        raise FileNotFoundError(f"Tool bash not found: {bash_path}")

    aaseq = (tool_args.get("aaseq") or "").strip()
    if not aaseq:
        ref = (tool_args.get("sequence_ref") or "query").strip()
        if ref not in seq_map:
            raise KeyError(f"sequence_ref '{ref}' not found")
        aaseq = seq_map[ref]

    query_id = (tool_args.get("sequence_ref") or "query") if tool_args.get("sequence_ref") else "query"

    call_id = f"{tool_name}_{time.strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    run_dir = out_root / call_id
    run_dir.mkdir(parents=True, exist_ok=True)

    cmd = ["bash", str(bash_path), "--aaseq", aaseq, "--query-id", str(query_id), "--out-root", str(run_dir)]

    def add_flag(name: str, flag: str, is_bool: bool = False):
        if name in tool_args and tool_args[name] is not None:
            v = tool_args[name]
            if is_bool:
                if bool(v):
                    cmd.append(flag)
            else:
                cmd.extend([flag, str(v)])

    if tool_name == "seq_basic_props":
        add_flag("membrane_run_thr", "--membrane-run-thr")
        add_flag("low_complexity_thr", "--low-complexity-thr")
        add_flag("include_cleaning_notes", "--include-cleaning-notes", is_bool=True)

    elif tool_name == "mmseqs2_besthit_uniprot":
        add_flag("threads", "--threads")
        add_flag("evalue", "--evalue")
        add_flag("min_seq_id", "--min-seq-id")
        add_flag("coverage", "--coverage")
        add_flag("target_db", "--target-db")
        add_flag("uniprot_dat", "--uniprot-dat")
        add_flag("no_seq_in_evidence", "--no-seq-in-evidence", is_bool=True)

    elif tool_name == "pfam_hmmscan":
        add_flag("pfam_hmm", "--pfam-hmm")
        add_flag("cpu", "--cpu")
        add_flag("incE", "--incE")
        add_flag("domE", "--domE")
        add_flag("select_evalue", "--select-evalue")
        add_flag("select_min_cov", "--select-min-cov")
        add_flag("topk", "--topk")
        add_flag("no_seq_in_evidence", "--no-seq-in-evidence", is_bool=True)

    elif tool_name == "tmbed_predict":
        add_flag("model_dir", "--model-dir")
        add_flag("use_gpu", "--use-gpu")
        add_flag("threads", "--threads")
        add_flag("batch_size", "--batch-size")
        add_flag("cuda_visible_devices", "--cuda-visible-devices")
        add_flag("no_seq_in_evidence", "--no-seq-in-evidence", is_bool=True)

    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    stdout = p.stdout.decode("utf-8", errors="replace").strip()
    stderr = p.stderr.decode("utf-8", errors="replace").strip()

    meta = {
        "cmd": cmd,
        "returncode": p.returncode,
        "run_dir": str(run_dir),
        "stderr": stderr[:4000],
    }

    if p.returncode != 0:
        raise RuntimeError(
            f"Tool failed ({p.returncode}): {bash_path}\nCMD: {' '.join(cmd)}\nSTDOUT:\n{stdout}\nSTDERR:\n{stderr}"
        )

    try:
        obj = json.loads(stdout)
    except Exception as e:
        raise ValueError(f"Cannot parse tool stdout JSON for {tool_name}: {e}\n---STDOUT---\n{stdout[:4000]}")
    return obj, meta


def get_model_texts(msg: Any) -> Tuple[str, str, str]:
    """
    Returns:
      - raw_content: msg.content (may be empty)
      - raw_reasoning: msg.reasoning_content (may be empty)
      - merged_for_sft: reasoning + content (preferred for training to ensure non-empty)
    """
    raw_content = getattr(msg, "content", None)
    raw_reasoning = getattr(msg, "reasoning_content", None)

    c = raw_content if isinstance(raw_content, str) else ""
    r = raw_reasoning if isinstance(raw_reasoning, str) else ""

    if r and c:
        merged = r + "\n\n" + c
    elif r:
        merged = r
    else:
        merged = c

    return c, r, merged


def llm_call(client: OpenAI, model: str, messages: List[Dict[str, Any]], temperature: float):
    return client.chat.completions.create(
        model=model,
        messages=messages,
        tools=TOOLS_SPEC,
        tool_choice="auto",
        temperature=temperature,
        max_tokens=2048,
        stream=False,
    )


def collect_one_sample(
    client: OpenAI,
    model: str,
    instruction: str,
    aaseq: str,
    sample_tool_out_root: Path,
    max_rounds: int,
    temperature: float,
    gold_output: Optional[str] = None,
    gold_metadata: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    seq_map = {"query": aaseq}

    user_collect = build_user_collection(instruction, aaseq)
    user_sft = build_user_sft(instruction, aaseq)

    api_messages: List[Dict[str, Any]] = [
        {"role": "system", "content": "You are a helpful assistant."},
        {"role": "user", "content": user_collect},
    ]

    sft_messages: List[Dict[str, Any]] = [{"role": "user", "content": user_sft, "tool_calls": None}]

    for _ in range(max_rounds):
        resp = llm_call(client, model, api_messages, temperature=temperature)
        msg = resp.choices[0].message

        raw_content, raw_reasoning, merged_for_sft = get_model_texts(msg)
        tool_calls = getattr(msg, "tool_calls", None)

        if tool_calls:
            sft_messages.append(
                {
                    "role": "assistant",
                    "content": merged_for_sft.strip(),
                    "tool_calls": normalize_tool_calls(tool_calls),
                    "reasoning_content": raw_reasoning if raw_reasoning else None,
                }
            )

            api_messages.append(
                {
                    "role": "assistant",
                    "content": raw_content,
                    "reasoning_content": raw_reasoning,
                    "tool_calls": [
                        {
                            "id": tc.id,
                            "type": tc.type,
                            "function": {"name": tc.function.name, "arguments": tc.function.arguments},
                        }
                        for tc in tool_calls
                    ],
                }
            )

            for tc in tool_calls:
                try:
                    args = json.loads(tc.function.arguments or "{}")
                except Exception:
                    args = {"_raw_arguments": tc.function.arguments or ""}

                tool_json, tool_meta = run_tool_via_bash(
                    tool_name=tc.function.name,
                    tool_args=args,
                    seq_map=seq_map,
                    out_root=sample_tool_out_root,
                )
                tool_resp_text = json.dumps(tool_json, ensure_ascii=False)

                sft_messages.append({"role": "tool", "content": tool_resp_text, "tool_calls": None})

                api_messages.append({"role": "tool", "tool_call_id": tc.id, "content": tool_resp_text})

            continue

        sft_messages.append(
            {
                "role": "assistant",
                "content": merged_for_sft.strip(),
                "tool_calls": None,
                "reasoning_content": raw_reasoning if raw_reasoning else None,
            }
        )

        api_messages.append(
            {
                "role": "assistant",
                "content": raw_content,
                "reasoning_content": raw_reasoning,
            }
        )

        if is_final_answer(merged_for_sft):
            break

    sft_record = {
        "messages": sft_messages,
        "tools": TOOLS_SPEC,
        "gold_output": gold_output,
        "gold_metadata": gold_metadata or {},
    }
    return sft_record


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-json", required=True)
    ap.add_argument("--out-sft-jsonl", required=True)
    ap.add_argument("--max-items", type=int, default=0)
    ap.add_argument("--tool-out-root", required=True)
    ap.add_argument("--max-rounds", type=int, default=12)
    ap.add_argument("--temperature", type=float, default=0.0)
    args = ap.parse_args()

    api_key = os.environ.get("DEEPSEEK_API_KEY", "").strip()
    base_url = os.environ.get("DEEPSEEK_BASE_URL", "https://api.deepseek.com").strip()
    model = os.environ.get("DEEPSEEK_MODEL_NAME", "deepseek-reasoner").strip()
    if not api_key:
        raise ValueError("Please export DEEPSEEK_API_KEY.")

    client = OpenAI(api_key=api_key, base_url=base_url)

    data = json.loads(Path(args.input_json).read_text(encoding="utf-8"))
    if args.max_items and args.max_items > 0:
        data = data[: args.max_items]

    out_sft = Path(args.out_sft_jsonl)
    out_sft.parent.mkdir(parents=True, exist_ok=True)
    out_sft.touch(exist_ok=True)

    tool_out_root = Path(args.tool_out_root)
    tool_out_root.mkdir(parents=True, exist_ok=True)

    for i, item in enumerate(data):
        try:
            instruction = (item.get("instruction") or "").strip()
            inp = (item.get("input") or "").strip()
            if not instruction or not inp:
                print(f"[WARN] skip item {i}: missing instruction/input")
                continue

            aaseq = extract_aaseq(inp)
            sample_dir = tool_out_root / f"sample_{i:06d}"
            sample_dir.mkdir(parents=True, exist_ok=True)

            sft_record = collect_one_sample(
                client=client,
                model=model,
                instruction=instruction,
                aaseq=aaseq,
                sample_tool_out_root=sample_dir,
                max_rounds=args.max_rounds,
                temperature=args.temperature,
                gold_output=item.get("output"),
                gold_metadata=item.get("metadata", {}),
            )

            with out_sft.open("a", encoding="utf-8") as f:
                f.write(json.dumps(sft_record, ensure_ascii=False) + "\n")

            print(f"✅ done {i+1}/{len(data)}")

        except Exception as e:
            print(f"❌ item {i} failed: {e}")
            continue


if __name__ == "__main__":
    main()
