# PFUA
Implementation and tooling setup notes for the paper: **Interleaved Tool-Call Reasoning for Protein Function Understanding** (https://arxiv.org/abs/2601.03604)


---

## Environment

### 1) Basic environment

This project follows the **verl** base environment. (https://github.com/volcengine/verl)

```bash
git clone https://github.com/volcengine/verl.git
cd verl

conda create -n verl python==3.12
conda activate verl

USE_MEGATRON=0 bash scripts/install_vllm_sglang_mcore.sh
pip install --no-deps -e .
```

### 2) Extra dependencies

```bash
# MMseqs2
conda install -c bioconda -c conda-forge mmseqs2

# TMbed
python -m pip install git+https://github.com/BernhoferM/TMbed.git

# HMMER (for Pfam)
conda install -c bioconda hmmer -y
```

### 3) Optional reference

We provide an `environment.yml` in this repo.

---

## Tool Setup
### 1) MMseqs2: build `targetDB` from Swiss-Prot FASTA

Download Swiss-Prot FASTA:

```bash
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
```

Build the MMseqs2 database:

```bash
cd /your/path/to/target_db
gunzip uniprot_sprot.fasta.gz
mmseqs createdb uniprot_sprot.fasta targetDB
```
### 2) Pfam: download `Pfam-A.hmm` and run `hmmpress`

Download and unpack:

```bash
wget -c https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
ls -lh Pfam-A.hmm.gz
gunzip -k Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```


### 3) TMbed

Install:

```bash
python -m pip install git+https://github.com/BernhoferM/TMbed.git
```

**First run note:** TMbed will download **ProtT5 (~2.25GB)** automatically on first execution.

---

## Protein Database Note

The protein database used for alignment/search is based on **UniProtKB Swiss-Prot**: (https://www.uniprot.org/uniprotkb?query=reviewed:true)
- Download cutoff date: **year 2025**
- Total annotations: **573,661** entries

---

## Usage
```bash
bash /your/path/to/tool_reasoning/run_api_qwen.sh
bash /your/path/to/tool_reasoning/run_api_kimi.sh
bash /your/path/to/tool_reasoning/run_api_deepseek.sh
```
---

## Citation
```bibtex
@misc{fan2026interleavedtoolcallreasoningprotein,
      title={Interleaved Tool-Call Reasoning for Protein Function Understanding},
      author={Chuanliu Fan and Zicheng Ma and Huanran Meng and Aijia Zhang and Wenjie Du and Jun Zhang and Yi Qin Gao and Ziqiang Cao and Guohong Fu},
      year={2026},
      eprint={2601.03604},
      archivePrefix={arXiv},
      primaryClass={cs.AI},
      url={https://arxiv.org/abs/2601.03604},
}
```

---

## Acknowledgements
- verl, MMseqs2, HMMER, Pfam, TMbed, Mol-Instructions, ProtTeX, and the UniProt Consortium for tooling and datasets.
