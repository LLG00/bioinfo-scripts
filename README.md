# hmm_domain_classifier

> A four-step pipeline for large-scale protein domain annotation and architectural classification using HMMER.  
> Developed as part of a research project during the author's Master's degree in Bioinformatics.

## Description

`hmm_domain_classifier` is a Bash + Python pipeline that annotates protein domains across large sets of genomes using Hidden Markov Models (HMM), classifies proteins by their domain architecture, and produces a per-genome summary table.

It computes:

- Domain presence per protein (from HMMER hmmsearch)
- Architectural classification (canonical vs non-canonical domain combinations)
- Per-genome protein counts per architecture type
- A filtered set of positive genomes for downstream analysis

## Author

- **Luiza Lima Galli**
- Developed during the author's M.Sc. in Bioinformatics (2025)
- Contact: luiza.lima.galli@gmail.com

## Features

- Downloads HMM profiles automatically from InterPro/Pfam and NCBI CDD
- Runs HMMER in parallel across thousands of proteomes using GNU parallel
- Classifies proteins by domain architecture using a user-defined JSON config
- Architecture definitions are fully configurable — no hardcoded domain sets
- Outputs a clean CSV table with per-genome counts per architecture type
- Copies positive genome FASTAs to a dedicated output directory

---

## The problem it solves

Simple sequence similarity searches (e.g., BLAST) often miss divergent protein variants or misclassify multidomain proteins. HMMER-based domain annotation uses profile Hidden Markov Models to detect conserved functional signatures with higher sensitivity. By combining domain co-occurrence into architectural profiles, this pipeline classifies proteins beyond simple presence/absence — distinguishing canonical multi-domain proteins from alternative single-domain or atypical variants.

---

## Pipeline overview

```
HMM profiles (Pfam/CDD)  +  Proteome FASTA files  +  architectures.json
         │
         ▼
01_download_hmm_profiles.sh   →  hmm_profiles/toxin_panel.hmm
         │
         ▼
02_run_hmmsearch.sh           →  results_hmm/*.txt + *.domtab
         │
         ▼
03_classify_architectures.py  →  positive_proteomes/ (filtered FASTAs)
         │
         ▼
04_count_per_genome.py        →  domain_counts_per_genome.csv
```

---

## Configuration

Domain architectures are defined in `architectures.json` — **edit this file to match your protein family of interest**:

```json
{
  "canonical_architectures": [
    {"name": "ArchType_A", "domains": ["PF00001", "PF00002", "PF00003"]},
    {"name": "ArchType_B", "domains": ["PF00001", "PF00002"]}
  ],
  "non_canonical_architectures": [
    {"name": "AltType_X", "domains": ["PF00005"]},
    {"name": "AltType_Y", "domains": ["PF00006"]}
  ]
}
```

Use Pfam accessions (e.g., `PF00001`) matching the profiles you downloaded in step 1.

---

## Usage

### Step 1 — Download HMM profiles

Edit the profile lists in the script to match your domain set, then run:

```bash
bash 01_download_hmm_profiles.sh
```

### Step 2 — Run HMMER in parallel

```bash
# Place your proteome .fa files in the proteomes/ directory, then:
bash 02_run_hmmsearch.sh
```

Edit `N_JOBS` in the script to match your available CPU cores.

### Step 3 — Classify by architecture

```bash
python 03_classify_architectures.py
```

### Step 4 — Count proteins per genome

```bash
python 04_count_per_genome.py
```

---

## Input files

| File/Directory | Description |
|----------------|-------------|
| `proteomes/*.fa` | Proteome FASTA files (one per genome) |
| `architectures.json` | Domain architecture definitions (user-configured) |

---

## Output files

| File/Directory | Description |
|----------------|-------------|
| `hmm_profiles/toxin_panel.hmm` | Consolidated and indexed HMM database |
| `results_hmm/*.txt` | HMMER full output per genome |
| `results_hmm/*.domtab` | HMMER domain table per genome |
| `positive_proteomes/` | FASTAs of genomes with at least one matching protein |
| `domain_counts_per_genome.csv` | Summary table: genomes × architecture counts |

---

## Requirements

- Python 3.7+
- HMMER ≥ 3.3
- GNU parallel
- wget

```bash
conda create -n hmmer_env -c conda-forge -c bioconda hmmer pandas biopython
conda activate hmmer_env
conda install -c conda-forge parallel
```

---

## Acknowledgments

- Some code components were developed with support from AI tools (e.g., ChatGPT) and reviewed/tested by the author.

## License

- MIT License.
- You are free to use, modify, and distribute.
- If you use this tool in your work, please cite the author (not obligatory, but appreciated).
