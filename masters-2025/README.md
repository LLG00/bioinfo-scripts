# codon_usage_pipeline

> A three-step pipeline to calculate Relative Synonymous Codon Usage (RSCU) from transcriptomic data.  
> Developed as part of a research project during the author's Master's degree in Bioinformatics.

---

## Description

`codon_usage_pipeline` is a Python + Bash pipeline that identifies the most highly expressed genes in a transcriptomic dataset and calculates their codon usage preferences. It computes:

- RSCU (Relative Synonymous Codon Usage) per codon
- Raw codon counts
- Codon frequency per 1000 coding codons

The pipeline produces a CSV table suitable for codon optimization of heterologous genes.

## Author

- **Luiza Lima Galli**
- Developed during the author's second semester of M.Sc. in Bioinformatics (2025)
- Contact: luiza.lima.galli@gmail.com

## Features

- Extracts top N most expressed genes from a Kallisto TPM matrix
- Handles common FASTA header format variations (versioned IDs, pipe-delimited headers)
- Filters CDS sequences to complete coding frames only
- Calculates RSCU using standard codon table (NCBI table 1)
- Outputs a clean, annotated CSV table

---

## The problem it solves

When expressing a foreign gene in a new host organism, rare codons can dramatically reduce translation efficiency. This pipeline identifies the host's most highly expressed genes (using TPM values from Kallisto), extracts their CDS sequences, and computes RSCU — a normalized metric that reveals which codons the organism prefers for each amino acid.

---

## Pipeline overview

```
Kallisto TPM matrix  +  CDS FASTA
         │
         ▼
01_extract_top_expressed.py   →  top_genes.cds.fasta
         │
         ▼
filter_complete_cds.sh        →  top_genes_complete.cds.fasta
         │
         ▼
02_calculate_rscu.py          →  codon_usage_rscu.csv
```

---

## Input files

| File | Description |
|------|-------------|
| `example_tpm.tsv` | Kallisto TPM matrix (genes × samples, tab-separated) |
| `example_cds.fasta` | CDS sequences in FASTA format |

---

## Usage

### Step 1 — Extract top expressed genes

```bash
pip install -r requirements.txt
python 01_extract_top_expressed.py
```

Edit `TOP_N` and `SAMPLE_COLUMNS` at the top of the script to match your data. Outputs `top_genes.cds.fasta`.

### Step 2 — Filter complete CDS sequences

```bash
bash filter_complete_cds.sh
```

Requires [seqkit](https://bioinf.shenwei.me/seqkit/). Outputs `top_genes_complete.cds.fasta`, keeping only sequences annotated as `complete` in their FASTA header.

### Step 3 — Calculate RSCU

```bash
python 02_calculate_rscu.py
```

Outputs `codon_usage_rscu.csv`.

---

## Output

A CSV file with one row per codon:

| Column | Description |
|--------|-------------|
| `AminoAcid` | Single-letter amino acid code |
| `Codon` | Trinucleotide codon |
| `Count` | Raw count across all input sequences |
| `RSCU` | Relative Synonymous Codon Usage (1.0 = neutral) |
| `Frequency_per_1000` | Codon frequency per 1000 coding codons |

**Interpreting RSCU:**
- `RSCU = 1.0` → neutral (codon used as expected by chance)
- `RSCU > 1.0` → preferred codon
- `RSCU < 1.0` → avoided codon

---

## ID matching strategy

The extraction script handles common FASTA header format variations automatically:

1. **Exact match** — `gene_001` matches `gene_001`
2. **First-token match** — `gene_001 description text` → tries `gene_001`
3. **Normalised match** — strips version numbers (`.1`) and pipe-delimited prefixes, lowercases

---

## Requirements

- Python 3.7+
- biopython ≥ 1.79
- pandas ≥ 1.3.0
- seqkit (for `filter_complete_cds.sh` only)

```bash
pip install -r requirements.txt
```

---

## Acknowledgments

- Some code components were developed with support from AI tools (e.g., ChatGPT) and reviewed/tested by the author.

## License

- MIT License.
- You are free to use, modify, and distribute.
- If you use this tool in your work, please cite the author (not obligatory, but appreciated).
