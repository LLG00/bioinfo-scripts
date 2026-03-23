# blast_coverage_parser

> A BLAST XML parser that calculates real query and subject coverage by merging HSP intervals.  
> Developed as part of a research project during the author's Master's degree in Bioinformatics.

## Description

`blast_coverage_parser` is a Python script designed to process BLAST XML output and calculate accurate coverage metrics for query and subject sequences. It merges overlapping HSPs before computing coverage, and produces a tab-separated summary with:

- Query coverage (%)
- Subject coverage (%)
- Identity
- Bit score
- E-value

## Author

- **Luiza Lima Galli**
- Developed during the author's first semester of M.Sc. in Bioinformatics (2025)
- Contact: luiza.lima.galli@gmail.com

## Features

- Parses BLAST XML output using Biopython
- Merges overlapping or adjacent HSP intervals to calculate true, non-redundant coverage
- Handles multiple hits and alignments per query
- Handles reverse-strand hits by normalising coordinate order
- Handles pipe-delimited GenBank-style hit headers automatically
- Outputs a clean, tab-delimited report

---

## The problem it solves

When a query aligns to a subject through multiple High-Scoring Pairs (HSPs), naive coverage calculation sums all HSP lengths independently, leading to **overestimated coverage** in regions of overlap. This script merges overlapping or adjacent intervals before computing coverage, producing accurate per-alignment coverage values.

---

## Input files

| File | Description |
|------|-------------|
| `query.fasta` | FASTA file with query sequences |
| `subject.fasta` | FASTA file with subject (database) sequences |
| `blast_results.xml` | BLAST output in XML format (`-outfmt 5`) |

### Generating BLAST XML output

```bash
blastp \
  -query query.fasta \
  -subject subject.fasta \
  -outfmt 5 \
  -out blast_results.xml
```

---

## Usage

```bash
# Install dependencies
pip install -r requirements.txt

# Run with example data
python blast_coverage_parser.py
```

By default, the script reads from `example_data/` and writes to `output_coverage.tsv`. Edit the file paths at the top of the script to point to your own data.

---

## Output

A tab-separated file (`output_coverage.tsv`) with one row per alignment:

| Column | Description |
|--------|-------------|
| `query_id` | Query sequence ID |
| `query_length` | Full length of query sequence (aa or nt) |
| `hit_id` | Subject/hit sequence ID |
| `hit_length` | Full length of subject sequence |
| `identity` | Fraction of identical residues across HSPs |
| `bit_score` | Bit score of the alignment |
| `e_value` | E-value of the alignment |
| `query_start` / `query_end` | Alignment coordinates on query |
| `hit_start` / `hit_end` | Alignment coordinates on subject |
| `query_coverage(%)` | Non-redundant query coverage (%) |
| `subject_coverage(%)` | Non-redundant subject coverage (%) |

### Example output

```
query_id  query_length  hit_id  hit_length  identity   bit_score  e_value   query_start  query_end  hit_start  hit_end  query_coverage(%)  subject_coverage(%)
seq1      138           hit1    138         0.978261   280.00     1.00e-98  1            138        1          138      100.00             100.00
seq1      138           hit2    138         0.833333   150.00     3.00e-45  5            100        10         105      69.57              69.57
```

---

## Requirements

- Python 3.7+
- biopython ≥ 1.79

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
