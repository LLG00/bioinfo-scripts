# pep2cds

> A codon-based CDS alignment tool guided by PEP protein alignments.  
> Developed during an academic exchange at Technische Universität Braunschweig (2024),  
> as part of the Final Project at the [Python for Life Scientists Course (BB34)](https://www.tu-braunschweig.de/en/ifp/pbb/teaching/pythoncourseprojects).

## Description

`pep2cds` is a Python script that reconstructs codon-preserving CDS alignments from protein (PEP) alignments. It maps codons from a nucleotide FASTA file using the corresponding translated protein sequences, and aligns them according to the input protein alignment.

This tool is useful when codon-aware alignments are required for downstream analyses such as:

- Selection detection (dN/dS)
- Evolutionary and codon substitution modeling
- Nucleotide-based phylogenetics (e.g., IQ-TREE, HyPhy)

## Author

- **Luiza Lima Galli**
- Developed during an academic exchange at Technische Universität Braunschweig (2024)
- Contact: luiza.lima.galli@gmail.com

## Features

- Back-translates a PEP alignment into a codon-based CDS alignment
- Maps each amino acid to its corresponding codon(s) from the original CDS sequences
- Expands gap characters (`-`) in the PEP alignment to codon gaps (`---`)
- Handles ambiguous or overlapping codon assignments via synchronisation
- Performs built-in sanity checks on input files
- No external libraries required — pure Python standard library

---

## The problem it solves

Tools like MUSCLE or MAFFT align sequences at the amino acid level, which is more accurate than aligning raw nucleotides. However, downstream analyses (e.g., dN/dS estimation, codon-model phylogenetics with IQ-TREE or HyPhy) require nucleotide alignments. This script back-translates a PEP alignment into codons, keeping gaps as `---` triplets, so the codon alignment is fully consistent with the protein alignment.

---

## Input files

| File | Description |
|------|-------------|
| `example_CDS.fasta` | Unaligned nucleotide CDS sequences (must start with ATG) |
| `example_PEP.fasta` | Unaligned amino acid sequences corresponding to the CDS |
| `example_PEP_align.fasta` | Aligned amino acid sequences (output from MUSCLE, MAFFT, etc.) |

> **Important:** sequence names must match exactly across all three files.

### Generating the PEP alignment (example with MUSCLE)

```bash
muscle -align example_PEP.fasta -output example_PEP_align.fasta
```

---

## Usage

```bash
# No installation needed — pure Python standard library
python pep2cds.py
```

By default, the script reads from `example_data/` and writes to `output_cds_align.fasta`. Edit the file paths at the bottom of the script to point to your own data.

---

## Output

A FASTA file (`output_cds_align.fasta`) where each sequence is the codon-based alignment corresponding to the input PEP alignment:

```
>seqA
ATGAAAGCGTTCGCC---TTTGTTGGCGGT...
>seqB
ATGAAAGCGTTTGCC---TTTGTTGGCGGT...
>seqC
ATGAAAGCGTTCGCGTTCGTTGGCGGCACG...
```

Gap columns (`-`) in the PEP alignment are expanded to codon gaps (`---`).

---

## Sanity checks

The script automatically warns if:

- A sequence present in the PEP file has no corresponding CDS entry
- A CDS sequence length is not exactly 3× its PEP length
- A CDS sequence does not start with `ATG`
- An amino acid in the alignment has no mapped codon

---

## Requirements

- Python 3.7+
- No external libraries required

---

## Acknowledgments

- Some parts of this script were developed with support from AI tools (e.g., ChatGPT) and reviewed/tested by the author.

## License

- MIT License.
- You are free to use, modify, and distribute.
- If you use this tool in your work, please cite the author (not obligatory, but appreciated).
