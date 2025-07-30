# pep2cds

> A codon-based CDS alignment tool based on PEP protein alignments.  
> Developed during an academic exchange at Technische Universität Braunschweig (2024),  
> as part of the Python Course Project at the [Institute of Biochemical Engineering (IFP)](https://www.tu-braunschweig.de/en/ifp/pbb/teaching/pythoncourseprojects).

## Description

`pep2cds` is a Python script that reconstructs codon-preserving CDS alignments from protein (PEP) alignments. It maps codons from a nucleotide FASTA file using the corresponding translated protein sequences, and aligns them according to the input protein alignment.

This tool is helpful when codon-aware alignments are required for downstream analyses such as selection detection, evolutionary modeling, or nucleotide-based phylogenetics.

## Author

- **Luiza Lima Galli**
- Developed during an academic internship in 2024 at TU Braunschweig, Germany
- Contact: luiza.lima.galli@gmail.com

## Features

- Validates CDS and PEP sequence lengths
- Maps each amino acid to its codon(s)
- Reconstructs a CDS alignment from a protein alignment, preserving gaps
- Handles inconsistent or incomplete input with clear warnings

## Input files

- `example_CDS.fasta` — CDS sequences (nucleotide)
- `example_PEP.fasta` — Protein sequences (translated from CDS)
- `example_PEP_align.tx

## Output file

- `example_CDS_align.txt` — Reconstructed codon-aware alignment in nucleotide format

## How to run

```bash
python pep2cds.py
