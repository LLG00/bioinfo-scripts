# blast_coverage_parser

> A BLAST XML parser that calculates real query and subject coverage by merging HSP intervals.  
> Developed as part of a research project during the author's Master's degree in Bioinformatics.

## Description

`blast_coverage_parser` is a Python script designed to process BLAST XML output and calculate more accurate coverage metrics for query and subject sequences. It merges overlapping HSPs and computes:

- Query coverage (%)
- Subject coverage (%)
- Identity
- Bit score
- E-value

The script produces a tab-separated summary of all relevant alignment data.

## Author

- **Luiza Lima Galli**
- Developed during the author's first year's M.Sc. in Bioinformatics
- Contact: luiza.lima.galli@gmail.com

## Features

- Parses BLAST XML output using Biopython
- Merges overlapping HSPs to calculate true coverage
- Handles multiple hits and alignments per query
- Outputs a clean and customizable tabular report

## Input files

- `example_BLAST.xml` — BLAST output file in XML format
- `query_FASTA_example.fasta` — Original query sequences (FASTA format)
- `subject_FASTA_example.fasta` — Subject (hit) sequences (FASTA format)

## Output

- `example_OUTPUT.txt` — Tab-delimited table with all relevant alignment statistics

## Acknowledgments
- Some code components were developed with support from AI tools (e.g., ChatGPT) and reviewed/tested by the author.

## License
- MIT License.
- You are free to use, modify, and distribute.
