#!/usr/bin/env bash
# Author: Luiza Lima Galli
# Date: 2025
# Description: Filters a CDS FASTA file to retain only complete coding sequences (containing "complete" and not "incomplete" in the header), then selects the top N sequences for codon usage analysis.
# Usage: bash filter_complete_cds.sh
# Dependencies: seqkit (https://bioinf.shenwei.me/seqkit/)
# License: MIT

set -euo pipefail

INPUT_FASTA="top_genes.cds.fasta" 
FILTERED_FASTA="top_genes_complete.cds.fasta"
TOP_N=200 

echo "Filtering complete CDS sequences..."
awk '/^>/{p=($0 ~ /complete/ && $0 !~ /incomplete/)} p' \
    "$INPUT_FASTA" > "$FILTERED_FASTA"

N_COMPLETE=$(grep -c "^>" "$FILTERED_FASTA" || true)
echo "  Complete sequences found: $N_COMPLETE"

echo "Selecting top $TOP_N sequences..."
seqkit head -n "$TOP_N" "$FILTERED_FASTA" > tmp_head.fasta
mv tmp_head.fasta "$FILTERED_FASTA"

echo "Done → $FILTERED_FASTA"
