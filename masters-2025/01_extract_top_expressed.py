#!/usr/bin/env python3
# Author: Luiza Lima Galli
# Date: 2025
# Description: Extracts the top N most highly expressed CDS sequences from a
#              Kallisto TPM matrix, using multiple ID-matching strategies to
#              handle common header format variations.
# Note: Some code components were generated with the help of AI tools (e.g., ChatGPT)
#       and adapted/tested by the author.
# License: MIT

import pandas as pd
from Bio import SeqIO
import re
from collections import defaultdict

# ---------------------------------------------------------------------------
# PARAMETERS — edit these to match your data
# ---------------------------------------------------------------------------
TPM_FILE       = "example_data/example_tpm.tsv"   # Kallisto TPM matrix (tab-separated)
CDS_FASTA      = "example_data/example_cds.fasta" # CDS sequences (FASTA)
OUTPUT_FASTA   = "top_genes.cds.fasta"             # Output FASTA with top expressed genes
TOP_N          = 10                                # Number of top genes to extract
SAMPLE_COLUMNS = None                              # List of sample column names to average;
                                                   # set to None to use all columns


# ---------------------------------------------------------------------------
# ID normalisation helpers
# Handles common header variations: versioned IDs (gene.1), pipe-delimited
# headers (sp|P12345|GENE_HUMAN), and mixed-case identifiers
# ---------------------------------------------------------------------------

def strip_version(x):
    """Remove trailing version numbers (e.g. 'gene.1' → 'gene')."""
    return re.sub(r'\.\d+$', '', str(x))

def first_token(x):
    """Return first whitespace-delimited token."""
    return str(x).split()[0]

def before_pipe(x):
    """Return portion before first pipe character."""
    return str(x).split("|")[0]

def normalize_id(x):
    """Apply all normalisation steps and lowercase."""
    return strip_version(first_token(before_pipe(str(x)))).lower()


# ---------------------------------------------------------------------------
# Load TPM matrix and select top expressed genes
# ---------------------------------------------------------------------------

tpm = pd.read_csv(TPM_FILE, sep="\t", index_col=0)

if SAMPLE_COLUMNS:
    tpm = tpm[SAMPLE_COLUMNS].copy()

tpm["mean_expression"] = tpm.mean(axis=1)
top_genes_df = tpm.sort_values("mean_expression", ascending=False).head(TOP_N)
top_genes     = list(top_genes_df.index.astype(str))
top_genes_set = set(top_genes)

print(f"Selected top {len(top_genes)} genes by mean TPM.")

# Build normalised lookup table for flexible matching
norm_to_original = defaultdict(list)
for g in top_genes:
    norm_to_original[normalize_id(g)].append(g)


# ---------------------------------------------------------------------------
# Match TPM gene IDs to CDS FASTA sequences
# Strategy: exact → first-token → normalised
# ---------------------------------------------------------------------------

found_records = {}

with open(CDS_FASTA) as fasta_in:
    for record in SeqIO.parse(fasta_in, "fasta"):
        rec_id      = record.id
        rec_id_norm = normalize_id(rec_id)
        matched     = None

        if rec_id in top_genes_set:
            matched = rec_id                                      # exact match

        if matched is None:
            tok = first_token(rec_id)
            if tok in top_genes_set:
                matched = tok                                     # first-token match

        if matched is None and rec_id_norm in norm_to_original:
            matched = norm_to_original[rec_id_norm][0]           # normalised match

        if matched and matched not in found_records:
            found_records[matched] = record


# ---------------------------------------------------------------------------
# Write output FASTA preserving original TPM rank order
# ---------------------------------------------------------------------------

with open(OUTPUT_FASTA, "w") as out_f:
    written = 0
    for g in top_genes:
        if g in found_records:
            SeqIO.write(found_records[g], out_f, "fasta")
            written += 1

print(f"Extracted {written}/{len(top_genes)} sequences → {OUTPUT_FASTA}")
if written < len(top_genes):
    missing = [g for g in top_genes if g not in found_records]
    print(f"  Could not match: {missing}")
