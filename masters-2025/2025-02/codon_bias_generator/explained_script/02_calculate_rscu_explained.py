#!/usr/bin/env python3
# Author: Luiza Lima Galli
# Date: 2025
# Description: Calculates Relative Synonymous Codon Usage (RSCU) and codon frequency (per 1000 codons) from a set of CDS sequences.
# Note: Some code components were generated with the help of AI tools (e.g., ChatGPT) and adapted/tested by the author.
# License: MIT

from Bio import SeqIO
from collections import Counter, defaultdict
import csv

# ---------------------------------------------------------------------------
# PARAMETERS — edit these to match your data
# ---------------------------------------------------------------------------
CDS_FILE   = "top_genes_complete.cds.fasta"  # Input: complete CDS sequences
OUTPUT_CSV = "codon_usage_rscu.csv"           # Output: RSCU table

# ---------------------------------------------------------------------------
# Standard codon table (NCBI translation table 1)
# ---------------------------------------------------------------------------
CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","AGT":"S","AGC":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"Stop","TAG":"Stop","TGA":"Stop",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGG":"W",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}

# ---------------------------------------------------------------------------
# Count codons across all CDS sequences
# ---------------------------------------------------------------------------
all_codons = []

for record in SeqIO.parse(CDS_FILE, "fasta"):
    seq = str(record.seq).replace("\n", "").replace(" ", "").upper()

    # Trim to nearest codon boundary if needed
    if len(seq) % 3 != 0:
        seq = seq[:-(len(seq) % 3)]

    all_codons.extend(seq[i:i+3] for i in range(0, len(seq), 3))

codon_counts = Counter(all_codons)
print(f"Total codons counted: {sum(codon_counts.values())}")

# ---------------------------------------------------------------------------
# Group codons by amino acid (excluding stop codons)
# ---------------------------------------------------------------------------
aa_codons = defaultdict(list)
for codon, aa in CODON_TABLE.items():
    if aa != "Stop":
        aa_codons[aa].append(codon)

# ---------------------------------------------------------------------------
# Calculate RSCU
# RSCU(codon) = observed_count / expected_count
# expected    = total_aa_count / number_of_synonymous_codons
# RSCU = 1.0  → neutral usage
# RSCU > 1.0  → preferred codon
# RSCU < 1.0  → avoided codon
# ---------------------------------------------------------------------------
rscu = {}
for aa, codons in aa_codons.items():
    total    = sum(codon_counts[c] for c in codons)
    n_codons = len(codons)
    for c in codons:
        expected = total / n_codons if total > 0 else 0
        rscu[c]  = codon_counts[c] / expected if expected > 0 else 0.0

# ---------------------------------------------------------------------------
# Write output CSV
# ---------------------------------------------------------------------------
total_codons = sum(codon_counts[c] for codons in aa_codons.values() for c in codons)

with open(OUTPUT_CSV, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["AminoAcid", "Codon", "Count", "RSCU", "Frequency_per_1000"])

    for aa in sorted(aa_codons.keys()):
        for c in sorted(aa_codons[aa]):
            freq = (codon_counts[c] / total_codons) * 1000 if total_codons > 0 else 0
            writer.writerow([aa, c, codon_counts[c], round(rscu[c], 3), round(freq, 2)])

print(f"RSCU table saved to '{OUTPUT_CSV}'")
print(f"Total coding codons analysed: {total_codons}")
