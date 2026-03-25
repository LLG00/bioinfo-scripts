#!/usr/bin/env python3
# Author: Luiza Lima Galli
# Date: 2026
# Description: Classifies proteins by domain architecture using HMMER results, identifies positive genomes, and copies their FASTA files to an output directory. Architectures are loaded from an external JSON config file.
# Note: Some code components were generated with the help of AI tools (e.g., ChatGPT) and adapted/tested by the author.
# License: MIT

import os
import json
import shutil
from collections import Counter

# ---------------------------------------------------------------------------
# PARAMETERS — edit these to match your environment
# ---------------------------------------------------------------------------
RESULTS_DIR    = "results_hmm"          # Directory with HMMER .txt output files
PROTEOMES_DIR  = "proteomes"            # Directory with original proteome .fa files
OUTPUT_DIR     = "positive_proteomes"   # Output directory for positive genomes
ARCH_CONFIG    = "architectures.json"   # Domain architecture config file

# ---------------------------------------------------------------------------
# Load architecture definitions from config file
# ---------------------------------------------------------------------------
with open(ARCH_CONFIG) as f:
    config = json.load(f)

canonical_architectures     = [set(a["domains"]) for a in config["canonical_architectures"]]
non_canonical_architectures = {a["name"]: set(a["domains"]) for a in config["non_canonical_architectures"]}

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Parse HMMER results and classify proteins
# ---------------------------------------------------------------------------
validated_proteomes_count = 0
detailed_stats = Counter()

for filename in os.listdir(RESULTS_DIR):
    if not filename.endswith(".txt"):
        continue

    proteins = {}
    current_pfam_id = None

    # Parse HMMER output: map domain accessions to each protein
    with open(os.path.join(RESULTS_DIR, filename)) as f:
        for line in f:
            if line.startswith("Accession:"):
                current_pfam_id = line.split()[1].split('.')[0]  # Strip version suffix

            if line.startswith(">>"):
                prot_id = line.split()[1]
                if prot_id not in proteins:
                    proteins[prot_id] = set()
                if current_pfam_id:
                    proteins[prot_id].add(current_pfam_id)

    # Classify each protein by its domain set
    is_positive_genome = False

    for domain_set in proteins.values():
        if domain_set in canonical_architectures:
            detailed_stats["TOTAL_CANONICAL"] += 1
            is_positive_genome = True
        else:
            for name, fset in non_canonical_architectures.items():
                if domain_set == fset:
                    detailed_stats["TOTAL_NON_CANONICAL"] += 1
                    is_positive_genome = True
                    break

    # Copy positive genome FASTA to output directory
    if is_positive_genome:
        target_fasta = filename.replace(".pep.all.txt", ".fa")
        src  = os.path.join(PROTEOMES_DIR, target_fasta)
        dest = os.path.join(OUTPUT_DIR, target_fasta)

        if os.path.exists(src):
            shutil.copy2(src, dest)
            validated_proteomes_count += 1

print(f"Positive genomes identified: {validated_proteomes_count}")
print(f"Canonical proteins:          {detailed_stats['TOTAL_CANONICAL']}")
print(f"Non-canonical proteins:      {detailed_stats['TOTAL_NON_CANONICAL']}")
