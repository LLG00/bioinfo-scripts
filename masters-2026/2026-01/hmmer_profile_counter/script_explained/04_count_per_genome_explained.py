#!/usr/bin/env python3
# Author: Luiza Lima Galli
# Date: 2025
# Description: Counts proteins per domain architecture per genome from HMMER results,
#              and exports a summary CSV table. Architectures are loaded from an
#              external JSON config file.
# Note: Some code components were generated with the help of AI tools (e.g., ChatGPT)
#       and adapted/tested by the author.
# License: MIT

import os
import csv
import json
from collections import Counter


# ---------------------------------------------------------------------------
# PARAMETERS — edit these to match your environment
# ---------------------------------------------------------------------------
RESULTS_DIR = "results_hmm"             # Directory with HMMER .txt output files
ARCH_CONFIG = "architectures.json"      # Domain architecture config file
OUTPUT_CSV  = "domain_counts_per_genome.csv"  # Output summary table


# ---------------------------------------------------------------------------
# Load architecture definitions from config file
# ---------------------------------------------------------------------------
with open(ARCH_CONFIG) as f:
    config = json.load(f)

canonical_architectures = {
    a["name"]: set(a["domains"]) for a in config["canonical_architectures"]
}
non_canonical_architectures = {
    a["name"]: set(a["domains"]) for a in config["non_canonical_architectures"]
}

all_arch_names = list(canonical_architectures.keys()) + list(non_canonical_architectures.keys())
all_columns    = ["Genome"] + all_arch_names
results_data   = []


# ---------------------------------------------------------------------------
# Parse HMMER results and count proteins per architecture per genome
# ---------------------------------------------------------------------------
for filename in os.listdir(RESULTS_DIR):
    if not filename.endswith(".txt"):
        continue

    genome_name     = filename.replace(".pep.all.txt", "")
    proteins        = {}
    current_pfam_id = None

    # Map domain accessions to each protein
    with open(os.path.join(RESULTS_DIR, filename)) as f:
        for line in f:
            if line.startswith("Accession:"):
                current_pfam_id = line.split()[1].split('.')[0]
            if line.startswith(">>"):
                prot_id = line.split()[1]
                if prot_id not in proteins:
                    proteins[prot_id] = set()
                if current_pfam_id:
                    proteins[prot_id].add(current_pfam_id)

    # Count matches per architecture
    genome_counts = Counter()

    for domain_set in proteins.values():
        match_found = False

        for name, arch_set in canonical_architectures.items():
            if domain_set == arch_set:
                genome_counts[name] += 1
                match_found = True
                break

        if not match_found:
            for name, arch_set in non_canonical_architectures.items():
                if domain_set == arch_set:
                    genome_counts[name] += 1
                    break

    # Record only genomes with at least one match
    if sum(genome_counts.values()) > 0:
        row = {"Genome": genome_name}
        for col in all_arch_names:
            row[col] = genome_counts.get(col, 0)
        results_data.append(row)


# ---------------------------------------------------------------------------
# Export results to CSV
# ---------------------------------------------------------------------------
with open(OUTPUT_CSV, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=all_columns)
    writer.writeheader()
    for data in sorted(results_data, key=lambda x: x['Genome']):
        writer.writerow(data)

print(f"Output saved to: {OUTPUT_CSV}")
print(f"Genomes with matches: {len(results_data)}")
