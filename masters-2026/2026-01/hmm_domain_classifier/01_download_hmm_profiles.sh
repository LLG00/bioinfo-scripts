#!/usr/bin/env bash
# Author: Luiza Lima Galli
# Date: 2025
# Description: Downloads HMM profiles from InterPro/Pfam and NCBI CDD, and consolidates them into a single indexed database for HMMER search.
# Usage: bash 01_download_hmm_profiles.sh
# Dependencies: wget, hmmer, gunzip
# License: MIT

set -euo pipefail

OUTPUT_DB="hmm_profiles/toxin_panel.hmm"
PROFILES_DIR="hmm_profiles"

mkdir -p "$PROFILES_DIR"

echo "Downloading Pfam profiles from InterPro..."

declare -A PFAM_PROFILES=(
    ["domain_A.hmm"]="PF00001"
    ["domain_B.hmm"]="PF00002"
    ["domain_C.hmm"]="PF00003"
    ["domain_D.hmm"]="PF00004"
    ["domain_E.hmm"]="PF00005"
)

for filename in "${!PFAM_PROFILES[@]}"; do
    accession="${PFAM_PROFILES[$filename]}"
    echo "  Downloading $accession → $filename"
    wget -q -O "$PROFILES_DIR/$filename" \
        "https://www.ebi.ac.uk/interpro/wwwapi/v1/entry/pfam/${accession}/?format=hmm"
done

echo "Downloading CDD profiles from NCBI..."

declare -A CDD_PROFILES=(
    ["domain_F.hmm"]="cd00001"
    ["domain_G.hmm"]="cd00002"
)

for filename in "${!CDD_PROFILES[@]}"; do
    accession="${CDD_PROFILES[$filename]}"
    echo "  Downloading $accession → $filename"
    wget -q "https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/hmm_raw/${accession}.hmm.gz"
    gunzip "${accession}.hmm.gz"
    mv "${accession}.hmm" "$PROFILES_DIR/$filename"
done

echo "Consolidating profiles into $OUTPUT_DB..."
cat "$PROFILES_DIR"/*.hmm > "$OUTPUT_DB"

echo "Indexing with hmmpress..."
hmmpress "$OUTPUT_DB"

echo "Done! HMM database ready at $OUTPUT_DB"
