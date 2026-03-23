#!/usr/bin/env bash
# Author: Luiza Lima Galli
# Date: 2025
# Description: Runs HMMER hmmsearch in parallel against a set of proteome FASTA files, producing per-genome domain annotation tables (.domtab) and full reports (.txt).
# Usage: bash 02_run_hmmsearch.sh
# Dependencies: hmmer, GNU parallel
# License: MIT

set -euo pipefail

HMM_DB="hmm_profiles/toxin_panel.hmm"
PROTEOMES_DIR="proteomes"
RESULTS_DIR="results_hmm"
N_JOBS=16
LOG_FILE="hmm_search.log"

mkdir -p "$RESULTS_DIR"

echo "Starting parallel HMMER search with $N_JOBS jobs..."
echo "Input:  $PROTEOMES_DIR/*.fa"
echo "Output: $RESULTS_DIR/"

nohup parallel --jobs "$N_JOBS" \
    "hmmsearch \
        --domtblout $RESULTS_DIR/{/.}.domtab \
        $HMM_DB {} > $RESULTS_DIR/{/.}.txt" \
    ::: "$PROTEOMES_DIR"/*.fa > "$LOG_FILE" 2>&1 &

echo "Jobs submitted. Monitor progress with: tail -f $LOG_FILE"
echo "PID: $!"
