#!/bin/bash

# This script extracts per-gene counts from featureCounts output files
# for multiple SRR samples listed in SRR_LIST.
# For each SRR it:
#   - reads the "counts_unique.txt" file produced by featureCounts,
#   - keeps only selected columns (gene_id, length, and counts columns),
#   - writes a clean per-gene table to OUTPUT_DIR as <SRR>_gene_counts_unique.tsv.
#
# No specific conda environment is required; standard Unix tools (head, tail, cut) are used.

SRR_LIST="input/srr_ids.txt"                     # List of SRR IDs (one per line)
COUNTS_ROOT="./GSE240542"                        # Root folder containing per-SRR featureCounts results
OUTPUT_DIR="./counts_unique"                     # Directory to store per-SRR gene count tables

mkdir -p "$OUTPUT_DIR"

while read -r SRR_ID; do
  [[ -z "$SRR_ID" ]] && continue

  COUNTS_FILE="${COUNTS_ROOT}/${SRR_ID}/featureCounts_2/counts_unique.txt"
  OUTPUT_FILE="${OUTPUT_DIR}/${SRR_ID}_gene_counts_unique.tsv"

  if [[ ! -f "$COUNTS_FILE" ]]; then
    echo "[$SRR_ID] counts_unique.txt not found, skipping"
    continue
  fi

  echo "[$SRR_ID] Processing $COUNTS_FILE"

  {
    # Header: second line (column names), select columns 1,6,7,8
    head -n 2 "$COUNTS_FILE" | tail -n 1 | cut -f 1,6,7,8
    # Data rows: from line 3 onward, same columns
    tail -n +3 "$COUNTS_FILE" | cut -f 1,6,7,8
  } > "$OUTPUT_FILE"

  echo "[$SRR_ID] Saved to $OUTPUT_FILE"

done < "$SRR_LIST"
