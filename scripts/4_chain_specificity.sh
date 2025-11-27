#!/bin/bash

# This script uses RSeQC's infer_experiment.py to determine library strandedness
# for each SRR ID listed in SRR_LIST. For each sample it:
#   - runs infer_experiment.py with a BED annotation and sorted BAM file,
#   - parses the fractions of reads explained by sense/antisense models,
#   - classifies the library as stranded, reversely stranded, or unstranded,
#   - writes separate SRR lists per category into OUTPUT_DIR,
#   - logs detailed output to LOG_FILE.
#
# Conda environment required:
#   define_chain_specificity.yml  (must contain infer_experiment.py, and 'bc' should be available in PATH)

if [[ -z "${CONDA_DEFAULT_ENV:-}" ]] || [[ "$CONDA_DEFAULT_ENV" != "rseqc" ]]; then
  echo "Activating conda environment 'rseqc'..."
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate rseqc || { echo "Error: conda environment 'rseqc' not found"; exit 1; }
fi

SRR_LIST="input/srr_ids.txt"                     # List of SRR IDs (one per line)
BED_FILE="./refs/gencode.v48.basic.annotation.bed"  # BED annotation used by infer_experiment.py
LOG_FILE="logs/chain.log"
OUTPUT_DIR="output"
BAM_ROOT="./GSE240542"                          # Root directory containing per-SRR folders
SAMPLE_SIZE=500000
UNSTRANDED_EPS=0.1

mkdir -p "$OUTPUT_DIR"
mkdir -p "$(dirname "$LOG_FILE")"
# Clear log file
> "$LOG_FILE"

# Counters
sense_count=0
antisense_count=0
unstranded_count=0

while read -r SRR; do
  [[ -z "$SRR" ]] && continue

  BAM="${BAM_ROOT}/${SRR}/STAR_results_2/${SRR}_Aligned.sortedByCoord.out.bam"
  if [[ ! -f "$BAM" ]]; then
    echo "[$SRR] BAM file not found: $BAM" >> "$LOG_FILE"
    continue
  fi

  echo "Processing $SRR" >> "$LOG_FILE"

  # Capture command output
  output=$(infer_experiment.py -r "$BED_FILE" -i "$BAM" -s "$SAMPLE_SIZE" 2>/dev/null)

  # Extract lines of interest
  pair_info=$(echo "$output" | grep "This is")
  failed=$(echo "$output" | grep "Fraction of reads failed to determine")
  sense=$(echo "$output" | grep 'Fraction of reads explained by "1++,1--,2+-,2-+')
  antisense=$(echo "$output" | grep 'Fraction of reads explained by "1+-,1-+,2++,2--"')

  # Extract numeric values
  sense_val=$(echo "$sense" | grep -oE "[0-9]+\.[0-9]+")
  antisense_val=$(echo "$antisense" | grep -oE "[0-9]+\.[0-9]+")

  # Log raw lines
  echo "$pair_info" >> "$LOG_FILE"
  echo "$failed" >> "$LOG_FILE"
  echo "$sense" >> "$LOG_FILE"
  echo "$antisense" >> "$LOG_FILE"

  # Determine library type
  diff=$(echo "$sense_val - $antisense_val" | bc -l)
  abs_diff=$(echo "${diff#-}" | bc -l)

  if (( $(echo "$abs_diff < $UNSTRANDED_EPS" | bc -l) )); then
    echo "$SRR - unstranded" >> "$LOG_FILE"
    ((unstranded_count++))
  elif (( $(echo "$sense_val > $antisense_val" | bc -l) )); then
    echo "$SRR - stranded (sense strand)" >> "$LOG_FILE"
    ((sense_count++))
  else
    echo "$SRR - reversely stranded (antisense)" >> "$LOG_FILE"
    ((antisense_count++))
  fi

  echo "" >> "$LOG_FILE"

done < "$SRR_LIST"

# --- POST-PROCESSING: build SRR lists per category ---

declare -a SRR_UNSTRANDED=()
declare -a SRR_STRANDED=()
declare -a SRR_REVERSELYSTRANDED=()

# Parse SRR IDs and categories from the log file
while read -r line; do
  if [[ "$line" =~ ^SRR[0-9]+ ]]; then
    SRR=$(echo "$line" | awk '{print $1}')
    if   [[ "$line" == *"unstranded"* ]]; then
      SRR_UNSTRANDED+=("$SRR")
    elif [[ "$line" == *"sense strand"* ]]; then
      SRR_STRANDED+=("$SRR")
    elif [[ "$line" == *"reversely stranded"* ]]; then
      SRR_REVERSELYSTRANDED+=("$SRR")
    fi
  fi
done < "$LOG_FILE"

# Create files only if the corresponding list is non-empty
if ((${#SRR_UNSTRANDED[@]})); then
  printf "%s\n" "${SRR_UNSTRANDED[@]}" > "${OUTPUT_DIR}/srr_unstranded.txt"
fi

if ((${#SRR_STRANDED[@]})); then
  printf "%s\n" "${SRR_STRANDED[@]}" > "${OUTPUT_DIR}/srr_stranded.txt"
fi

if ((${#SRR_REVERSELYSTRANDED[@]})); then
  printf "%s\n" "${SRR_REVERSELYSTRANDED[@]}" > "${OUTPUT_DIR}/srr_reverselystranded.txt"
fi

echo "SRR lists saved to directory: ${OUTPUT_DIR} (unstranded threshold = ${UNSTRANDED_EPS})" | tee -a "$LOG_FILE"

# Final statistics
{
  echo "Strandedness summary:"
  echo "Stranded (sense strand): $sense_count"
  echo "Reversely stranded (antisense): $antisense_count"
  echo "Unstranded: $unstranded_count"
} >> "$LOG_FILE"
