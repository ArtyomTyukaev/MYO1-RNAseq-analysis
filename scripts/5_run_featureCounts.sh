#!/bin/bash

# This script runs featureCounts for multiple SRR samples listed in SRR_LIST.
# For each SRR it:
#   - copies the sorted BAM to an SSD-based temporary directory,
#   - runs featureCounts twice:
#       1) unique counts only (for DE / gene-level counts)
#       2) with multimapping reads and fractional assignment (for TPM / QC),
#   - saves outputs into per-SRR featureCounts directories,
#   - uses file locks to avoid concurrent processing of the same SRR,
#   - logs all operations to LOG_FILE.
#
# Conda environment required:
#   RNAseq_preprocessing  (must contain subread/featureCounts)

# Activate conda environment if needed
if [[ -z "${CONDA_DEFAULT_ENV:-}" ]] || [[ "$CONDA_DEFAULT_ENV" != "RNAseq" ]]; then
  echo "Activating conda environment: RNAseq..."
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate RNAseq || { echo "Error: conda environment 'RNAseq' not found"; exit 1; }
fi

### S_PARAM
# 0 = unstranded
# 1 = stranded
# 2 = reversely stranded
S_PARAM=2

# Paths
GTF_FILE="./refs/gencode.v48.basic.annotation.gtf"
SRR_LIST="input/srr_remaining.txt"
LOG_FILE="./logs/featureCounts_r.log"
SSD_TMP="./tmp_featurecounts"
MAX_PARALLEL=4    # Number of SRRs processed in parallel
FC_THREADS=5      # Threads per featureCounts run
LOCK_DIR=".fc_locks"

# Checks
if [[ ! -f "$SRR_LIST" ]]; then
  echo "Error: SRR list file not found: $SRR_LIST" | tee -a "$LOG_FILE"
  exit 1
fi
if [[ ! -f "$GTF_FILE" ]]; then
  echo "Error: GTF file not found at: $GTF_FILE" | tee -a "$LOG_FILE"
  exit 1
fi

# Initialization
mkdir -p "$SSD_TMP"
mkdir -p "$LOCK_DIR"
mkdir -p "$(dirname "$LOG_FILE")"
> "$LOG_FILE"

echo "=== featureCounts run started ===" | tee -a "$LOG_FILE"
echo "Start date: $(date)" | tee -a "$LOG_FILE"
echo "Parameters: --fracOverlap 0.1, -T $FC_THREADS, -s $S_PARAM" | tee -a "$LOG_FILE"

log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

process_srr() {
  local SRR_ID="$1"
  local LOCK_FILE="$LOCK_DIR/$SRR_ID.lock"
  local BAM_FILE="./GSE240542/${SRR_ID}/STAR_results_2/${SRR_ID}_Aligned.sortedByCoord.out.bam"
  local RESULT_DIR="./GSE240542/${SRR_ID}/featureCounts_2"
  local TMP_DIR="$SSD_TMP/$SRR_ID"

  if [[ -f "$LOCK_FILE" ]]; then
    log "[$SRR_ID] Already being processed by another process — skipping"
    return 0
  fi
  if ! touch "$LOCK_FILE"; then
    log "[$SRR_ID] Failed to create lock file"
    return 1
  fi

  if [[ ! -f "$BAM_FILE" ]]; then
    log "[$SRR_ID] BAM file not found"
    rm -f "$LOCK_FILE"
    return 1
  fi

  mkdir -p "$TMP_DIR"
  if ! cp "$BAM_FILE" "$TMP_DIR/"; then
    log "[$SRR_ID] Failed to copy BAM to SSD temporary directory"
    rm -rf "$TMP_DIR"
    rm -f "$LOCK_FILE"
    return 1
  fi

  mkdir -p "$RESULT_DIR"

  log "[$SRR_ID] Running featureCounts"
  START_TIME=$(date +%s)

  # Run 1: unique counts only
  featureCounts \
    -a "$GTF_FILE" \
    -o "$TMP_DIR/counts_unique.txt" \
    -g gene_id \
    --extraAttributes gene_name \
    --countReadPairs \
    -t exon \
    -p -B -C \
    -Q 10 \
    --fracOverlap 0.1 \
    -T "$FC_THREADS" \
    -s "$S_PARAM" \
    "$TMP_DIR/$(basename "$BAM_FILE")" >> "$LOG_FILE" 2>&1
  FC_EXIT_CODE_1=$?

  # Run 2: with multimapping reads, fractional assignment
  featureCounts \
    -a "$GTF_FILE" \
    -o "$TMP_DIR/counts_fraction.txt" \
    -g gene_id \
    --extraAttributes gene_name \
    --countReadPairs \
    -t exon \
    -p -B -C \
    -M --fraction \
    -Q 0 \
    --fracOverlap 0.1 \
    -T "$FC_THREADS" \
    -s "$S_PARAM" \
    "$TMP_DIR/$(basename "$BAM_FILE")" >> "$LOG_FILE" 2>&1
  FC_EXIT_CODE_2=$?

  END_TIME=$(date +%s)

  if [[ "$FC_EXIT_CODE_1" -eq 0 && "$FC_EXIT_CODE_2" -eq 0 ]]; then
    log "[$SRR_ID] Both featureCounts runs finished successfully in $((END_TIME - START_TIME)) seconds"
    mv "$TMP_DIR/counts_unique.txt" "$RESULT_DIR/counts_unique.txt"
    mv "$TMP_DIR/counts_unique.txt.summary" "$RESULT_DIR/${SRR_ID}_unique.summary"
    mv "$TMP_DIR/counts_fraction.txt" "$RESULT_DIR/counts_fraction.txt"
    mv "$TMP_DIR/counts_fraction.txt.summary" "$RESULT_DIR/${SRR_ID}_fraction.summary"
  else
    log "[$SRR_ID] Error in one of the featureCounts runs"
  fi
  
  rm -rf "$TMP_DIR"
  rm -f "$LOCK_FILE"
}

export -f process_srr log
export GTF_FILE LOG_FILE SSD_TMP FC_THREADS LOCK_DIR S_PARAM

log "Parallel processing started"
cat "$SRR_LIST" | parallel -j "$MAX_PARALLEL" process_srr {}
log "Processing finished"

# Notes on strand-specific counting (-s / isStrandSpecific):
# A single integer value (applied to all input files) or a string of comma-separated
# values (applied to each corresponding input file) should be provided.
# Possible values include: 0 (unstranded), 1 (stranded) and 2 (reversely stranded).
# Default value is 0 (i.e. unstranded read counting). For paired-end reads, strand of
# the first read is taken as the strand of the whole fragment. FLAG field is used to
# tell if a read is first or second in a pair.

