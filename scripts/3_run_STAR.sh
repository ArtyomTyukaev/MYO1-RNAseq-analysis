#!/usr/bin/env bash
set -euo pipefail

# This script runs STAR alignment for multiple SRR samples listed in SRR_LIST.
# For each SRR:
#   - uses trimmed FASTQ files from SRR*/trimmed_fastq/
#   - copies them to a temporary directory (prefer SSD, fallback to HDD)
#   - runs STAR with a given genome index and GTF annotation
#   - writes outputs into SRR*/STAR_results_2/
#   - sorts the resulting BAM with samtools
#   - logs all operations to LOGFILE
#
# Conda environment required:
#   RNAseq_preprocessing  (must contain STAR and samtools)

### ============ Settings ============

# Working directory containing SRR* folders
WORKDIR="./data/GSE240542"
cd "$WORKDIR" || { echo "Failed to cd into $WORKDIR"; exit 1; }

# STAR genome index and annotation
GENOME_DIR="./refs/star_index"
GTF_FILE="./refs/annotation.gtf"

# List of SRR IDs to process (one per line)
SRR_LIST="./input/srr_ids_remaining.txt"

# Logging and temporary directories
LOGFILE="./logs/STAR.log"
SSD_TMP="./tmp_star_ssd"
HDD_TMP="./tmp_star_hdd"

THREADS=15
REQUIRED_RAM_GB=50

### ============ Conda activation (if needed) ============
if [[ -z "${CONDA_DEFAULT_ENV:-}" || "$CONDA_DEFAULT_ENV" != "RNAseq" ]]; then
  echo "Activating conda environment: RNAseq..."
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate RNAseq || { echo "Error: environment 'RNAseq' not found"; exit 1; }
fi

mkdir -p "$SSD_TMP" "$HDD_TMP"
mkdir -p "$(dirname "$LOGFILE")"
> "$LOGFILE"

log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"
}

wait_for_memory() {
  while true; do
    FREE_MEM=$(free -g | awk '/Mem:/ {print $7}')
    if [[ "$FREE_MEM" -ge "$REQUIRED_RAM_GB" ]]; then
      return
    fi
    log "Not enough free memory (${FREE_MEM} GB available, need >= ${REQUIRED_RAM_GB} GB). Waiting 5 minutes..."
    sleep 300
  done
}

run_star_with_tmpdir() {
  local TMP_BASE="$1"
  local SRR_ID="$2"
  local R1="$3"
  local R2="$4"
  local TMP_DIR="$TMP_BASE/$SRR_ID"

  mkdir -p "$TMP_DIR"
  cp "$R1" "$TMP_DIR/" && cp "$R2" "$TMP_DIR/" || { rm -rf "$TMP_DIR"; return 1; }

  local TMP_R1="$TMP_DIR/${SRR_ID}_1.trim.fastq.gz"
  local TMP_R2="$TMP_DIR/${SRR_ID}_2.trim.fastq.gz"
  local OUT_DIR="$WORKDIR/$SRR_ID/STAR_results_2"
  local SAMPLE_LOG="$OUT_DIR/${SRR_ID}.log"
  mkdir -p "$OUT_DIR"

  wait_for_memory
  log "[$SRR_ID] Running STAR with temporary base directory: $TMP_BASE"

  STAR \
    --runThreadN "$THREADS" \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$TMP_R1" "$TMP_R2" \
    --readFilesCommand zcat \
    --outSAMtype BAM Unsorted \
    --sjdbGTFfile "$GTF_FILE" \
    --outFileNamePrefix "$OUT_DIR/${SRR_ID}_" \
    >> "$SAMPLE_LOG" 2>&1

  local status=$?
  rm -rf "$TMP_DIR"

  if [[ $status -ne 0 ]]; then
    log "[$SRR_ID] STAR exited with non-zero status (likely memory-related). Skipping this sample."
    return 1
  fi
}

for SRR_ID in $(cat "$SRR_LIST"); do
  log "[$SRR_ID] Processing started"

  SRR_DIR="$WORKDIR/$SRR_ID"
  TRIM_DIR="$SRR_DIR/trimmed_fastq"
  R1="$TRIM_DIR/${SRR_ID}_1.trim.fastq.gz"
  R2="$TRIM_DIR/${SRR_ID}_2.trim.fastq.gz"

  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    log "[$SRR_ID] Trimmed FASTQ files not found"
    continue
  fi

  if run_star_with_tmpdir "$SSD_TMP" "$SRR_ID" "$R1" "$R2"; then
    log "[$SRR_ID] STAR completed successfully on SSD temp"
  elif run_star_with_tmpdir "$HDD_TMP" "$SRR_ID" "$R1" "$R2"; then
    log "[$SRR_ID] STAR completed successfully on HDD temp"
  else
    log "[$SRR_ID] STAR failed in both SSD and HDD temp locations — skipping this sample"
    continue
  fi

  # BAM sorting
  wait_for_memory

  UNSORTED_BAM="$SRR_DIR/STAR_results_2/${SRR_ID}_Aligned.out.bam"
  SORTED_BAM="$SRR_DIR/STAR_results_2/${SRR_ID}_Aligned.sortedByCoord.out.bam"

  if [[ -f "$UNSORTED_BAM" ]]; then
    log "[$SRR_ID] Sorting BAM with samtools..."

    samtools sort -@ 8 -o "$SORTED_BAM" "$UNSORTED_BAM" >> "$LOGFILE" 2>&1
    if [[ $? -eq 0 ]]; then
      log "[$SRR_ID] BAM successfully sorted"
      rm "$UNSORTED_BAM"
    else
      log "[$SRR_ID] BAM sorting failed (likely due to insufficient memory). Skipping this sample"
      continue
    fi
  else
    log "[$SRR_ID] Unsorted BAM file not found after STAR"
  fi

done

log "All SRR samples processed"

