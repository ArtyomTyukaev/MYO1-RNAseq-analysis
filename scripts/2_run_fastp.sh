#!/usr/bin/env bash
set -euo pipefail

# This script runs fastp trimming for all SRR* folders inside WORKDIR.
# It expects each SRR directory to contain paired-end FASTQ files:
#   SRRXXXXXXX_1.fastq.gz
#   SRRXXXXXXX_2.fastq.gz
#
# Output:
#   - trimmed FASTQs written into SRR*/trimmed_fastq/
#   - fastp HTML/JSON reports written into REPORTS_DIR
#   - a combined MultiQC report (optional)
#
# Conda environment required:
#   RNAseq_preprocessing  (must contain fastp and optionally multiqc)

### ============ Settings ============

# Location of SRR* directories
WORKDIR="./data"

# Output locations
REPORTS_DIR="./reports/fastp"
LOGFILE="./logs/fastp_run.log"

# Threads for fastp
FP_THREADS=8

### ============ Conda activation (if needed) ============
if [[ -z "${CONDA_DEFAULT_ENV:-}" ]] || [[ "$CONDA_DEFAULT_ENV" != "RNAseq" ]]; then
  echo "Activating conda environment: RNAseq..."
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate RNAseq || { echo "Error: environment 'RNAseq' not found"; exit 1; }
fi

mkdir -p "$REPORTS_DIR"
mkdir -p "$(dirname "$LOGFILE")"
: > "$LOGFILE"

log(){ echo "[$(date +'%F %T')] $*" | tee -a "$LOGFILE"; }

### ============ Process all SRR directories ============
shopt -s nullglob
for SRR_DIR in "$WORKDIR"/SRR*; do
  [[ -d "$SRR_DIR" ]] || continue
  SRR_ID="$(basename "$SRR_DIR")"

  R1="$SRR_DIR/${SRR_ID}_1.fastq.gz"
  R2="$SRR_DIR/${SRR_ID}_2.fastq.gz"

  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    log "[$SRR_ID] FASTQ files not found — skipping"
    continue
  fi

  # Trimmed output inside each SRR directory
  TRIM_DIR="$SRR_DIR/trimmed_fastq"
  mkdir -p "$TRIM_DIR"
  O1="$TRIM_DIR/${SRR_ID}_1.trim.fastq.gz"
  O2="$TRIM_DIR/${SRR_ID}_2.trim.fastq.gz"

  # fastp reports in a shared reports directory
  HTML="$REPORTS_DIR/${SRR_ID}.fastp.html"
  JSON="$REPORTS_DIR/${SRR_ID}.fastp.json"

  # Skip if everything already exists
  if [[ -s "$O1" && -s "$O2" && -s "$HTML" && -s "$JSON" ]]; then
    log "[$SRR_ID] Trimmed files and reports already exist — skipping"
    continue
  fi

  log "[$SRR_ID] fastp started"
  fastp -w "$FP_THREADS" \
    -i "$R1" -I "$R2" \
    -o "$O1" -O "$O2" \
    --detect_adapter_for_pe \
    --cut_front \
    --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 22 \
    --length_required 30 \
    --html "$HTML" \
    --json "$JSON" \
    >> "$LOGFILE" 2>&1

  status=$?
  if [[ $status -eq 0 ]]; then
    log "[$SRR_ID] fastp completed"
  else
    log "[$SRR_ID] fastp failed with exit code $status"
  fi
done

### ============ Optional MultiQC summary ============
if command -v multiqc >/dev/null 2>&1; then
  log "Generating MultiQC summary from fastp reports..."
  multiqc "$REPORTS_DIR" -o "$(dirname "$REPORTS_DIR")/fastp_multiqc" >> "$LOGFILE" 2>&1 || true
fi

log "Done"

