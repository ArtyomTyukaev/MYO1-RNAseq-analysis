#!/bin/bash

# This script performs parallel processing of .sra files using fasterq-dump and pigz.
# It:
#  - uses RNAseq_preprocessing conda enviroment
#  - reads SRR identifiers from a text file (one SRR per line),
#  - uses per-SRR locks to avoid concurrent processing of the same SRR,
#  - runs fasterq-dump to extract FASTQ files into a temporary directory (e.g. SSD),
#  - compresses the FASTQ files with pigz,
#  - moves the compressed FASTQ files back into the SRR directory,
#  - logs all operations.
#
# Requirements:
#  - SRA Toolkit (fasterq-dump) in PATH
#  - pigz in PATH
#  - GNU parallel in PATH
#  - conda/virtualenv/etc. should be activated manually before running, if needed.

LOGFILE="pipeline.log"
SRR_LIST="download_remained.txt"   # Text file with SRR IDs, one per line
MAX_PARALLEL=2                     # Number of SRRs processed in parallel
FQ_THREADS=10                      # Threads for fasterq-dump
PIGZ_THREADS=4                     # Threads for pigz
LOCK_DIR=".locks"
SSD_TMP="./fasterq_tmp"            # Temporary directory for fasterq-dump output

# Initialization
mkdir -p "$LOCK_DIR"
mkdir -p "$SSD_TMP"
> "$LOGFILE"

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"
}

process_srr() {
    local SRR_ID="$1"
    local LOCK_FILE="$LOCK_DIR/$SRR_ID.lock"
    local NFS_DIR="$PWD/$SRR_ID"
    local TMP_DIR="$SSD_TMP/$SRR_ID"

    # Check for existing lock
    if [[ -f "$LOCK_FILE" ]]; then
        log "[$SRR_ID] Already being processed — skipping"
        return 0
    fi

    # Create lock
    if ! touch "$LOCK_FILE"; then
        log "[$SRR_ID] Failed to create lock"
        return 1
    fi

    log "[$SRR_ID] Processing started"

    if [[ ! -d "$NFS_DIR" ]]; then
        log "[$SRR_ID] Directory not found — skipping"
        rm -f "$LOCK_FILE"
        return 1
    fi

    if [[ ! -f "$NFS_DIR/${SRR_ID}.sra" ]]; then
        log "[$SRR_ID] .sra file missing"
        rm -f "$LOCK_FILE"
        return 1
    fi

    if [[ -f "$NFS_DIR/${SRR_ID}_1.fastq.gz" && -f "$NFS_DIR/${SRR_ID}_2.fastq.gz" ]]; then
        log "[$SRR_ID] Already processed — skipping"
        rm -f "$LOCK_FILE"
        return 0
    fi

    mkdir -p "$TMP_DIR"

    # Run fasterq-dump to temporary directory
    log "[$SRR_ID] Running fasterq-dump to temporary directory $TMP_DIR"
    START_TIME=$(date +%s)
    if ! fasterq-dump --split-files --threads "$FQ_THREADS" \
        --outdir "$TMP_DIR" --temp "$TMP_DIR/tmp" \
        "$NFS_DIR/${SRR_ID}.sra" >> "$LOGFILE" 2>&1; then
        log "[$SRR_ID] fasterq-dump error"
        rm -rf "$TMP_DIR"
        rm -f "$LOCK_FILE"
        return 1
    fi
    END_TIME=$(date +%s)
    log "[$SRR_ID] fasterq-dump completed in $((END_TIME - START_TIME)) seconds"

    # Compress with pigz and move back
    for fq in "$TMP_DIR"/${SRR_ID}_*.fastq; do
        if [[ -f "$fq" ]]; then
            log "[$SRR_ID] Compressing $fq with pigz"
            if pigz -p "$PIGZ_THREADS" "$fq" >> "$LOGFILE" 2>&1; then
                mv "${fq}.gz" "$NFS_DIR/"
            else
                log "[$SRR_ID] Compression failed for $fq"
            fi
        fi
    done

    rm -rf "$TMP_DIR"
    rm -f "$LOCK_FILE"
    log "[$SRR_ID] Processing completed successfully"
    return 0
}

export -f process_srr log
export LOGFILE FQ_THREADS PIGZ_THREADS LOCK_DIR SSD_TMP

log "Pipeline started"
log "Using total CPU cores: $((MAX_PARALLEL * FQ_THREADS)) ($MAX_PARALLEL parallel SRRs × $FQ_THREADS threads each)"

if ! command -v parallel >/dev/null 2>&1; then
    log "GNU parallel is not installed"
    exit 1
fi

if [[ ! -f "$SRR_LIST" ]]; then
    log "Input list $SRR_LIST not found"
    exit 1
fi

cat "$SRR_LIST" | parallel -j "$MAX_PARALLEL" process_srr {}

log "All tasks completed"

