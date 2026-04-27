#!/usr/bin/env bash
# =============================================================================
# Winnowmap2 ONT Alignment Pipeline with Methylation Tag Preservation
# =============================================================================
# Version: 1.1.0
# Author:  Venkat
# License: MIT
#
# Usage:
#   bash winnowmap_align.sh -i <ubam_or_dir> -r <reference.fa> -o <output_dir> [options]
#
# Required:
#   -i  Input uBAM file OR directory containing multiple uBAM chunks
#   -r  Reference FASTA (e.g. hg38.fa — must be samtools faidx indexed)
#   -o  Output directory (created if it doesn't exist)
#
# Optional:
#   -s  Sample name           [default: derived from uBAM filename]
#   -t  Threads               [default: 16]
#   -k  k-mer size            [default: 15]
#   -w  Meryl database dir    [default: <ref_dir>/<ref_base>.meryl]
#   -f  High-freq k-mer file  [default: <ref_dir>/<ref_base>_repetitive_k<k>.txt]
#   -m  Skip meryl steps if k-mer file already exists [flag]
#   -c  Cleanup intermediates after completion [flag]
#   -h  Show this help message
#
# Examples:
#   # Single uBAM
#   bash winnowmap_align.sh \
#     -i /data/ubam/sample.bam \
#     -r /data/reference/hg38.fa \
#     -o /data/output \
#     -s MySample -t 16 -c
#
#   # Directory of uBAM chunks (auto-merge)
#   bash winnowmap_align.sh \
#     -i /data/ubam_chunks/ \
#     -r /data/reference/hg38.fa \
#     -o /data/output \
#     -s MySample -t 16 -c
#
#   # Reuse existing k-mer file (skip meryl)
#   bash winnowmap_align.sh \
#     -i /data/ubam/sample.bam \
#     -r /data/reference/hg38.fa \
#     -o /data/output \
#     -f /data/reference/hg38_repetitive_k15.txt \
#     -m -c
# =============================================================================

set -euo pipefail

# ── Version ───────────────────────────────────────────────────────────────────
VERSION="1.1.0"

# ── Defaults ──────────────────────────────────────────────────────────────────
THREADS=16
KMER=15
SKIP_MERYL=false
CLEANUP=false
SAMPLE=""
MERYL_DIR=""
KMER_FILE=""
INPUT=""

# ── Colors ────────────────────────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

log_info()    { echo -e "${BLUE}[INFO]${NC}  $(date '+%H:%M:%S') $*"; }
log_success() { echo -e "${GREEN}[OK]${NC}    $(date '+%H:%M:%S') $*"; }
log_warn()    { echo -e "${YELLOW}[WARN]${NC}  $(date '+%H:%M:%S') $*"; }
log_error()   { echo -e "${RED}[ERROR]${NC} $(date '+%H:%M:%S') $*"; exit 1; }
log_step()    { echo -e "${CYAN}[STEP]${NC}  $(date '+%H:%M:%S') $*"; }

# ── Help ──────────────────────────────────────────────────────────────────────
usage() {
    sed -n '2,50p' "$0" | sed 's/^# \{0,1\}//'
    exit 0
}

# ── Parse arguments ───────────────────────────────────────────────────────────
while getopts ":i:r:o:s:t:k:w:f:mch" opt; do
    case $opt in
        i) INPUT="$OPTARG"    ;;
        r) REF="$OPTARG"      ;;
        o) OUTDIR="$OPTARG"   ;;
        s) SAMPLE="$OPTARG"   ;;
        t) THREADS="$OPTARG"  ;;
        k) KMER="$OPTARG"     ;;
        w) MERYL_DIR="$OPTARG";;
        f) KMER_FILE="$OPTARG";;
        m) SKIP_MERYL=true    ;;
        c) CLEANUP=true       ;;
        h) usage              ;;
        :) log_error "Option -$OPTARG requires an argument." ;;
        \?) log_error "Unknown option: -$OPTARG" ;;
    esac
done

# ── Validate required arguments ───────────────────────────────────────────────
[[ -z "${INPUT:-}" ]] && log_error "Input uBAM or directory (-i) is required."
[[ -z "${REF:-}"   ]] && log_error "Reference FASTA (-r) is required."
[[ -z "${OUTDIR:-}"]] && log_error "Output directory (-o) is required."
[[ ! -f "$REF"     ]] && log_error "Reference FASTA not found: $REF"

# ── Check dependencies ────────────────────────────────────────────────────────
log_info "Checking dependencies..."
for tool in winnowmap meryl samtools python3; do
    if command -v "$tool" &>/dev/null; then
        log_success "$tool: $(command -v $tool)"
    else
        log_error "$tool not found. Please install and re-run."
    fi
done
python3 -c "import pysam" 2>/dev/null || log_error "pysam not found. Install: pip install pysam"
log_success "pysam: OK"

# ── Derive defaults ───────────────────────────────────────────────────────────
REF_DIR=$(dirname "$REF")
REF_BASE=$(basename "$REF")
REF_BASE="${REF_BASE%.fa}"
REF_BASE="${REF_BASE%.fasta}"

[[ -z "$MERYL_DIR" ]] && MERYL_DIR="${REF_DIR}/${REF_BASE}.meryl"
[[ -z "$KMER_FILE" ]] && KMER_FILE="${REF_DIR}/${REF_BASE}_repetitive_k${KMER}.txt"

mkdir -p "$OUTDIR"

# ── Handle input: single BAM or directory of BAM chunks ──────────────────────
UBAM=""
MERGED_CHUNKS=false

if [[ -f "$INPUT" ]]; then
    # Single BAM file
    UBAM="$INPUT"
    [[ -z "$SAMPLE" ]] && SAMPLE=$(basename "$UBAM" .bam)
    log_info "Input mode: single uBAM → $UBAM"

elif [[ -d "$INPUT" ]]; then
    # Directory of BAM chunks — auto-merge
    log_info "Input mode: directory of uBAM chunks → $INPUT"
    mapfile -t BAM_CHUNKS < <(find "$INPUT" -name "*.bam" | sort)
    CHUNK_COUNT=${#BAM_CHUNKS[@]}

    [[ $CHUNK_COUNT -eq 0 ]] && log_error "No BAM files found in directory: $INPUT"
    log_info "Found ${CHUNK_COUNT} BAM chunk(s)"

    if [[ $CHUNK_COUNT -eq 1 ]]; then
        UBAM="${BAM_CHUNKS[0]}"
        [[ -z "$SAMPLE" ]] && SAMPLE=$(basename "$UBAM" .bam)
        log_warn "Only 1 chunk found — using directly without merge"
    else
        [[ -z "$SAMPLE" ]] && SAMPLE=$(basename "$INPUT" | sed 's/[^a-zA-Z0-9_-]/_/g')
        UBAM="${OUTDIR}/${SAMPLE}_merged_ubam.bam"

        if [[ -f "$UBAM" ]]; then
            log_warn "Merged uBAM already exists — skipping merge: $UBAM"
        else
            log_step "Merging ${CHUNK_COUNT} uBAM chunks..."
            MERGE_START=$(date +%s)

            CHUNK_LIST="${OUTDIR}/${SAMPLE}_chunk_list.txt"
            printf '%s\n' "${BAM_CHUNKS[@]}" > "$CHUNK_LIST"

            samtools merge -@ "$THREADS" -f -b "$CHUNK_LIST" "$UBAM"
            samtools index "$UBAM"

            MERGE_END=$(date +%s)
            MERGED_READS=$(samtools flagstat "$UBAM" | awk 'NR==1{print $1}')
            log_success "Merged ${CHUNK_COUNT} chunks → ${MERGED_READS} total reads in $((MERGE_END - MERGE_START))s"
            MERGED_CHUNKS=true
        fi
    fi
else
    log_error "Input (-i) must be a BAM file or a directory: $INPUT"
fi

# ── Print run summary ─────────────────────────────────────────────────────────
echo ""
echo -e "${BLUE}======================================================${NC}"
echo -e "${BLUE}   Winnowmap2 ONT Alignment Pipeline v${VERSION}${NC}"
echo -e "${BLUE}======================================================${NC}"
echo "  Input        : $INPUT"
echo "  uBAM used    : $UBAM"
echo "  Reference    : $REF"
echo "  Output dir   : $OUTDIR"
echo "  Sample name  : $SAMPLE"
echo "  Threads      : $THREADS"
echo "  k-mer size   : $KMER"
echo "  Meryl dir    : $MERYL_DIR"
echo "  k-mer file   : $KMER_FILE"
echo "  Skip meryl   : $SKIP_MERYL"
echo "  Cleanup      : $CLEANUP"
echo -e "${BLUE}======================================================${NC}"
echo ""

PIPELINE_START=$(date +%s)

# ── STEP 1 & 2: Meryl k-mer counting and extraction ──────────────────────────
if [[ "$SKIP_MERYL" == true ]] && [[ -f "$KMER_FILE" ]]; then
    log_warn "Skipping meryl — reusing existing k-mer file: $KMER_FILE"
elif [[ -f "$KMER_FILE" ]]; then
    log_warn "k-mer file already exists — skipping meryl: $KMER_FILE"
else
    log_step "STEP 1/6: Counting ${KMER}-mers from reference with meryl..."
    STEP_START=$(date +%s)
    meryl count k="$KMER" "$REF" output "$MERYL_DIR"
    STEP_END=$(date +%s)
    log_success "meryl count complete in $((STEP_END - STEP_START))s"

    log_step "STEP 2/6: Extracting high-frequency k-mers (distinct > 0.9998)..."
    STEP_START=$(date +%s)
    meryl print greater-than distinct=0.9998 "$MERYL_DIR" > "$KMER_FILE"
    KMER_COUNT=$(wc -l < "$KMER_FILE")
    STEP_END=$(date +%s)
    log_success "Extracted ${KMER_COUNT} repetitive k-mers in $((STEP_END - STEP_START))s → $KMER_FILE"
fi

# ── STEP 3: uBAM → FASTQ ─────────────────────────────────────────────────────
FASTQ="${OUTDIR}/${SAMPLE}_reads.fastq"
log_step "STEP 3/6: Converting uBAM to FASTQ (preserving MM/ML methylation tags)..."
STEP_START=$(date +%s)

samtools fastq -T MM,ML "$UBAM" > "$FASTQ"

READ_COUNT=$(grep -c "^@" "$FASTQ" || true)
STEP_END=$(date +%s)
log_success "Extracted ${READ_COUNT} reads in $((STEP_END - STEP_START))s → $FASTQ"

# ── STEP 4: Winnowmap2 alignment ──────────────────────────────────────────────
SAM="${OUTDIR}/${SAMPLE}_aligned.sam"
log_step "STEP 4/6: Running Winnowmap2 alignment (k=${KMER}, threads=${THREADS})..."
STEP_START=$(date +%s)

winnowmap \
    -W "$KMER_FILE" \
    -ax map-ont \
    -k "$KMER" \
    -t "$THREADS" \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ONT" \
    --MD \
    --eqx \
    -o "$SAM" \
    "$REF" \
    "$FASTQ"

STEP_END=$(date +%s)
log_success "Alignment complete in $((STEP_END - STEP_START))s → $SAM"

# ── STEP 5: Sort and index aligned BAM ───────────────────────────────────────
ALIGNED_SORTED="${OUTDIR}/${SAMPLE}_aligned_sorted.bam"
log_step "STEP 5/6: Sorting and indexing aligned BAM..."
STEP_START=$(date +%s)

samtools sort -@ "$THREADS" -o "$ALIGNED_SORTED" "$SAM"
samtools index "$ALIGNED_SORTED"

STEP_END=$(date +%s)
log_success "Sorted BAM ready in $((STEP_END - STEP_START))s → $ALIGNED_SORTED"

# ── STEP 6: Methylation tag transplant ───────────────────────────────────────
log_step "STEP 6/6: Transplanting MM/ML methylation tags from uBAM..."
STEP_START=$(date +%s)

UBAM_NAMESORTED="${OUTDIR}/${SAMPLE}_ubam_namesorted.bam"
ALIGNED_NAMESORTED="${OUTDIR}/${SAMPLE}_aligned_namesorted.bam"
METH_UNSORTED="${OUTDIR}/${SAMPLE}_with_meth_unsorted.bam"
FINAL_BAM="${OUTDIR}/${SAMPLE}_aligned_meth_sorted.bam"

log_info "  Name-sorting uBAM and aligned BAM..."
samtools sort -@ "$THREADS" -n "$UBAM"           -o "$UBAM_NAMESORTED"
samtools sort -@ "$THREADS" -n "$ALIGNED_SORTED" -o "$ALIGNED_NAMESORTED"

python3 - <<PYEOF
import pysam
import sys

ubam_path    = "${UBAM_NAMESORTED}"
aligned_path = "${ALIGNED_NAMESORTED}"
out_path     = "${METH_UNSORTED}"

print("  Loading MM/ML tags from uBAM...")
ubam = pysam.AlignmentFile(ubam_path, "rb", check_sq=False)
meth_tags = {}
for read in ubam:
    tags = {}
    if read.has_tag("MM"): tags["MM"] = read.get_tag("MM")
    if read.has_tag("ML"): tags["ML"] = read.get_tag("ML")
    if tags:
        meth_tags[read.query_name] = tags
ubam.close()
print(f"  Loaded tags for {len(meth_tags):,} reads")

aligned = pysam.AlignmentFile(aligned_path, "rb")
out     = pysam.AlignmentFile(out_path, "wb", header=aligned.header)

tagged, untagged = 0, 0
for read in aligned:
    if read.query_name in meth_tags:
        for tag, val in meth_tags[read.query_name].items():
            read.set_tag(tag, val)
        tagged += 1
    else:
        untagged += 1
    out.write(read)

aligned.close()
out.close()
print(f"  Tagged: {tagged:,} | Untagged (no meth data): {untagged:,}")
if untagged > 0:
    pct = (untagged / (tagged + untagged)) * 100
    print(f"  Warning: {pct:.1f}% of reads had no MM/ML tags in uBAM")
PYEOF

log_info "  Coordinate-sorting final BAM..."
samtools sort -@ "$THREADS" -o "$FINAL_BAM" "$METH_UNSORTED"
samtools index "$FINAL_BAM"

STEP_END=$(date +%s)
log_success "Methylation tags transplanted in $((STEP_END - STEP_START))s → $FINAL_BAM"

# ── QC: Flagstat ──────────────────────────────────────────────────────────────
echo ""
log_info "Flagstat summary:"
samtools flagstat "$FINAL_BAM" | tee "${OUTDIR}/${SAMPLE}_flagstat.txt"

# ── QC: Coverage per chromosome ───────────────────────────────────────────────
echo ""
log_info "Coverage per canonical chromosome (top 10):"
samtools idxstats "$FINAL_BAM" | \
    grep -E "^chr[0-9XY]+\s" | \
    awk '{print $1, $3}' | sort -k2 -rn | head -10 | \
    awk '{printf "  %-10s %d reads\n", $1, $2}'

# ── QC: Verify MM tags ────────────────────────────────────────────────────────
echo ""
log_info "Verifying MM/ML methylation tags in final BAM..."
MM_CHECK=$(samtools view "$FINAL_BAM" | head -50 | grep -o "MM:Z:[^ ]*" | head -3 || true)
if [[ -n "$MM_CHECK" ]]; then
    log_success "MM tags confirmed present:"
    echo "$MM_CHECK" | sed 's/^/    /'
else
    log_warn "MM tags not detected — check if uBAM contained methylation data (requires Dorado with --modified-bases)"
fi

# ── Cleanup ───────────────────────────────────────────────────────────────────
if [[ "$CLEANUP" == true ]]; then
    echo ""
    log_info "Cleaning up intermediate files..."
    rm -f "$FASTQ" \
          "$SAM" \
          "$ALIGNED_SORTED" "${ALIGNED_SORTED}.bai" \
          "$UBAM_NAMESORTED" \
          "$ALIGNED_NAMESORTED" \
          "$METH_UNSORTED"
    [[ "$MERGED_CHUNKS" == true ]] && rm -f "${OUTDIR}/${SAMPLE}_chunk_list.txt"
    log_success "Intermediate files removed"
fi

# ── Final summary ─────────────────────────────────────────────────────────────
PIPELINE_END=$(date +%s)
TOTAL=$((PIPELINE_END - PIPELINE_START))
MINS=$((TOTAL / 60))
SECS=$((TOTAL % 60))

echo ""
echo -e "${GREEN}======================================================${NC}"
echo -e "${GREEN}   Pipeline Complete! v${VERSION}${NC}"
echo -e "${GREEN}======================================================${NC}"
echo "  Final BAM    : $FINAL_BAM"
echo "  Flagstat     : ${OUTDIR}/${SAMPLE}_flagstat.txt"
echo "  Total time   : ${MINS}m ${SECS}s"
echo -e "${GREEN}======================================================${NC}"
