#!/bin/bash

###############################################################################
# seq_qual_metrics.sh
#
# Per-window read-depth and mapping-quality metrics per pooled sample, on the same
# window grid as grenedalf interval windows (start = 1 + k*step; the last window per
# chromosome is the first whose end reaches the chromosome end; end is not clipped),
# so collate can join them onto diversity/FST windows by (chr, start, end, sample).
# Coordinates are 1-based inclusive.
#
# One pass per sample, split into genome chunks that run in parallel (--threads):
#   - samtools view   -> aligned bases (CIGAR M/=/X) per bin: all, MAPQ >= --min-mapq,
#                        MAPQ 0, MAPQ-weighted; primary reads, pairs, proper pairs,
#                        mates on another chromosome
#   - samtools depth  -> per-base depth at MAPQ >= --min-mapq: bases with depth >= --min-depth
# Bins are gcd(all window and step sizes) wide, so every --window-size/--step-size scale
# is built from the same pass. Reads flagged unmapped, secondary, QC-fail or duplicate are
# ignored (samtools depth defaults); supplementary alignments add bases but not read counts.
#
# Author: ssol-poolseq
# Usage: See below or run with --help
###############################################################################

set -euo pipefail

BAM_FILES=()
SAMPLE_INFO_CSV=""
OUTPUT_DIR=""
REFERENCE_GENOME=""
WINDOW_SIZES=()
STEP_SIZES=()
THREADS=1
MAX_JOBS=1
MIN_MAPQ=20
MIN_DEPTH=1
DRY_RUN=false
VERBOSE=false
SAMPLE_NAMES=()
USE_UNFILTERED=true

# Log file (set after output-dir is known)
SCRIPT_NAME=$(basename "$0")
LOG_FILE=""
LOG_REDIRECT_ACTIVE=false

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Input (choose one):
  -b, --bam FILE              BAM file(s) - can be specified multiple times
  -i, --sample-info FILE      CSV with columns: sample_name, bam_file (or bam_file path)

Required:
  -o, --output-dir DIR        Output directory for TSV files

Optional:
  --reference-genome FILE     Reference FASTA; its .fai sets chromosome order in output (default: BAM header order)
  -w, --window-size N         Window size in bp (can be specified multiple times; default: 10000)
  --step-size N               Step size in bp (can be specified multiple times; default: window/2)
  --min-mapq N                MAPQ threshold for mean_coverage_hq and breadth (default: $MIN_MAPQ)
  --min-depth N               Depth (at MAPQ >= --min-mapq) for a base to count toward breadth (default: $MIN_DEPTH)
  -t, --threads N             Genome chunks processed in parallel per sample (default: 1);
                              each chunk runs samtools view + samtools depth (~2 cores)
  -j, --max-jobs N            Run up to N samples in parallel (default: 1); total cores ~ 2 x threads x max-jobs
  --no-unfiltered             Use the listed BAMs as given. By default, X.unfiltered.bam is used in place
                              of X.bam when it exists (kept by process_poolseq.sh --bam-min-mapq /
                              --bam-primary-only), since a MAPQ-filtered BAM has no MAPQ-0 reads to measure.
  --dry-run                   Preview commands without executing
  --verbose                   Verbose diagnostics
  -h, --help                  Show this help

Output (per sample, per window/step combination):
  {output_dir}/{sample}/seq_qual_metrics_w{W}_s{S}.tsv
  Columns:
    chr, start, end, sample      window (1-based inclusive; grenedalf grid)
    mean_coverage                mean read depth over window bases in the chromosome (all MAPQ)
    mean_mapping_quality         mean MAPQ weighted by aligned bases
    mean_coverage_hq             mean depth from reads with MAPQ >= --min-mapq
    frac_mapq0                   share of aligned bases from MAPQ 0 reads (multi-mapping / repeats)
    breadth                      fraction of window bases with depth >= --min-depth at MAPQ >= --min-mapq
    rel_coverage                 mean_coverage / sample median mean_coverage (windows with coverage > 0):
                                 ~1 single copy, ~2 collapsed duplication, << 1 deletion / poor mappability
    n_reads                      primary mapped reads starting in the window
    frac_proper_pair             proper pairs / paired reads starting in the window
    frac_mate_other_chr          paired reads whose mate maps to another chromosome
  Logs: {output_dir}/log/seq_qual_metrics_YYYYmmdd_HHMMSS.log

Example:
  $0 --bam sample1.bam --bam sample2.bam -o out --window-size 10000 --step-size 5000 -t 6
  $0 --sample-info sample_info.csv -o out --window-size 10000 --reference-genome /path/to/ref.fa -t 4 -j 2
  $0 --bam s.bam -o out --window-size 1000 --window-size 10000 --step-size 500 --step-size 5000
EOF
}

log() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} ${LOG_PREFIX:-}$1"
    if [[ -n "$LOG_FILE" && -n "$1" && "$LOG_REDIRECT_ACTIVE" == false ]]; then
        echo "[$(date +'%Y-%m-%d %H:%M:%S')] ${LOG_PREFIX:-}$1" >> "$LOG_FILE"
    fi
}
log_error() {
    echo -e "${RED}[ERROR]${NC} ${LOG_PREFIX:-}$1" >&2
    if [[ -n "$LOG_FILE" && -n "$1" && "$LOG_REDIRECT_ACTIVE" == false ]]; then
        echo "[ERROR] ${LOG_PREFIX:-}$1" >> "$LOG_FILE"
    fi
}
log_warn() {
    echo -e "${YELLOW}[WARN]${NC} ${LOG_PREFIX:-}$1"
    if [[ -n "$LOG_FILE" && -n "$1" && "$LOG_REDIRECT_ACTIVE" == false ]]; then
        echo "[WARN] ${LOG_PREFIX:-}$1" >> "$LOG_FILE"
    fi
}
log_dry_run() {
    echo -e "${YELLOW}[DRY-RUN]${NC} $1"
    if [[ -n "$LOG_FILE" && -n "$1" && "$LOG_REDIRECT_ACTIVE" == false ]]; then
        echo "[DRY-RUN] $1" >> "$LOG_FILE"
    fi
}

check_command() {
    command -v "$1" &>/dev/null || { log_error "$1 not found"; exit 1; }
}

gcd() {
    local a=$1 b=$2 t
    while (( b > 0 )); do t=$((a % b)); a=$b; b=$t; done
    echo "$a"
}

# Parse args
while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--bam) BAM_FILES+=("$2"); shift 2 ;;
        -i|--sample-info) SAMPLE_INFO_CSV="$2"; shift 2 ;;
        -o|--output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --reference-genome) REFERENCE_GENOME="$2"; shift 2 ;;
        -w|--window-size) WINDOW_SIZES+=("$2"); shift 2 ;;
        --step-size) STEP_SIZES+=("$2"); shift 2 ;;
        --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
        --min-depth) MIN_DEPTH="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        -j|--max-jobs) MAX_JOBS="$2"; shift 2 ;;
        --no-unfiltered) USE_UNFILTERED=false; shift ;;
        --dry-run) DRY_RUN=true; shift ;;
        --verbose) VERBOSE=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "Unknown option: $1"; usage; exit 1 ;;
    esac
done

if [[ -z "$OUTPUT_DIR" ]]; then
    log_error "Missing --output-dir"
    usage
    exit 1
fi
for pair in "max-jobs:$MAX_JOBS" "threads:$THREADS"; do
    if ! [[ "${pair#*:}" =~ ^[1-9][0-9]*$ ]]; then
        log_error "Invalid --${pair%%:*} (must be a positive integer): ${pair#*:}"
        exit 1
    fi
done
if ! [[ "$MIN_MAPQ" =~ ^[0-9]+$ && "$MIN_DEPTH" =~ ^[1-9][0-9]*$ ]]; then
    log_error "--min-mapq must be >= 0 and --min-depth >= 1 (integers)"
    exit 1
fi

# Initialize log file
mkdir -p "$OUTPUT_DIR"
LOG_DIR="${OUTPUT_DIR}/log"
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +'%Y%m%d_%H%M%S')
LOG_FILE="${LOG_DIR}/${SCRIPT_NAME%.*}_${TIMESTAMP}.log"
if [[ "$DRY_RUN" == false ]]; then
    touch "$LOG_FILE"
    {
        echo "=========================================="
        echo "$SCRIPT_NAME"
        echo "Started: $(date +'%Y-%m-%d %H:%M:%S')"
        echo "Log file: $LOG_FILE"
        echo "=========================================="
    } >> "$LOG_FILE"
    log "Log file: $LOG_FILE"
    exec 3>&1 4>&2
    exec 1> >(tee -a "$LOG_FILE" >&3)
    exec 2>&1
    LOG_REDIRECT_ACTIVE=true
fi

# Resolve BAMs and sample names from CSV if provided
if [[ -n "$SAMPLE_INFO_CSV" ]]; then
    if [[ ! -f "$SAMPLE_INFO_CSV" ]]; then
        log_error "Sample info file not found: $SAMPLE_INFO_CSV"
        exit 1
    fi
    csv_dir=$(dirname "$SAMPLE_INFO_CSV")
    while IFS=',' read -r col1 col2 col3 col4 col5 rest; do
        # Skip header row (first column is sample_name, Sample, bam_file, etc.)
        case "$col1" in
            [Ss]ample_name|[Ss]ample|[Bb]am_file|"# Sample") continue ;;
        esac
        name="${col1}"
        bam="${col5}"
        [[ -z "$bam" ]] && bam="${col2}"
        [[ "$bam" != /* ]] && bam="${csv_dir}/${bam}"
        if [[ -f "$bam" ]]; then
            SAMPLE_NAMES+=("$name")
            BAM_FILES+=("$bam")
        else
            log_warn "BAM not found for sample $name: $bam (skipped)"
        fi
    done < <(grep -v '^#' "$SAMPLE_INFO_CSV" | tr -d '"\r')
fi

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    log_error "No BAM files specified. Use --bam or --sample-info."
    usage
    exit 1
fi

# process_poolseq.sh keeps the pre-filter BAM beside a MAPQ/primary-filtered one; the metrics
# (frac_mapq0, mean_coverage vs mean_coverage_hq) need the reads the filter removed.
if [[ "$USE_UNFILTERED" == true ]]; then
    for i in "${!BAM_FILES[@]}"; do
        unf="${BAM_FILES[$i]%.bam}.unfiltered.bam"
        if [[ "${BAM_FILES[$i]}" != *.unfiltered.bam && -f "$unf" ]]; then
            log "Using unfiltered BAM $unf in place of ${BAM_FILES[$i]}"
            BAM_FILES[$i]="$unf"
        fi
    done
fi

check_command samtools

# Build window/step arrays (default: one scale 10000/5000; else step = window/2 when not specified)
declare -a WINDOW_SIZE_ARRAY
declare -a STEP_SIZE_ARRAY
if [[ ${#WINDOW_SIZES[@]} -eq 0 ]]; then
    WINDOW_SIZE_ARRAY=(10000)
    STEP_SIZE_ARRAY=(5000)
    log "Using default window=10000, step=5000 (step = window/2)"
else
    WINDOW_SIZE_ARRAY=("${WINDOW_SIZES[@]}")
    STEP_SIZE_ARRAY=()
    for i in "${!WINDOW_SIZE_ARRAY[@]}"; do
        if [[ ${#STEP_SIZES[@]} -gt i && -n "${STEP_SIZES[$i]:-}" ]]; then
            STEP_SIZE_ARRAY[$i]="${STEP_SIZES[$i]}"
        else
            STEP_SIZE_ARRAY[$i]=$((${WINDOW_SIZE_ARRAY[$i]} / 2))
        fi
        log "Window/step scale $((i+1)): window=${WINDOW_SIZE_ARRAY[$i]}, step=${STEP_SIZE_ARRAY[$i]}"
    done
fi
BIN=0
for i in "${!WINDOW_SIZE_ARRAY[@]}"; do
    for v in "${WINDOW_SIZE_ARRAY[$i]}" "${STEP_SIZE_ARRAY[$i]}"; do
        if ! [[ "$v" =~ ^[1-9][0-9]*$ ]]; then
            log_error "Window and step sizes must be positive integers (got: $v)"
            exit 1
        fi
        BIN=$(gcd "$v" "$BIN")
    done
done
WINDOWS_CSV=$(IFS=,; echo "${WINDOW_SIZE_ARRAY[*]}")
STEPS_CSV=$(IFS=,; echo "${STEP_SIZE_ARRAY[*]}")
log "Bin size: ${BIN} bp (gcd of window and step sizes); MAPQ threshold: ${MIN_MAPQ}; breadth depth: ${MIN_DEPTH}"
log "Parallelism: ${THREADS} chunk(s) per sample, ${MAX_JOBS} sample(s) at a time"
if (( BIN < 50 )); then
    log_warn "Bin size ${BIN} bp is small (window/step sizes share few factors); binning will use more memory and time"
fi

# --- awk programs ----------------------------------------------------------------------

# samtools view lines -> per-bin sums. Regions file (idx, chr, start, end) is read first; aligned
# bases are clipped to the chunk's region so reads spanning chunk edges are counted once.
read -r -d '' AWK_VIEW << 'EOF' || true
BEGIN { OFS = "\t" }
FNR == NR { RI[$2] = $1; RBEG[$2] = $3; REND[$2] = $4; next }
function addseg(p, n,    b, room, k) {
    if (p < rs) { n -= rs - p; p = rs }
    if (p + n - 1 > re) n = re - p + 1
    while (n > 0) {
        b = int((p - 1) / g); room = (b + 1) * g - p + 1; k = (n < room) ? n : room
        key = ci SUBSEP b; seen[key] = 1
        A[key] += k; MW[key] += k * mq
        if (mq >= q) H[key] += k
        if (mq == 0) Z[key] += k
        p += k; n -= k
    }
}
{
    chr = $3; if (!(chr in RI)) next
    ci = RI[chr]; rs = RBEG[chr]; re = REND[chr]
    flag = $2 + 0; pos = $4 + 0; mq = $5 + 0; cig = $6
    if (cig == "*") next
    p = pos
    while (match(cig, /^[0-9]+[MIDNSHP=X]/)) {
        n = substr(cig, 1, RLENGTH - 1) + 0; op = substr(cig, RLENGTH, 1); cig = substr(cig, RLENGTH + 1)
        if (op == "M" || op == "=" || op == "X") { addseg(p, n); p += n }
        else if (op == "D" || op == "N") p += n
    }
    if (pos >= rs && pos <= re && int(flag / 2048) % 2 == 0) {
        key = ci SUBSEP int((pos - 1) / g); seen[key] = 1
        R[key]++
        if (flag % 2 == 1) {
            P[key]++
            if (int(flag / 2) % 2 == 1) PP[key]++
            if (int(flag / 8) % 2 == 0 && $7 != "=" && $7 != chr) MO[key]++
        }
    }
}
END {
    for (key in seen) {
        split(key, kk, SUBSEP)
        print kk[1], kk[2], "V", A[key] + 0, H[key] + 0, Z[key] + 0, MW[key] + 0, R[key] + 0, P[key] + 0, PP[key] + 0, MO[key] + 0
    }
}
EOF

# samtools depth -a lines -> per-bin counts of bases with depth >= md
read -r -d '' AWK_DEPTH << 'EOF' || true
BEGIN { OFS = "\t" }
FNR == NR { RI[$2] = $1; next }
{
    key = RI[$1] SUBSEP int(($2 - 1) / g)
    if ($3 >= md) C[key]++
    else if (!(key in C)) C[key] = 0
}
END { for (key in C) { split(key, kk, SUBSEP); print kk[1], kk[2], "D", C[key] } }
EOF

# Sorted bin lines (idx, bin, type, values) -> one window table per scale, per chromosome via
# prefix sums. Chromosome table (idx, chr, length) is read first.
read -r -d '' AWK_WINDOWS << 'EOF' || true
BEGIN {
    OFS = "\t"
    ns = split(wins, W, ","); split(steps, S, ",")
    for (j = 1; j <= ns; j++) out[j] = prefix "_" j ".tsv"
    cur = ""
}
FNR == NR { NAME[$1] = $2; LEN[$1] = $3; next }
function fmt(x) { return (x == "NA") ? "NA" : sprintf("%.4g", x) }
function flush(ci,    L, nb, b, f, j, w, s, start, end, e, b0, b1, nbases, a, h, z, mw, r, pp, pr, mo, c) {
    if (ci == "") return
    L = LEN[ci]; nb = int((L + g - 1) / g)
    for (f = 1; f <= 9; f++) cum[f, -1] = 0
    for (b = 0; b < nb; b++)
        for (f = 1; f <= 9; f++) cum[f, b] = cum[f, b - 1] + V[f, b]
    for (j = 1; j <= ns; j++) {
        w = W[j]; s = S[j]
        for (start = 1; ; start += s) {
            end = start + w - 1; e = (end < L) ? end : L
            b0 = int((start - 1) / g); b1 = int((e + g - 1) / g) - 1
            nbases = e - start + 1
            a = cum[1, b1] - cum[1, b0 - 1]; h = cum[2, b1] - cum[2, b0 - 1]; z = cum[3, b1] - cum[3, b0 - 1]
            mw = cum[4, b1] - cum[4, b0 - 1]; r = cum[5, b1] - cum[5, b0 - 1]; pr = cum[6, b1] - cum[6, b0 - 1]
            pp = cum[7, b1] - cum[7, b0 - 1]; mo = cum[8, b1] - cum[8, b0 - 1]; c = cum[9, b1] - cum[9, b0 - 1]
            print ci, NAME[ci], start, end, fmt(a / nbases), fmt(a > 0 ? mw / a : "NA"), fmt(h / nbases), \
                  fmt(a > 0 ? z / a : "NA"), fmt(c / nbases), r, fmt(pr > 0 ? pp / pr : "NA"), fmt(pr > 0 ? mo / pr : "NA") > out[j]
            if (end >= L) break
        }
    }
    delete V; delete cum
}
{
    if ($1 != cur) { flush(cur); cur = $1 }
    if ($3 == "V") { for (f = 1; f <= 8; f++) V[f, $2] += $(f + 3) }
    else V[9, $2] += $4
}
END { flush(cur) }
EOF

# --- per-sample processing ---------------------------------------------------------------

# Process one sample: index BAM if needed, compute bins over parallel chunks, build all scales.
process_one_sample() {
    local bam="$1" sample_name="$2"
    local sample_out tmp chroms total target n_chunks i j w s out_tsv median
    sample_out="${OUTPUT_DIR}/${sample_name}"

    if [[ "$DRY_RUN" == true ]]; then
        log_dry_run "Sample $sample_name: $bam"
        log_dry_run "  per chunk (${THREADS} in parallel): samtools view -F 0x704 <bam> <regions> | awk (bins of ${BIN} bp)"
        log_dry_run "  per chunk: samtools depth -a -Q ${MIN_MAPQ} -r <region> <bam> | awk (bases with depth >= ${MIN_DEPTH})"
        for i in "${!WINDOW_SIZE_ARRAY[@]}"; do
            log_dry_run "  -> ${sample_out}/seq_qual_metrics_w${WINDOW_SIZE_ARRAY[$i]}_s${STEP_SIZE_ARRAY[$i]}.tsv"
        done
        return 0
    fi

    mkdir -p "$sample_out"
    tmp=$(mktemp -d "${sample_out}/.tmp_seq_qual.XXXXXX")
    trap 'kill $(jobs -rp) 2>/dev/null || true; rm -rf "'"$tmp"'"' EXIT

    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" && ! -f "${bam}.csi" ]]; then
        log "Indexing BAM: $bam"
        samtools index -@ "$THREADS" "$bam" || { log_error "samtools index failed for $bam"; return 1; }
    fi

    # Chromosome table: idx, name, length (reference .fai order when given, else BAM header order)
    samtools view -H "$bam" | awk -F'\t' '$1 == "@SQ" { for (i = 2; i <= NF; i++) { if ($i ~ /^SN:/) n = substr($i, 4); if ($i ~ /^LN:/) l = substr($i, 4) } print n "\t" l }' > "${tmp}/header.tsv"
    if [[ -n "$REFERENCE_GENOME" && -f "${REFERENCE_GENOME}.fai" ]]; then
        # Window ends come from these lengths; they must match the reference grenedalf uses
        local mismatch
        mismatch=$(awk -F'\t' 'NR == FNR { len[$1] = $2; next } ($1 in len) && len[$1] != $2 { n++; if (n <= 3) ex = ex " " $1 ":" $2 "!=" len[$1] } END { if (n) print n ex }' \
            "${REFERENCE_GENOME}.fai" "${tmp}/header.tsv")
        if [[ -n "$mismatch" ]]; then
            log_error "$bam: ${mismatch%% *} chromosome length(s) differ from ${REFERENCE_GENOME}.fai (BAM!=fai:${mismatch#* }); was it aligned to another assembly?"
            return 1
        fi
        awk -F'\t' -v OFS='\t' 'NR == FNR { ord[$1] = NR; next } { print ($1 in ord) ? ord[$1] : 1000000 + FNR, $1, $2 }' \
            "${REFERENCE_GENOME}.fai" "${tmp}/header.tsv" | sort -k1,1n > "${tmp}/chroms.tsv"
    else
        awk -F'\t' -v OFS='\t' '{ print NR, $1, $2 }' "${tmp}/header.tsv" > "${tmp}/chroms.tsv"
    fi
    chroms=$(wc -l < "${tmp}/chroms.tsv")
    [[ "$chroms" -gt 0 ]] || { log_error "No @SQ lines in $bam"; return 1; }
    total=$(awk '{ s += $3 } END { print s }' "${tmp}/chroms.tsv")

    # Chunks: chromosomes split at bin boundaries into pieces of ~total/(4*threads) (>= 1 Mb);
    # small chromosomes are packed together. Each chunk lists each chromosome at most once.
    target=$(( total / (4 * THREADS) ))
    (( target < 1000000 )) && target=1000000
    target=$(( (target + BIN - 1) / BIN * BIN ))
    awk -F'\t' -v OFS='\t' -v t="$target" -v dir="$tmp" '
        { idx = $1; chr = $2; L = $3
          for (s = 1; s <= L; s += t) {
              e = s + t - 1; if (e > L) e = L
              if (acc + (e - s + 1) > t && acc > 0) { n++; acc = 0 }
              print idx, chr, s, e > (dir "/chunk_" n ".regions")
              acc += e - s + 1
              if (acc >= t) { n++; acc = 0 }
          } }
        END { }' "${tmp}/chroms.tsv"
    n_chunks=$(ls "${tmp}"/chunk_*.regions | wc -l)
    log "Sample $sample_name: ${chroms} chromosome(s), ${total} bp, ${n_chunks} chunk(s)"

    run_chunk() {
        local reg="$1" base="${1%.regions}"
        while IFS=$'\t' read -r _ c s e; do
            samtools view -F 0x704 "$bam" "{${c}}:${s}-${e}"
        done < "$reg" | awk -v g="$BIN" -v q="$MIN_MAPQ" "$AWK_VIEW" "$reg" - > "${base}.bins" || return 1
        while IFS=$'\t' read -r _ c s e; do
            samtools depth -a -Q "$MIN_MAPQ" -r "{${c}}:${s}-${e}" "$bam"
        done < "$reg" | awk -v g="$BIN" -v md="$MIN_DEPTH" "$AWK_DEPTH" "$reg" - >> "${base}.bins" || return 1
    }

    local reg
    for reg in "${tmp}"/chunk_*.regions; do
        while (( $(jobs -rp | wc -l) >= THREADS )); do
            wait -n || true
        done
        { run_chunk "$reg" || touch "${reg%.regions}.failed"; } &
    done
    wait
    if compgen -G "${tmp}/chunk_*.failed" > /dev/null; then
        log_error "$(ls "${tmp}"/chunk_*.failed | wc -l) chunk(s) failed for $sample_name ($bam)"
        return 1
    fi

    cat "${tmp}"/chunk_*.bins | sort -t$'\t' -k1,1n -k2,2n -S 1G > "${tmp}/bins.sorted"
    awk -F'\t' -v g="$BIN" -v wins="$WINDOWS_CSV" -v steps="$STEPS_CSV" -v prefix="${tmp}/scale" \
        "$AWK_WINDOWS" "${tmp}/chroms.tsv" "${tmp}/bins.sorted"

    for i in "${!WINDOW_SIZE_ARRAY[@]}"; do
        j=$((i + 1)); w="${WINDOW_SIZE_ARRAY[$i]}"; s="${STEP_SIZE_ARRAY[$i]}"
        out_tsv="${sample_out}/seq_qual_metrics_w${w}_s${s}.tsv"
        median=$(awk -F'\t' '$5 > 0 { print $5 }' "${tmp}/scale_${j}.tsv" | sort -g | awk '{ v[NR] = $1 } END { if (NR == 0) print "NA"; else if (NR % 2) print v[(NR + 1) / 2]; else print (v[NR / 2] + v[NR / 2 + 1]) / 2 }')
        {
            printf 'chr\tstart\tend\tsample\tmean_coverage\tmean_mapping_quality\tmean_coverage_hq\tfrac_mapq0\tbreadth\trel_coverage\tn_reads\tfrac_proper_pair\tfrac_mate_other_chr\n'
            awk -F'\t' -v OFS='\t' -v sn="$sample_name" -v med="$median" '{
                rel = (med == "NA" || med == 0) ? "NA" : sprintf("%.4g", $5 / med)
                print $2, $3, $4, sn, $5, $6, $7, $8, $9, rel, $10, $11, $12 }' "${tmp}/scale_${j}.tsv"
        } > "$out_tsv"
        log "  Wrote $out_tsv ($(( $(wc -l < "$out_tsv") - 1 )) windows; median coverage ${median})"
        if [[ "$VERBOSE" == true ]]; then
            head -n 3 "$out_tsv" | while IFS= read -r line; do log "    $line"; done
        fi
    done
}

# Start processing
if [[ "$DRY_RUN" == true ]]; then
    log_dry_run "DRY-RUN MODE: Commands will be previewed but not executed"
    log_dry_run "Samples to process: ${#BAM_FILES[@]}; window/step combinations: ${#WINDOW_SIZE_ARRAY[@]}"
fi

STATUS_DIR="${LOG_DIR}/.status_${TIMESTAMP}"
mkdir -p "$STATUS_DIR"
declare -a RUN_SAMPLES=()
for i in "${!BAM_FILES[@]}"; do
    bam="${BAM_FILES[$i]}"
    if [[ ${#SAMPLE_NAMES[@]} -gt i && -n "${SAMPLE_NAMES[$i]:-}" ]]; then
        sample_name="${SAMPLE_NAMES[$i]}"
    else
        sample_name=$(basename "$bam" .bam)
        sample_name="${sample_name%.unfiltered}"
        sample_name="${sample_name%_All_seq.dedup}"
        sample_name="${sample_name%.dedup}"
    fi

    if [[ "$DRY_RUN" == true ]]; then
        process_one_sample "$bam" "$sample_name"
        continue
    fi
    while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
        wait -n || true
    done
    # errexit stays active in the subshell; the status file is written only when the sample completes
    ( export LOG_PREFIX="[$sample_name] "; process_one_sample "$bam" "$sample_name"; touch "${STATUS_DIR}/${i}.ok" ) &
    RUN_SAMPLES[$i]="$sample_name"
done
wait
FAILED_SAMPLES=()
for i in "${!RUN_SAMPLES[@]}"; do
    [[ -f "${STATUS_DIR}/${i}.ok" ]] || FAILED_SAMPLES+=("${RUN_SAMPLES[$i]}")
done
rm -rf "$STATUS_DIR"

if [[ "$DRY_RUN" == true ]]; then
    log_dry_run "DRY-RUN complete. Run without --dry-run to execute the analysis."
elif [[ ${#FAILED_SAMPLES[@]} -gt 0 ]]; then
    log_error "Failed samples: ${FAILED_SAMPLES[*]} (see $LOG_FILE)"
    exit 1
else
    log "Done. Output under $OUTPUT_DIR"
fi
