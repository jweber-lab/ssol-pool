#!/bin/bash

###############################################################################
# build_mask.sh
#
# Build one genome mask BED that combines repeat annotation with windows flagged by
# alignment-quality metrics (seq_qual_metrics.sh), for --filter-mask-total-bed in
# calculate_pi_theta.sh, calculate_fst.sh, snp_call.sh and indel_call.sh.
#
# Repeat sources (any number): BED files, RepeatMasker .out files, soft-masked FASTA
# (lowercase runs). Seq-qual source: per-sample seq_qual_metrics_w{W}_s{S}.tsv at one
# scale (overlapping windows, step <= w/4, recommended); a window is flagged when any sample
# (or at least --min-sample-fraction of samples) fails a threshold (high MAPQ-0 share, collapsed
# duplication or depleted coverage relative to the sample median, low high-MAPQ breadth), and a
# base is masked when at least --min-window-vote of the windows covering it are flagged.
#
# Outputs (in --output-dir):
#   mask.bed                    combined mask: sorted (reference order), merged, 0-based half-open
#   mask_sources.bed            per-source merged intervals; 4th column repeat or seq_qual
#   seq_qual_flagged_windows.tsv  chr, start, end, n_flagged, n_samples, reasons (when --seq-qual-dir)
#   mask_summary.tsv            per chromosome and genome: length, repeat_bp, seq_qual_bp, masked_bp, frac_masked
###############################################################################

set -euo pipefail

REFERENCE=""
OUTPUT_DIR=""
REPEAT_BEDS=()
REPEATMASKER_OUTS=()
SOFTMASK_FASTAS=()
SEQ_QUAL_DIR=""
SEQ_QUAL_SCALE=""
MAX_FRAC_MAPQ0=0.25
MAX_REL_COVERAGE=2
MIN_REL_COVERAGE=0.2
MIN_BREADTH=0.5
MIN_MEAN_MAPQ=0
MIN_SAMPLE_FRACTION=0
MIN_WINDOW_VOTE=0.5
MERGE_DISTANCE=0
MIN_LENGTH=1
DRY_RUN=false

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

log() { echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }

usage() {
    cat << EOF
Usage: $0 --reference REF -o DIR [sources] [options]

Required:
  -r, --reference FILE        Reference FASTA (with .fai) or the .fai itself; sets chromosome
                              names, lengths and output order
  -o, --output-dir DIR        Output directory

Repeat sources (repeatable; any combination):
  --repeat-bed FILE           BED (0-based half-open; extra columns ignored), e.g. from repeatmasker_bed.sh
  --repeatmasker-out FILE     RepeatMasker .out (converted with repeatmasker_bed.sh)
  --softmask-fasta FILE       Soft-masked FASTA; lowercase runs are masked

Alignment-quality source (seq_qual_metrics.sh output):
  --seq-qual-dir DIR          Directory searched recursively for seq_qual_metrics_w*_s*.tsv
  --seq-qual-scale wW_sS      Scale to use (default: the smallest step, then smallest window; overlapping
                              windows with step <= w/4, e.g. -w 1000 --step-size 250, are recommended)
  --max-frac-mapq0 X          Flag when frac_mapq0 > X (default: $MAX_FRAC_MAPQ0)
  --max-rel-coverage X        Flag when rel_coverage > X, collapsed duplications (default: $MAX_REL_COVERAGE)
  --min-rel-coverage X        Flag when rel_coverage < X, deletions / unmappable (default: $MIN_REL_COVERAGE)
  --min-breadth X             Flag when breadth < X (default: $MIN_BREADTH)
  --min-mean-mapq X           Flag when mean_mapping_quality < X (default: $MIN_MEAN_MAPQ = off)
  --min-sample-fraction X     Flag a window when >= X of samples fail a threshold, and at least one
                              (default: $MIN_SAMPLE_FRACTION = any failing sample flags the window for all).
                              Raise it (e.g. 0.5) with many pools so one pool's dropout is not masked everywhere.
  --min-window-vote X         Mask a base when >= X of the windows covering it are flagged (default: $MIN_WINDOW_VOTE).
                              Metrics are window averages, so overlapping windows catch problems that
                              straddle a boundary; the vote keeps the mask near the problem instead of
                              masking every flagged window in full. Use a small value (e.g. 0.01) for
                              that union instead.
                              Any threshold can be disabled with a value outside the metric range
                              (e.g. --max-frac-mapq0 1, --min-breadth 0).

Options:
  --merge-distance N          Merge intervals closer than N bp (default: $MERGE_DISTANCE = touching/overlapping only)
  --min-length N              Drop combined intervals shorter than N bp after merging (default: $MIN_LENGTH)
  --dry-run                   Show inputs and settings without writing
  -h, --help                  Show this help

Example:
  $0 --reference ref.unmasked.fa --repeatmasker-out ref.fa.out \\
     --seq-qual-dir seq_qual_mask --output-dir mask
  # seq_qual_mask from: seq_qual_metrics.sh ... -w 1000 --step-size 250 -o seq_qual_mask
  # then: calculate_pi_theta.sh / calculate_fst.sh / snp_call.sh / indel_call.sh --filter-mask-total-bed mask/mask.bed
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--reference) REFERENCE="$2"; shift 2 ;;
        -o|--output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --repeat-bed) REPEAT_BEDS+=("$2"); shift 2 ;;
        --repeatmasker-out) REPEATMASKER_OUTS+=("$2"); shift 2 ;;
        --softmask-fasta) SOFTMASK_FASTAS+=("$2"); shift 2 ;;
        --seq-qual-dir) SEQ_QUAL_DIR="$2"; shift 2 ;;
        --seq-qual-scale) SEQ_QUAL_SCALE="$2"; shift 2 ;;
        --max-frac-mapq0) MAX_FRAC_MAPQ0="$2"; shift 2 ;;
        --max-rel-coverage) MAX_REL_COVERAGE="$2"; shift 2 ;;
        --min-rel-coverage) MIN_REL_COVERAGE="$2"; shift 2 ;;
        --min-breadth) MIN_BREADTH="$2"; shift 2 ;;
        --min-mean-mapq) MIN_MEAN_MAPQ="$2"; shift 2 ;;
        --min-sample-fraction) MIN_SAMPLE_FRACTION="$2"; shift 2 ;;
        --min-window-vote) MIN_WINDOW_VOTE="$2"; shift 2 ;;
        --merge-distance) MERGE_DISTANCE="$2"; shift 2 ;;
        --min-length) MIN_LENGTH="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "Unknown option: $1"; usage; exit 1 ;;
    esac
done

[[ -n "$REFERENCE" && -n "$OUTPUT_DIR" ]] || { log_error "--reference and --output-dir are required"; usage; exit 1; }
if [[ ${#REPEAT_BEDS[@]} -eq 0 && ${#REPEATMASKER_OUTS[@]} -eq 0 && ${#SOFTMASK_FASTAS[@]} -eq 0 && -z "$SEQ_QUAL_DIR" ]]; then
    log_error "Give at least one source: --repeat-bed, --repeatmasker-out, --softmask-fasta or --seq-qual-dir"
    exit 1
fi

if [[ "$REFERENCE" == *.fai ]]; then FAI="$REFERENCE"; else FAI="${REFERENCE}.fai"; fi
[[ -f "$FAI" ]] || { log_error "Reference index not found: $FAI (run: samtools faidx $REFERENCE)"; exit 1; }
for f in "${REPEAT_BEDS[@]}" "${REPEATMASKER_OUTS[@]}" "${SOFTMASK_FASTAS[@]}"; do
    [[ -f "$f" ]] || { log_error "Input not found: $f"; exit 1; }
done
for v in "$MAX_FRAC_MAPQ0" "$MAX_REL_COVERAGE" "$MIN_REL_COVERAGE" "$MIN_BREADTH" "$MIN_MEAN_MAPQ" "$MIN_SAMPLE_FRACTION" "$MIN_WINDOW_VOTE"; do
    [[ "$v" =~ ^[0-9]*\.?[0-9]+$ ]] || { log_error "Thresholds must be non-negative numbers (got: $v)"; exit 1; }
done
[[ "$MERGE_DISTANCE" =~ ^[0-9]+$ && "$MIN_LENGTH" =~ ^[0-9]+$ ]] || { log_error "--merge-distance and --min-length must be integers"; exit 1; }

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

# Seq-qual inputs at one scale
SQ_FILES=()
if [[ -n "$SEQ_QUAL_DIR" ]]; then
    [[ -d "$SEQ_QUAL_DIR" ]] || { log_error "Seq-qual directory not found: $SEQ_QUAL_DIR"; exit 1; }
    mapfile -t scales < <(find "$SEQ_QUAL_DIR" -name 'seq_qual_metrics_w*_s*.tsv' -printf '%f\n' \
        | sed -n 's/^seq_qual_metrics_\(w[0-9]*_s[0-9]*\)\.tsv$/\1/p' | sort -u)
    [[ ${#scales[@]} -gt 0 ]] || { log_error "No seq_qual_metrics_w*_s*.tsv under $SEQ_QUAL_DIR"; exit 1; }
    if [[ -z "$SEQ_QUAL_SCALE" ]]; then
        SEQ_QUAL_SCALE=$(printf '%s\n' "${scales[@]}" | awk -F'[ws_]+' '{ print $3, $2, $0 }' | sort -k1,1n -k2,2n | head -1 | cut -d' ' -f3)
    elif ! printf '%s\n' "${scales[@]}" | grep -qx "$SEQ_QUAL_SCALE"; then
        log_error "Scale $SEQ_QUAL_SCALE not found in $SEQ_QUAL_DIR (found: ${scales[*]})"
        exit 1
    fi
    sq_w=$(sed 's/^w\([0-9]*\)_s.*/\1/' <<< "$SEQ_QUAL_SCALE"); sq_s=$(sed 's/^w[0-9]*_s//' <<< "$SEQ_QUAL_SCALE")
    if (( sq_s >= sq_w )); then
        log_warn "Seq-qual scale $SEQ_QUAL_SCALE has no overlap: each base is judged by one window, so problems straddling a window boundary can be missed (e.g. -w 1000 --step-size 250 is safer)"
    elif (( sq_w < 3 * sq_s )); then
        log_warn "Seq-qual scale $SEQ_QUAL_SCALE: each base is in only $(( (sq_w + sq_s - 1) / sq_s )) windows, so with --min-window-vote ${MIN_WINDOW_VOTE} any flagged window may mask its full ${sq_w} bp (a step of w/4 or smaller gives a tighter mask)"
    fi
    mapfile -t SQ_FILES < <(find "$SEQ_QUAL_DIR" -name "seq_qual_metrics_${SEQ_QUAL_SCALE}.tsv" | sort)
fi

log "Reference index: $FAI ($(wc -l < "$FAI") sequences)"
log "Repeat sources: ${#REPEAT_BEDS[@]} BED, ${#REPEATMASKER_OUTS[@]} RepeatMasker .out, ${#SOFTMASK_FASTAS[@]} soft-masked FASTA"
if [[ ${#SQ_FILES[@]} -gt 0 ]]; then
    log "Seq-qual source: ${#SQ_FILES[@]} sample file(s) at ${SEQ_QUAL_SCALE}; flag if frac_mapq0 > ${MAX_FRAC_MAPQ0}, rel_coverage > ${MAX_REL_COVERAGE} or < ${MIN_REL_COVERAGE}, breadth < ${MIN_BREADTH}, mean MAPQ < ${MIN_MEAN_MAPQ}; flag a window if >= ${MIN_SAMPLE_FRACTION} of samples (at least one); mask a base if >= ${MIN_WINDOW_VOTE} of the windows covering it are flagged"
fi
if [[ "$DRY_RUN" == true ]]; then
    log "DRY-RUN: would write ${OUTPUT_DIR}/{mask.bed,mask_sources.bed,mask_summary.tsv}"
    exit 0
fi

mkdir -p "$OUTPUT_DIR"
TMP=$(mktemp -d "${OUTPUT_DIR}/.build_mask.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

# Sort BED by reference order, clip to chromosome bounds, merge within distance d.
# Reads chromosome table (fai) first; intervals on unknown chromosomes are counted and dropped.
merge_bed() {
    local d="$1"; shift
    awk -F'\t' -v OFS='\t' 'NR == FNR { ord[$1] = NR; len[$1] = $2; next }
        /^(#|track|browser)/ || NF < 3 { next }
        !($1 in ord) { unk[$1]++; next }
        { s = $2 + 0; e = $3 + 0; if (s < 0) s = 0; if (e > len[$1]) e = len[$1]; if (e > s) print ord[$1], $1, s, e }
        END { n = 0; for (c in unk) { n++; if (n <= 3) ex = ex (n > 1 ? "," : "") c }
              if (n) printf "%d sequence name(s) not in the reference index were dropped (e.g. %s)\n", n, ex > "/dev/stderr" }' "$FAI" "$@" \
    | sort -t$'\t' -k1,1n -k3,3n -k4,4n \
    | awk -F'\t' -v OFS='\t' -v d="$d" '
        { if ($2 == c && $3 <= e + d) { if ($4 > e) e = $4; next }
          if (c != "") print c, s, e
          c = $2; s = $3; e = $4 }
        END { if (c != "") print c, s, e }'
}

# --- repeat sources ----------------------------------------------------------------
: > "${TMP}/repeat_raw.bed"
for f in "${REPEAT_BEDS[@]}"; do
    cut -f1-3 "$f" >> "${TMP}/repeat_raw.bed"
done
for f in "${REPEATMASKER_OUTS[@]}"; do
    bash "${SCRIPT_DIR}/repeatmasker_bed.sh" -i "$f" -o "${TMP}/rm.bed" --no-merge > /dev/null
    cat "${TMP}/rm.bed" >> "${TMP}/repeat_raw.bed"
done
for f in "${SOFTMASK_FASTAS[@]}"; do
    # Lowercase runs per line; runs split across lines are joined by the merge step
    { if [[ "$f" == *.gz ]]; then gzip -dc "$f"; else cat "$f"; fi; } | awk -v OFS='\t' '
        /^>/ { chr = substr($1, 2); pos = 0; next }
        { line = $0; off = 0
          while (match(line, /[acgtnrykmswbdhv]+/)) {
              print chr, pos + off + RSTART - 1, pos + off + RSTART - 1 + RLENGTH
              off += RSTART + RLENGTH - 1; line = substr(line, RSTART + RLENGTH)
          }
          pos += length($0) }' >> "${TMP}/repeat_raw.bed"
done
merge_bed 0 "${TMP}/repeat_raw.bed" > "${TMP}/repeat.bed"

# --- seq-qual source ---------------------------------------------------------------
: > "${TMP}/seq_qual.bed"
if [[ ${#SQ_FILES[@]} -gt 0 ]]; then
    awk -F'\t' -v OFS='\t' -v mq0="$MAX_FRAC_MAPQ0" -v rhi="$MAX_REL_COVERAGE" -v rlo="$MIN_REL_COVERAGE" \
        -v br="$MIN_BREADTH" -v mmq="$MIN_MEAN_MAPQ" -v frac="$MIN_SAMPLE_FRACTION" -v out="${TMP}/flagged.tsv" '
        FNR == 1 {
            delete col
            for (i = 1; i <= NF; i++) col[$i] = i
            need = "chr start end sample frac_mapq0 rel_coverage breadth mean_mapping_quality"
            n = split(need, nm, " ")
            for (i = 1; i <= n; i++) if (!(nm[i] in col)) {
                printf "%s lacks column %s; re-run seq_qual_metrics.sh\n", FILENAME, nm[i] > "/dev/stderr"; bad = 1; exit 1 }
            next
        }
        {
            key = $col["chr"] SUBSEP $col["start"] SUBSEP $col["end"]
            if (!(key in order)) { order[key] = ++nwin; keys[nwin] = key }
            samples[$col["sample"]] = 1
            r = ""
            v = $col["frac_mapq0"];           if (v != "NA" && v + 0 > mq0) r = r ",mapq0"
            v = $col["rel_coverage"];         if (v != "NA" && v + 0 > rhi) r = r ",high_coverage"
            if (v != "NA" && v + 0 < rlo) r = r ",low_coverage"
            v = $col["breadth"];              if (v != "NA" && v + 0 < br) r = r ",low_breadth"
            v = $col["mean_mapping_quality"]; if (v != "NA" && v + 0 < mmq) r = r ",low_mapq"
            if (r != "") {
                nflag[key]++
                n = split(substr(r, 2), rr, ",")
                for (i = 1; i <= n; i++) if (!((key, rr[i]) in seen)) { seen[key, rr[i]] = 1; why[key] = why[key] (why[key] == "" ? "" : ",") rr[i] }
            }
        }
        END {
            if (bad) exit 1
            ns = 0; for (s in samples) ns++
            need_n = frac * ns; if (need_n < 1) need_n = 1
            print "chr", "start", "end", "n_flagged", "n_samples", "reasons" > out
            for (k = 1; k <= nwin; k++) {
                key = keys[k]
                split(key, p, SUBSEP)
                f = (nflag[key] + 0 >= need_n - 1e-9)
                if (f) print p[1], p[2], p[3], nflag[key], ns, why[key] > out
                print p[1], p[2] - 1, 1, f
                print p[1], p[3], -1, -f
            }
        }' "${SQ_FILES[@]}" > "${TMP}/seq_qual_events.tsv"
    # Sweep window starts/ends: a segment is masked when >= vote of the windows covering it are flagged
    sort -t$'\t' -k1,1 -k2,2n "${TMP}/seq_qual_events.tsv" \
    | awk -F'\t' -v OFS='\t' -v vote="$MIN_WINDOW_VOTE" '
        function emit(to) { if (cov > 0 && to > pos && flg >= vote * cov - 1e-9 && flg > 0) print chr, pos, to }
        $1 != chr { chr = $1; pos = $2; cov = 0; flg = 0 }
        $2 != pos { emit($2); pos = $2 }
        { cov += $3; flg += $4 }' > "${TMP}/seq_qual_raw.bed"
    merge_bed 0 "${TMP}/seq_qual_raw.bed" > "${TMP}/seq_qual.bed"
    cp "${TMP}/flagged.tsv" "${OUTPUT_DIR}/seq_qual_flagged_windows.tsv"
fi

# --- combine -----------------------------------------------------------------------
cat "${TMP}/repeat.bed" "${TMP}/seq_qual.bed" > "${TMP}/all.bed"
merge_bed "$MERGE_DISTANCE" "${TMP}/all.bed" | awk -F'\t' -v m="$MIN_LENGTH" '$3 - $2 >= m' > "${OUTPUT_DIR}/mask.bed"
{ awk -v OFS='\t' '{ print $0, "repeat" }' "${TMP}/repeat.bed"
  awk -v OFS='\t' '{ print $0, "seq_qual" }' "${TMP}/seq_qual.bed"; } \
    | awk -F'\t' -v OFS='\t' 'NR == FNR { ord[$1] = NR; next } { print ord[$1], $0 }' "$FAI" - \
    | sort -t$'\t' -k1,1n -k3,3n -k5,5 | cut -f2- > "${OUTPUT_DIR}/mask_sources.bed"

awk -F'\t' -v OFS='\t' '
    FILENAME == ARGV[1] { name[++n] = $1; len[$1] = $2; next }
    FILENAME == ARGV[2] { rep[$1] += $3 - $2; next }
    FILENAME == ARGV[3] { sq[$1] += $3 - $2; next }
    { all[$1] += $3 - $2 }
    END {
        print "chr", "length", "repeat_bp", "seq_qual_bp", "masked_bp", "frac_masked"
        for (i = 1; i <= n; i++) {
            c = name[i]; L += len[c]; R += rep[c]; S += sq[c]; A += all[c]
            print c, len[c], rep[c] + 0, sq[c] + 0, all[c] + 0, sprintf("%.4f", len[c] ? all[c] / len[c] : 0)
        }
        print "genome", L, R + 0, S + 0, A + 0, sprintf("%.4f", L ? A / L : 0)
    }' "$FAI" "${TMP}/repeat.bed" "${TMP}/seq_qual.bed" "${OUTPUT_DIR}/mask.bed" > "${OUTPUT_DIR}/mask_summary.tsv"

read -r _ glen grep_bp gsq_bp gall gfrac < <(tail -1 "${OUTPUT_DIR}/mask_summary.tsv")
log "Repeat mask: ${grep_bp} bp; seq-qual mask: ${gsq_bp} bp; combined: ${gall} of ${glen} bp (${gfrac})"
log "Wrote ${OUTPUT_DIR}/mask.bed ($(wc -l < "${OUTPUT_DIR}/mask.bed") intervals), mask_sources.bed, mask_summary.tsv"
[[ ${#SQ_FILES[@]} -gt 0 ]] && log "Wrote ${OUTPUT_DIR}/seq_qual_flagged_windows.tsv ($(( $(wc -l < "${OUTPUT_DIR}/seq_qual_flagged_windows.tsv") - 1 )) windows)"
exit 0
