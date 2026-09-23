#!/usr/bin/env bash
#
# Options:
#   -i PATH        input FASTQ (.fq/.fastq, optionally .gz)
#   -o TEMPLATE    output template; '{frac}' is replaced by each fraction
#                  (without '{frac}', '_<frac>' is inserted before the extension)
#   -f LIST        comma-separated fractions (0-1, or 0-100 if given as percent)
#   --from-end     sample the LATEST fraction(s) instead of earliest
#   --mode M       auto (default) | decompress | direct
#                  decompress: gunzip once to plain FASTQ, all passes read plain
#                  (ontime is ~15x faster on plain input); direct: read .gz as-is
#   --nested-below B   fractions < B use the nested cascade, >= B run in
#                  parallel from the full input (default 0.5; 0 = all parallel,
#                  1 = all nested)
#   --workdir DIR  directory for intermediates (default: input's directory)
#   --no-keep      delete the sorted-timestamp cache after the run
#   --force        overwrite existing outputs
#   -j, --threads N    threads for pigz/seqkit/sort (default: nproc)
#   -h             this help

set -euo pipefail

INPUT=""; OUT_TEMPLATE=""; FRACS_RAW=""
FROM_END=false; MODE="auto"; WORKDIR=""
KEEP=true; FORCE=false
THREADS="$(nproc 2>/dev/null || echo 4)"
NESTED_BELOW="0.5"

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }
usage() { sed -n '2,47p' "$0" | sed 's/^# \{0,1\}//'; exit 0; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) INPUT="$2"; shift 2 ;;
        -o) OUT_TEMPLATE="$2"; shift 2 ;;
        -f) FRACS_RAW="$2"; shift 2 ;;
        --from-end) FROM_END=true; shift ;;
        --mode) MODE="$2"; shift 2 ;;
        --nested-below) NESTED_BELOW="$2"; shift 2 ;;
        --workdir) WORKDIR="$2"; shift 2 ;;
        --keep-intermediates) KEEP=true; shift ;;
        --no-keep) KEEP=false; shift ;;
        --force) FORCE=true; shift ;;
        -j|--threads) THREADS="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) die "Unknown option: $1 (use -h for help)" ;;
    esac
done

[[ -n "$INPUT" && -n "$OUT_TEMPLATE" && -n "$FRACS_RAW" ]] || die "-i, -o and -f are required (use -h for help)"
[[ -f "$INPUT" ]] || die "Input not found: $INPUT"
[[ "$MODE" =~ ^(auto|decompress|direct)$ ]] || die "--mode must be auto|decompress|direct"
awk -v b="$NESTED_BELOW" 'BEGIN{exit !(b+0 >= 0 && b+0 <= 1)}' || die "--nested-below must be in [0,1]"

for tool in seqkit ontime sort sed awk wc; do
    command -v "$tool" >/dev/null 2>&1 || die "Missing dependency: $tool  (conda install -c bioconda seqkit ontime pigz)"
done
if command -v pigz >/dev/null 2>&1; then
    DCOMP=(pigz -dc -p "$THREADS"); COMP=(pigz -p "$THREADS")
    GZLIST=(pigz -l)
else
    log "NOTE: pigz not found, falling back to gzip (slower)"
    DCOMP=(gzip -dc); COMP=(gzip); GZLIST=(gzip -l)
fi

# ---- fractions (0-1 or percent; validated) ----
FRACS=()
while IFS= read -r v; do
    [[ -z "$v" ]] && continue
    norm=$(awk -v x="$v" 'BEGIN{ x=x+0; if (x > 1) x = x/100; if (x <= 0 || x >= 1) exit 1; printf "%g", x }') \
        || die "Invalid fraction: '$v' (must be in 0-1, or 0-100 if given as percent)"
    FRACS+=("$norm")
done < <(printf '%s\n' "$FRACS_RAW" | tr ',' '\n' | sed 's/ //g')
[[ ${#FRACS[@]} -gt 0 ]] || die "No valid fractions parsed from: $FRACS_RAW"
mapfile -t FRACS < <(printf '%s\n' "${FRACS[@]}" | sort -n | uniq)
log "Fractions (normalized): ${FRACS[*]}"

# ---- output paths ----
OUT_PATHS=()
for frac in "${FRACS[@]}"; do
    if [[ "$OUT_TEMPLATE" == *'{frac}'* ]]; then
        OUT_PATHS+=("${OUT_TEMPLATE//\{frac\}/$frac}")
    else
        case "$OUT_TEMPLATE" in
            *.fq.gz)     OUT_PATHS+=("${OUT_TEMPLATE%.fq.gz}_${frac}.fq.gz") ;;
            *.fastq.gz)  OUT_PATHS+=("${OUT_TEMPLATE%.fastq.gz}_${frac}.fastq.gz") ;;
            *.fq)        OUT_PATHS+=("${OUT_TEMPLATE%.fq}_${frac}.fq") ;;
            *.fastq)     OUT_PATHS+=("${OUT_TEMPLATE%.fastq}_${frac}.fastq") ;;
            *)           OUT_PATHS+=("${OUT_TEMPLATE}_${frac}") ;;
        esac
    fi
done
for p in "${OUT_PATHS[@]}"; do
    [[ -e "$p" && "$FORCE" != true ]] && die "Output exists: $p (use --force to overwrite)"
done
GZ_OUT=false; [[ "${OUT_PATHS[0]}" == *.gz ]] && GZ_OUT=true

# ---- workdir / mode ----
[[ -z "$WORKDIR" ]] && WORKDIR="$(dirname "$INPUT")"
[[ -d "$WORKDIR" && -w "$WORKDIR" ]] || die "Workdir not writable: $WORKDIR"

IS_GZ=false; [[ "$INPUT" == *.gz ]] && IS_GZ=true
PLAIN_INPUT="$INPUT"
DECOMPRESSED=""

if [[ "$MODE" == "auto" ]]; then
    if [[ "$IS_GZ" == true ]]; then
        EST_UNCOMP=$("${GZLIST[@]}" "$INPUT" 2>/dev/null | awk 'NR==2{print $2}')
        FREE=$(df -P -B1 "$WORKDIR" | awk 'NR==2{print $4}')
        if [[ -n "$EST_UNCOMP" && -n "$FREE" ]] && \
           awk -v e="$EST_UNCOMP" -v f="$FREE" 'BEGIN{exit !(f > e*1.3)}'; then
            MODE="decompress"
        else
            MODE="direct"
            log "auto mode: not enough free disk in $WORKDIR for decompressed copy -> direct"
        fi
    else
        MODE="direct"
    fi
fi

# ---- stages 1+2: timestamps (cached) ----
BASE="$(basename "${INPUT%.gz}")"
TS_SORTED="$WORKDIR/${BASE}.timestamps.sorted.txt"
TS_RAW="$WORKDIR/${BASE}.timestamps.raw.txt"

if [[ -s "$TS_SORTED" && "$TS_SORTED" -nt "$INPUT" ]]; then
    log "Reusing cached sorted timestamps: $TS_SORTED"
    if [[ "$MODE" == "decompress" && "$IS_GZ" == true ]]; then
        DECOMPRESSED="$WORKDIR/${BASE}.ontime-frac_decompressed.fq"
        if [[ -s "$DECOMPRESSED" && "$DECOMPRESSED" -nt "$INPUT" ]]; then
            log "Reusing decompressed copy: $DECOMPRESSED"
        else
            log "Decompressing input (one-time): $DECOMPRESSED"
            "${DCOMP[@]}" "$INPUT" > "$DECOMPRESSED"
        fi
        PLAIN_INPUT="$DECOMPRESSED"
    fi
else
    if [[ "$MODE" == "decompress" && "$IS_GZ" == true ]]; then
        # fused pass: one read of the .gz yields both the plain copy and timestamps
        DECOMPRESSED="$WORKDIR/${BASE}.ontime-frac_decompressed.fq"
        log "Decompressing + extracting timestamps (fused pass): $DECOMPRESSED"
        "${DCOMP[@]}" "$INPUT" \
            | tee "$DECOMPRESSED" \
            | seqkit seq -n -j "$THREADS" - 2>/dev/null \
            | sed -nE 's/.*st:Z:([^[:space:]]+).*/\1/p' > "$TS_RAW"
        PLAIN_INPUT="$DECOMPRESSED"
    else
        log "Stage 1: extracting st:Z timestamps with seqkit (-j $THREADS)"
        seqkit seq -n -j "$THREADS" "$PLAIN_INPUT" 2>/dev/null \
            | sed -nE 's/.*st:Z:([^[:space:]]+).*/\1/p' > "$TS_RAW"
    fi
    N_RAW=$(wc -l < "$TS_RAW")
    [[ "$N_RAW" -gt 0 ]] || die "No st:Z timestamps found in input"
    log "Stage 1 done: $N_RAW timestamps"
    log "Stage 2: sorting timestamps"
    sort -S 50% --parallel="$THREADS" "$TS_RAW" > "$TS_SORTED"
    rm -f "$TS_RAW"
    log "Stage 2 done: $TS_SORTED"
fi

TOTAL=$(wc -l < "$TS_SORTED")
log "Total reads with st:Z: $TOTAL"

# ---- cutoffs (absolute, from the full sorted list) ----
DIRECTION="--to";  [[ "$FROM_END" == true ]] && DIRECTION="--from"
MODE_WORD="earliest"; [[ "$FROM_END" == true ]] && MODE_WORD="latest"

declare -A CUTOFF_OF FINAL_OF
for i in "${!FRACS[@]}"; do
    frac="${FRACS[$i]}"
    IDX=$(awk -v t="$TOTAL" -v f="$frac" 'BEGIN{x=t*f; i=int(x); print (x>i)?i+1:i}')
    [[ "$IDX" -lt 1 ]] && IDX=1
    if [[ "$FROM_END" == true ]]; then POS=$(( TOTAL - IDX + 1 )); else POS=$IDX; fi
    [[ "$POS" -lt 1 ]] && POS=1; [[ "$POS" -gt "$TOTAL" ]] && POS=$TOTAL
    CUTOFF_OF[$frac]=$(sed -n "${POS}p" "$TS_SORTED")
    FINAL_OF[$frac]="${OUT_PATHS[$i]}"
    log "Cutoff [$frac]: $MODE_WORD $IDX/$TOTAL reads, st:Z ${CUTOFF_OF[$frac]}"
done

# ---- split into parallel (>= border) and nested (< border) groups ----
PAR=(); NEST=()
for frac in "${FRACS[@]}"; do
    if awk -v f="$frac" -v b="$NESTED_BELOW" 'BEGIN{exit !(f < b)}'; then
        NEST+=("$frac")
    else
        PAR+=("$frac")
    fi
done
# NEST descending (largest first)
if [[ ${#NEST[@]} -gt 0 ]]; then
    mapfile -t NEST < <(printf '%s\n' "${NEST[@]}" | sort -rn)
fi
[[ ${#PAR[@]} -gt 0 ]]  && log "Parallel group (from full input): ${PAR[*]}"
[[ ${#NEST[@]} -gt 0 ]] && log "Nested group (cascade): ${NEST[*]}"

# ---- stage 3a: parallel group ----
BRIDGE=""   # plain copy of the smallest parallel-group output, for the nested chain
PIDS=(); LOGS=()
if [[ ${#PAR[@]} -gt 0 ]]; then
    NPAR=${#PAR[@]}
    T_PER=$(awk -v t="$THREADS" -v n="$NPAR" 'BEGIN{x=int(t/n); print (x<1)?1:x}')
    # smallest PAR fraction bridges to the nested chain
    BRIDGE_FRAC=$(printf '%s\n' "${PAR[@]}" | sort -n | head -1)
    for frac in "${PAR[@]}"; do
        out="${FINAL_OF[$frac]}"
        lg="$WORKDIR/.ontime-frac_${frac}.log"; LOGS+=("$lg")
        if [[ "$GZ_OUT" == true ]]; then
            if [[ ${#NEST[@]} -gt 0 && "$frac" == "$BRIDGE_FRAC" ]]; then
                BRIDGE="$WORKDIR/.ontime-frac_bridge_${frac}.fq"
                ( ontime "$DIRECTION" "${CUTOFF_OF[$frac]}" "$PLAIN_INPUT" 2>"$lg" \
                    | tee "$BRIDGE" | pigz -p "$T_PER" -c > "$out" ) &
            else
                ( ontime "$DIRECTION" "${CUTOFF_OF[$frac]}" "$PLAIN_INPUT" 2>"$lg" \
                    | pigz -p "$T_PER" -c > "$out" ) &
            fi
        else
            ( ontime "$DIRECTION" "${CUTOFF_OF[$frac]}" "$PLAIN_INPUT" -o "$out" 2>"$lg" ) &
            if [[ ${#NEST[@]} -gt 0 && "$frac" == "$BRIDGE_FRAC" ]]; then
                BRIDGE="$out"   # plain final doubles as bridge
            fi
        fi
        PIDS+=($!)
    done
    FAILED=0
    for i in "${!PIDS[@]}"; do
        if ! wait "${PIDS[$i]}"; then
            log "Stage 3 [${PAR[$i]}]: FAILED — log follows"; cat "${LOGS[$i]}" >&2; FAILED=1
        fi
    done
    [[ "$FAILED" -eq 0 ]] || die "One or more ontime jobs failed"
    for i in "${!PAR[@]}"; do
        log "Stage 3 [${PAR[$i]}]: done -> ${FINAL_OF[${PAR[$i]}]}"
    done
fi

# ---- stage 3b: nested cascade (descending; each pass reads the previous output) ----
BGPIDS=()
if [[ ${#NEST[@]} -gt 0 ]]; then
    if [[ -n "$BRIDGE" ]]; then
        PREV="$BRIDGE"; PREV_KIND="bridge"
    else
        PREV="$PLAIN_INPUT"; PREV_KIND="input"
    fi
    LAST_IDX=$(( ${#NEST[@]} - 1 ))
    for j in "${!NEST[@]}"; do
        frac="${NEST[$j]}"
        out="${FINAL_OF[$frac]}"
        lg="$WORKDIR/.ontime-frac_${frac}.log"
        if [[ "$GZ_OUT" == true ]]; then
            if [[ $j -eq $LAST_IDX ]]; then
                # smallest nested fraction: pipe straight to final .gz, no temp
                log "Stage 3 [$frac]: ontime <- $(basename "$PREV")"
                ontime "$DIRECTION" "${CUTOFF_OF[$frac]}" "$PREV" 2>"$lg" \
                    | pigz -p "$THREADS" -c > "$out"
            else
                TMP="$WORKDIR/.ontime-frac_nest_${frac}.fq"
                log "Stage 3 [$frac]: ontime <- $(basename "$PREV")"
                ontime "$DIRECTION" "${CUTOFF_OF[$frac]}" "$PREV" -o "$TMP" 2>"$lg"
            fi
        else
            log "Stage 3 [$frac]: ontime <- $(basename "$PREV")"
            ontime "$DIRECTION" "${CUTOFF_OF[$frac]}" "$PREV" -o "$out" 2>"$lg"
        fi
        # PREV is consumed; finalize it
        if [[ "$PREV_KIND" == "bridge" && "$BRIDGE" != "${FINAL_OF[$BRIDGE_FRAC]:-}" ]]; then
            rm -f "$PREV"   # bridge temp; its .gz was written by tee
        elif [[ "$PREV_KIND" == "temp" ]]; then
            # compress previous intermediate to its final .gz in the background
            ( "${COMP[@]}" -c "$PREV" > "${FINAL_OF[$PREV_FRAC]}" && rm -f "$PREV" ) &
            BGPIDS+=($!)
        fi
        if [[ "$GZ_OUT" == true && $j -ne $LAST_IDX ]]; then
            PREV="$TMP"; PREV_KIND="temp"; PREV_FRAC="$frac"
        fi
        log "Stage 3 [$frac]: done -> $out"
    done
    for p in "${BGPIDS[@]:-}"; do [[ -n "$p" ]] && wait "$p"; done
fi

# ---- summary ----
log "Verifying output read counts..."
for frac in "${FRACS[@]}"; do
    out="${FINAL_OF[$frac]}"
    if [[ "$out" == *.gz ]]; then
        N=$("${DCOMP[@]}" "$out" | awk 'NR%4==1' | wc -l)
    else
        N=$(awk 'NR%4==1' "$out" | wc -l)
    fi
    EXPECTED=$(awk -v t="$TOTAL" -v f="$frac" 'BEGIN{x=t*f; i=int(x); print (x>i)?i+1:i}')
    log "  $out : $N reads (target $EXPECTED, fraction=$frac)"
done

rm -f "$WORKDIR"/.ontime-frac_*.log
if [[ "$KEEP" == true ]]; then
    log "Keeping timestamp cache for reuse: $TS_SORTED"
else
    rm -f "$TS_SORTED"
fi
[[ -n "$DECOMPRESSED" ]] && log "Decompressed copy retained for reuse: $DECOMPRESSED (delete manually if not needed)"
log "DONE"
