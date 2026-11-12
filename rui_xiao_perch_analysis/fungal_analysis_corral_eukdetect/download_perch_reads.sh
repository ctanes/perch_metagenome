#!/usr/bin/env bash
#
# download_perch_reads.sh (v3)
#
# v3 fixes vs v2:
#   1. "$srr_1.fastq" -> "${srr}_1.fastq"  (bash read srr_1 as a variable name;
#      under set -u every worker died instantly at that line)
#   2. Each worker runs as  ( convert_one ... || record failure )  so hard
#      crashes are captured in failed.list too
#   3. The final .sra sweep only runs after verifying ALL expected .fastq.gz
#      outputs exist - it can no longer run on a false "all good" signal
#
# Input : CSV with header  srr,sample_id   (default: ./perch_reads.csv)
# Output: /data/ruicatxiao/perch/reads/<sample_id>_1.fastq.gz
#         /data/ruicatxiao/perch/reads/<sample_id>_2.fastq.gz
#
# Usage : bash download_perch_reads.sh [path/to/perch_reads.csv]

set -euo pipefail

# ------------------------------- settings ----------------------------------
CSV="${1:-./perch_reads.csv}"
OUTDIR="/data/ruicatxiao/perch/reads"

PREFETCH_PAR=20
PREFETCH_BATCH=5
PREFETCH_MAX_KB=209715200   # 200 GB; prefetch's default -X cap is only 20G
CONV_PAR=4
DUMP_THREADS=8
PIGZ_THREADS=8
# Keep CONV_PAR * (DUMP_THREADS + PIGZ_THREADS) roughly <= your core count.

# ------------------------------ bookkeeping --------------------------------
LOGDIR="$OUTDIR/logs"
mkdir -p "$OUTDIR" "$LOGDIR"
PREFETCH_LOG="$LOGDIR/prefetch.log"
CONV_LOG="$LOGDIR/convert.log"
FAILED_LIST="$LOGDIR/failed.list"
PAIRS_TSV="$LOGDIR/srr_sampleid.tsv"
: > "$PREFETCH_LOG"; : > "$CONV_LOG"; : > "$FAILED_LIST"

# --------------------------- parse + validate ------------------------------
[[ -r "$CSV" ]] || { echo "Cannot read $CSV" >&2; exit 1; }

awk -F',' '
  NR==1 { next }
  { gsub(/\r/, "") }
  $1 != "" {
    if (NF < 2) { print "ERROR: line " NR " has no sample_id: " $0 > "/dev/stderr"; exit 1 }
    print $1 "\t" $2
  }
' "$CSV" > "$PAIRS_TSV" || { echo "Bad CSV: $CSV" >&2; exit 1; }

SRRS=(); SIDS=()
while IFS=$'\t' read -r srr sid; do
  SRRS+=("$srr"); SIDS+=("$sid")
done < "$PAIRS_TSV"

N=${#SRRS[@]}
(( N > 0 )) || { echo "No SRRs found in $CSV" >&2; exit 1; }

for s in "${SRRS[@]}"; do
  [[ "$s" =~ ^SRR[0-9]+$ ]] || { echo "Invalid SRR accession: $s" >&2; exit 1; }
done

dups=$(printf '%s\n' "${SIDS[@]}" | sort | uniq -d)
if [[ -n "$dups" ]]; then
  echo "ERROR: duplicate sample_id(s) would overwrite each other's fastqs:" >&2
  echo "$dups" >&2
  exit 1
fi

echo "Found $N SRR/sample pairs. Output dir: $OUTDIR"

# --------------------------- phase 1: prefetch ------------------------------
echo "== [1/3] Prefetching $N SRA files (parallelism $PREFETCH_PAR) ..."
printf '%s\n' "${SRRS[@]}" \
  | xargs -r -n "$PREFETCH_BATCH" -P "$PREFETCH_PAR" \
      prefetch -O "$OUTDIR" -X "$PREFETCH_MAX_KB" -p >> "$PREFETCH_LOG" 2>&1

# sra-tools 3.x + -O nests each SRA under $OUTDIR/<SRR>/<SRR>.sra.
# Normalize everything to $OUTDIR/<SRR>.sra so the rest of the pipeline works.
echo "      Normalizing SRA file locations ..."
normalize_sra() {
  local s="$1"
  local flat="$OUTDIR/$s.sra"
  [[ -s "$flat" ]] && return 0                    # already where we want it
  local found
  found=$(find "$OUTDIR" -type f -name "$s.sra" 2>/dev/null)
  [[ -z "$found" ]] && return 0                   # missing -> caught by check below
  local -a hits
  mapfile -t hits <<< "$found"
  local pick="${hits[0]}"
  if (( ${#hits[@]} > 1 )); then
    pick=$(ls -t -- "${hits[@]}" | head -n 1)     # newest wins
    echo "      WARN: multiple $s.sra found; using $(basename "$pick")" >&2
  fi
  mv -f "$pick" "$flat"
  rmdir "$(dirname "$pick")" 2>/dev/null || true  # drop the now-empty <SRR>/ dir
  echo "      moved: ${pick#"$OUTDIR"/} -> $s.sra"
}
for s in "${SRRS[@]}"; do
  normalize_sra "$s"
done

missing=0
for s in "${SRRS[@]}"; do
  if [[ ! -s "$OUTDIR/$s.sra" ]]; then
    echo "  missing/empty: $s.sra" >&2
    missing=1
  fi
done
if (( missing )); then
  echo "Prefetch incomplete - check $PREFETCH_LOG" >&2
  exit 1
fi
echo "All .sra files present in $OUTDIR"

# ---------------------- phase 2: split, rename, gzip ------------------------
echo "== [2/3] fasterq-dump -> rename to sample_id -> pigz (parallelism $CONV_PAR) ..."

export OUTDIR DUMP_THREADS PIGZ_THREADS CONV_LOG

convert_one() {
  set +e
  local srr="$1" sid="$2"
  local sra="$OUTDIR/$srr.sra"
  local g1="$OUTDIR/${sid}_1.fastq.gz"
  local g2="$OUTDIR/${sid}_2.fastq.gz"

  # Resume guard: final outputs already exist -> drop the .sra and move on
  if [[ -s "$g1" && -s "$g2" ]] || [[ -s "$OUTDIR/${sid}.fastq.gz" ]]; then
    echo "[skip] $srr ($sid) already done" >> "$CONV_LOG"
    rm -f "$sra"
    return 0
  fi

  # Clean any half-finished leftovers from a previous interrupted run
  rm -f "$OUTDIR/${srr}_1.fastq" "$OUTDIR/${srr}_2.fastq" \
        "$OUTDIR/$srr.fastq" "$OUTDIR/${sid}_1.fastq" "$OUTDIR/${sid}_2.fastq" \
        "$g1" "$g2" 2>> "$CONV_LOG"

  # fasterq-dump names outputs <SRR>_1.fastq / <SRR>_2.fastq
  if ! fasterq-dump -S -t "$DUMP_THREADS" -O "$OUTDIR" -f "$sra" >> "$CONV_LOG" 2>&1; then
    echo "[FAIL] fasterq-dump $srr" >> "$CONV_LOG"
    return 1
  fi

  local f1="$OUTDIR/${srr}_1.fastq" f2="$OUTDIR/${srr}_2.fastq"
  if [[ -s "$f1" && -s "$f2" ]]; then
    mv -f "$f1" "$OUTDIR/${sid}_1.fastq" &&
    mv -f "$f2" "$OUTDIR/${sid}_2.fastq" &&
    pigz -p "$PIGZ_THREADS" -f "$OUTDIR/${sid}_1.fastq" "$OUTDIR/${sid}_2.fastq" >> "$CONV_LOG" 2>&1
  elif [[ -s "$OUTDIR/$srr.fastq" ]]; then
    # unexpected single-end run - don't silently mix layouts
    echo "[WARN] $srr appears single-end; wrote ${sid}.fastq.gz" >> "$CONV_LOG"
    mv -f "$OUTDIR/$srr.fastq" "$OUTDIR/${sid}.fastq" &&
    pigz -p "$PIGZ_THREADS" -f "$OUTDIR/${sid}.fastq" >> "$CONV_LOG" 2>&1
  else
    echo "[FAIL] $srr: expected ${srr}_1.fastq/${srr}_2.fastq but found none" >> "$CONV_LOG"
    return 1
  fi

  if [[ -s "$g1" && -s "$g2" ]] || [[ -s "$OUTDIR/${sid}.fastq.gz" ]]; then
    rm -f "$sra"   # delete intermediate as soon as this sample is verified
    echo "[ok] $srr -> $sid" >> "$CONV_LOG"
  else
    echo "[FAIL] $srr: final .fastq.gz missing" >> "$CONV_LOG"
    return 1
  fi
}
export -f convert_one

# wait -n needs bash >= 4.3; older bash falls back to waiting on whole batches
HAS_WAIT_N=0
if (( BASH_VERSINFO[0] > 4 || (BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] >= 3) )); then
  HAS_WAIT_N=1
fi

active=0
for i in "${!SRRS[@]}"; do
  # The "|| echo" wrapper guarantees ANY failure - including hard crashes -
  # is recorded in failed.list instead of vanishing.
  ( convert_one "${SRRS[$i]}" "${SIDS[$i]}" \
      || echo "${SRRS[$i]} ${SIDS[$i]}" >> "$FAILED_LIST" ) \
      >> "$CONV_LOG" 2>&1 &
  active=$((active + 1))
  if (( active >= CONV_PAR )); then
    if (( HAS_WAIT_N )); then
      wait -n || true
      active=$((active - 1))
    else
      wait || true
      active=0
    fi
  fi
done
wait || true

# ---------------------- phase 3: verify, then cleanup -----------------------
# Defense in depth: even if failed.list were somehow empty, the .sra sweep
# only runs when every expected .fastq.gz has been verified on disk.
ok=0
for i in "${!SRRS[@]}"; do
  sid="${SIDS[$i]}"
  if [[ -s "$OUTDIR/${sid}_1.fastq.gz" && -s "$OUTDIR/${sid}_2.fastq.gz" ]] \
     || [[ -s "$OUTDIR/${sid}.fastq.gz" ]]; then
    ok=$((ok + 1))
  else
    echo "  missing outputs for ${SRRS[$i]} ($sid)" >&2
    echo "${SRRS[$i]} $sid" >> "$FAILED_LIST"
  fi
done

if (( ok != N )); then
  echo
  echo "ERROR: $((N - ok)) sample(s) have no verified .fastq.gz output - all .sra files are kept."
  echo "Failed/missing samples (srr sample_id):"
  sort -u "$FAILED_LIST" | sed 's/^/  /'
  echo "Inspect $CONV_LOG, fix the issue, and re-run: completed samples are skipped."
  exit 1
fi

echo "== [3/3] All $N samples verified - removing intermediate .sra files ..."
rm -f "$OUTDIR"/*.sra
echo "Done. $(find "$OUTDIR" -name '*.fastq.gz' | wc -l) .fastq.gz files in $OUTDIR"
echo "Disk usage: $(du -sh "$OUTDIR" | cut -f1)"