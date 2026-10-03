#!/usr/bin/env bash
# E2i sweeps of -task dc-megablast and -task blastn-short against NCBI (comparison only).
# Usage: sweeps.sh LOSAT NCBI_BIN_DIR OUTPUT_DIR WORK_DIR
# One TSV/log per command goes into OUTPUT_DIR, with the exit status appended.
set -u
N=$1
NCBI=$2
OUT=$3
WORK=$4
W=/mnt/c/Users/genom/GitHub/LOSAT-web-gui
E=$W/docs/evidence
JOBS=${JOBS:-6}
mkdir -p "$OUT" "$WORK"
step() { echo "$(date -u +%H:%M:%S) $*"; }
# The first commands run one after the other with up to 6 jobs; the batch and slice sweeps then
# run three at a time with 2 jobs each (at most 8 jobs in total).
step "scoring sweeps"
for f in 0 6 7; do
  python3 $E/losat_web_e2c/scoring_sweep.py --bin-dir "$NCBI" --losat "$N" --jobs $JOBS --outfmt $f \
    --tasks dc-megablast,blastn-short > "$OUT/scoring-sweep-fmt$f.tsv" 2>&1; echo "exit $?" >> "$OUT/scoring-sweep-fmt$f.tsv"
done
step "word size sweep"
python3 $E/losat_web_e2c/word_size_sweep.py --bin-dir "$NCBI" --losat "$N" --jobs $JOBS \
  --tasks dc-megablast,blastn-short > "$OUT/word-size-sweep.tsv" 2>&1; echo "exit $?" >> "$OUT/word-size-sweep.tsv"
step "check_inputs"
rm -rf "$WORK/inputs"
python3 $E/losat_web_e2i/check_inputs.py --bin-dir "$NCBI" --losat "$N" --work "$WORK/inputs" \
  > "$OUT/check-inputs.tsv" 2>&1; echo "exit $?" >> "$OUT/check-inputs.tsv"
step "batch sweeps (3 in parallel, 2 jobs each)"
for opts in "-task dc-megablast" "-task blastn-short" "-task blastn-short -evalue 1e-5"; do
  name=$(echo "$opts" | tr -c 'A-Za-z0-9\n' '-' | sed 's/^-*//;s/-*$//')
  (python3 $E/losat_web_e2f/batch_sweep.py --bin-dir "$NCBI" --losat "$N" --work "$WORK/batch-$name" --jobs 2 \
     --options "$opts" > "$OUT/batch-sweep-$name.tsv" 2>&1; echo "exit $?" >> "$OUT/batch-sweep-$name.tsv") &
done
wait
step "slice sweeps (3 in parallel, 2 jobs each)"
for opts in "-task dc-megablast" "-task blastn-short" "-task blastn-short -evalue 1e-5"; do
  name=$(echo "$opts" | tr -c 'A-Za-z0-9\n' '-' | sed 's/^-*//;s/-*$//')
  (python3 $E/losat_web_e2c/slice_sweep.py --bin-dir "$NCBI" --losat "$N" --work "$WORK/slice-$name" --jobs 2 \
     --options "$opts" > "$OUT/slice-sweep-$name.tsv" 2>&1; echo "exit $?" >> "$OUT/slice-sweep-$name.tsv") &
done
wait
step "done"
