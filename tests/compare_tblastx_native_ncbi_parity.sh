#!/usr/bin/env bash
set -euo pipefail

# Diagnostic comparison for docs/tblastx_native_ncbi_parity_plan.md Phase 0.
# NCBI BLAST+ is used only as a validation oracle here, never by LOSAT runtime
# code. Reference-only allowance: AGENTS.md "Testing Expectations".
#
# NCBI reference comments:
# - c++/src/algo/blast/core/blast_engine.c:1482-1541 defines the ungapped
#   TBLASTX post-search order ending in preliminary E-value reaping and bit
#   score assignment.
# - c++/src/objtools/align_format/format_flags.cpp:39-40 defines the default
#   outfmt 6 fields used by the classifier.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CRATE_DIR="$ROOT_DIR/LOSAT"
SCRATCH_DIR="${1:-"$ROOT_DIR/.tmp/tblastx_native_ncbi_parity_$(date +%Y%m%d_%H%M%S)"}"
case "$SCRATCH_DIR" in
  /*) ;;
  *) SCRATCH_DIR="$ROOT_DIR/$SCRATCH_DIR" ;;
esac

QUERY="${QUERY:-"$CRATE_DIR/tests/fasta/LC738874.fasta"}"
SUBJECT="${SUBJECT:-"$CRATE_DIR/tests/fasta/LC738875.fasta"}"
OUTFMT="${OUTFMT:-6}"
THREADS="${THREADS:-1}"
NATIVE_BIN="${NATIVE_BIN:-"$CRATE_DIR/target/release/LOSAT"}"
NCBI_TBLASTX_BIN="${NCBI_TBLASTX:-tblastx}"
PYTHON_BIN="${PYTHON:-python3}"

mkdir -p "$SCRATCH_DIR"

if [[ "$OUTFMT" != "6" ]]; then
  printf 'This Phase 0 classifier expects default outfmt 6; got OUTFMT=%s\n' "$OUTFMT" >&2
  exit 2
fi

if ! command -v "$NCBI_TBLASTX_BIN" >/dev/null 2>&1; then
  printf 'NCBI tblastx not found; set NCBI_TBLASTX to the oracle binary.\n' >&2
  exit 2
fi

NCBI_TBLASTX_PATH="$(command -v "$NCBI_TBLASTX_BIN")"
printf '%s\n' "$NCBI_TBLASTX_PATH" >"$SCRATCH_DIR/ncbi_tblastx.path.txt"
"$NCBI_TBLASTX_PATH" -version >"$SCRATCH_DIR/ncbi_tblastx.version.txt" 2>&1

if [[ "${LOSAT_SKIP_BUILD:-0}" != "1" ]]; then
  (cd "$CRATE_DIR" && cargo build --release)
fi

run_native() {
  "$NATIVE_BIN" tblastx -q "$QUERY" -s "$SUBJECT" --outfmt "$OUTFMT" -n "$THREADS" -o "$1"
}

run_ncbi() {
  "$NCBI_TBLASTX_PATH" -query "$QUERY" -subject "$SUBJECT" -outfmt "$OUTFMT" -num_threads "$THREADS" -out "$1"
}

summarize_output() {
  local label="$1"
  local raw="$2"
  local sorted="$SCRATCH_DIR/${label}.sorted.out"
  LC_ALL=C sort "$raw" >"$sorted"
  {
    printf 'label\t%s\n' "$label"
    printf 'raw\t%s\t' "$raw"
    sha256sum "$raw"
    printf 'sorted\t%s\t' "$sorted"
    sha256sum "$sorted"
    printf 'line_count\t'
    wc -l <"$raw"
  } >"$SCRATCH_DIR/${label}.summary.tsv"
}

{
  printf 'query\t%s\n' "$QUERY"
  printf 'subject\t%s\n' "$SUBJECT"
  printf 'outfmt\t%s\n' "$OUTFMT"
  printf 'threads\t%s\n' "$THREADS"
  printf 'native_command\t%s tblastx -q %s -s %s --outfmt %s -n %s -o %s\n' \
    "$NATIVE_BIN" "$QUERY" "$SUBJECT" "$OUTFMT" "$THREADS" "$SCRATCH_DIR/native.raw.out"
  printf 'ncbi_command\t%s -query %s -subject %s -outfmt %s -num_threads %s -out %s\n' \
    "$NCBI_TBLASTX_PATH" "$QUERY" "$SUBJECT" "$OUTFMT" "$THREADS" "$SCRATCH_DIR/ncbi.raw.out"
} >"$SCRATCH_DIR/commands.tsv"

run_native "$SCRATCH_DIR/native.raw.out"
run_ncbi "$SCRATCH_DIR/ncbi.raw.out"

summarize_output native "$SCRATCH_DIR/native.raw.out"
summarize_output ncbi "$SCRATCH_DIR/ncbi.raw.out"

comm -3 "$SCRATCH_DIR/ncbi.sorted.out" "$SCRATCH_DIR/native.sorted.out" \
  >"$SCRATCH_DIR/ncbi_vs_native.sorted.comm3"

"$PYTHON_BIN" "$ROOT_DIR/tests/classify_tblastx_outfmt6.py" \
  --ncbi "$SCRATCH_DIR/ncbi.raw.out" \
  --native "$SCRATCH_DIR/native.raw.out" \
  --out-dir "$SCRATCH_DIR"

{
  printf 'scratch_dir\t%s\n' "$SCRATCH_DIR"
  printf 'query\t%s\n' "$QUERY"
  printf 'subject\t%s\n' "$SUBJECT"
  printf 'outfmt\t%s\n' "$OUTFMT"
  printf 'threads\t%s\n' "$THREADS"
  printf 'ncbi_tblastx\t%s\n' "$NCBI_TBLASTX_PATH"
  printf 'ncbi_vs_native_sorted_comm3_lines\t'
  wc -l <"$SCRATCH_DIR/ncbi_vs_native.sorted.comm3"
  cat "$SCRATCH_DIR/classification.summary.tsv"
} >"$SCRATCH_DIR/README.tsv"

printf 'Wrote TBLASTX native-vs-NCBI parity artifacts to %s\n' "$SCRATCH_DIR"
