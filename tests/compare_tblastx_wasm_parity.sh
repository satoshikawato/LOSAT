#!/usr/bin/env bash
set -euo pipefail

# Diagnostic comparison for docs/tblastx_wasm_parity_plan.md.
# NCBI BLAST+ is used only as a validation oracle here, never by LOSAT runtime
# code. Reference-only allowance: AGENTS.md "Testing Expectations"; NCBI
# command-line options correspond to cmdline_flags.cpp:46-94.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CRATE_DIR="$ROOT_DIR/LOSAT"
SCRATCH_DIR="${1:-"$ROOT_DIR/.tmp/tblastx_wasm_parity_$(date +%Y%m%d_%H%M%S)"}"
case "$SCRATCH_DIR" in
  /*) ;;
  *) SCRATCH_DIR="$ROOT_DIR/$SCRATCH_DIR" ;;
esac

QUERY="$CRATE_DIR/tests/fasta/LC738874.fasta"
SUBJECT="$CRATE_DIR/tests/fasta/LC738875.fasta"
OUTFMT="6"
THREADS="1"
ALLOW_LOSAT_INTERNAL_DIFFS="${ALLOW_LOSAT_INTERNAL_DIFFS:-0}"
REQUIRE_NCBI_TBLASTX_PARITY="${REQUIRE_NCBI_TBLASTX_PARITY:-0}"

NATIVE_BIN="$CRATE_DIR/target/release/LOSAT"
WASM_BIN="$CRATE_DIR/target/wasm32-wasip1/release/LOSAT.wasm"
NODE_RUNNER="$SCRATCH_DIR/run_wasi.js"

mkdir -p "$SCRATCH_DIR"

cat >"$NODE_RUNNER" <<'JS'
const { WASI } = require('wasi');
const fs = require('fs');

const wasmPath = process.argv[2];
const args = process.argv.slice(3);
const wasi = new WASI({
  version: 'preview1',
  args: [wasmPath, ...args],
  env: process.env,
  preopens: { '/': '/' },
});

(async () => {
  const bytes = fs.readFileSync(wasmPath);
  const module = await WebAssembly.compile(bytes);
  const instance = await WebAssembly.instantiate(module, {
    wasi_snapshot_preview1: wasi.wasiImport,
  });
  wasi.start(instance);
})().catch((err) => {
  console.error(err);
  process.exit(1);
});
JS

run_native() {
  "$NATIVE_BIN" tblastx -q "$QUERY" -s "$SUBJECT" --outfmt "$OUTFMT" -n "$THREADS" -o "$1"
}

run_wasm() {
  node "$NODE_RUNNER" "$WASM_BIN" tblastx -q "$QUERY" -s "$SUBJECT" --outfmt "$OUTFMT" -n "$THREADS" -o "$1"
}

run_ncbi() {
  local tblastx_bin="${NCBI_TBLASTX:-tblastx}"
  "$tblastx_bin" -query "$QUERY" -subject "$SUBJECT" -outfmt "$OUTFMT" -num_threads "$THREADS" -out "$1"
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

compare_sorted() {
  local left="$1"
  local right="$2"
  comm -3 "$SCRATCH_DIR/${left}.sorted.out" "$SCRATCH_DIR/${right}.sorted.out" \
    >"$SCRATCH_DIR/${left}_vs_${right}.comm3"
}

# NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
# ```c
# while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
#     if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
#         continue;
#     }
# }
# ```
compare_raw() {
  local left="$1"
  local right="$2"
  diff -u "$SCRATCH_DIR/${left}.raw.out" "$SCRATCH_DIR/${right}.raw.out" \
    >"$SCRATCH_DIR/${left}_vs_${right}.raw.diff" || true
}

(
  cd "$CRATE_DIR"
  cargo build --release
  cargo build --release --target wasm32-wasip1 --no-default-features
)
run_native "$SCRATCH_DIR/native.raw.out"
run_wasm "$SCRATCH_DIR/wasm_simd.raw.out"

(
  cd "$CRATE_DIR"
  cargo build --release --target wasm32-wasip1 --no-default-features --features tblastx-wasm-scalar
)
run_wasm "$SCRATCH_DIR/wasm_scalar.raw.out"

if command -v "${NCBI_TBLASTX:-tblastx}" >/dev/null 2>&1; then
  run_ncbi "$SCRATCH_DIR/ncbi.raw.out"
else
  printf 'NCBI tblastx not found; set NCBI_TBLASTX to capture the oracle output.\n' \
    >"$SCRATCH_DIR/ncbi.SKIPPED"
fi

for label in native wasm_simd wasm_scalar; do
  summarize_output "$label" "$SCRATCH_DIR/${label}.raw.out"
done
if [[ -f "$SCRATCH_DIR/ncbi.raw.out" ]]; then
  summarize_output ncbi "$SCRATCH_DIR/ncbi.raw.out"
fi

compare_sorted native wasm_simd
compare_sorted native wasm_scalar
compare_sorted wasm_scalar wasm_simd
compare_raw native wasm_simd
compare_raw native wasm_scalar
compare_raw wasm_scalar wasm_simd
if [[ -f "$SCRATCH_DIR/ncbi.raw.out" ]]; then
  compare_sorted ncbi native
  compare_sorted ncbi wasm_simd
  compare_sorted ncbi wasm_scalar
  compare_raw ncbi native
  compare_raw ncbi wasm_simd
  compare_raw ncbi wasm_scalar
fi

{
  printf 'scratch_dir\t%s\n' "$SCRATCH_DIR"
  printf 'query\t%s\n' "$QUERY"
  printf 'subject\t%s\n' "$SUBJECT"
  printf 'outfmt\t%s\n' "$OUTFMT"
  printf 'threads\t%s\n' "$THREADS"
  printf 'native_vs_wasm_simd_diff_lines\t'
  wc -l <"$SCRATCH_DIR/native_vs_wasm_simd.comm3"
  printf 'native_vs_wasm_scalar_diff_lines\t'
  wc -l <"$SCRATCH_DIR/native_vs_wasm_scalar.comm3"
  printf 'wasm_scalar_vs_wasm_simd_diff_lines\t'
  wc -l <"$SCRATCH_DIR/wasm_scalar_vs_wasm_simd.comm3"
  printf 'native_vs_wasm_simd_raw_diff_lines\t'
  wc -l <"$SCRATCH_DIR/native_vs_wasm_simd.raw.diff"
  printf 'native_vs_wasm_scalar_raw_diff_lines\t'
  wc -l <"$SCRATCH_DIR/native_vs_wasm_scalar.raw.diff"
  printf 'wasm_scalar_vs_wasm_simd_raw_diff_lines\t'
  wc -l <"$SCRATCH_DIR/wasm_scalar_vs_wasm_simd.raw.diff"
  if [[ -f "$SCRATCH_DIR/ncbi.raw.out" ]]; then
    printf 'ncbi_vs_native_raw_diff_lines\t'
    wc -l <"$SCRATCH_DIR/ncbi_vs_native.raw.diff"
    printf 'ncbi_vs_wasm_simd_raw_diff_lines\t'
    wc -l <"$SCRATCH_DIR/ncbi_vs_wasm_simd.raw.diff"
    printf 'ncbi_vs_wasm_scalar_raw_diff_lines\t'
    wc -l <"$SCRATCH_DIR/ncbi_vs_wasm_scalar.raw.diff"
  fi
} >"$SCRATCH_DIR/README.tsv"

losat_internal_raw_diff_lines=$(
  awk -F'\t' '
    $1 == "native_vs_wasm_simd_raw_diff_lines" ||
    $1 == "native_vs_wasm_scalar_raw_diff_lines" ||
    $1 == "wasm_scalar_vs_wasm_simd_raw_diff_lines" { total += $2 }
    END { print total + 0 }
  ' "$SCRATCH_DIR/README.tsv"
)

if [[ "$losat_internal_raw_diff_lines" -ne 0 && "$ALLOW_LOSAT_INTERNAL_DIFFS" != "1" ]]; then
  printf 'LOSAT native/Wasm scalar/Wasm SIMD raw outputs differ; inspect %s/*.raw.diff\n' \
    "$SCRATCH_DIR" >&2
  exit 1
fi

if [[ "$REQUIRE_NCBI_TBLASTX_PARITY" == "1" && -f "$SCRATCH_DIR/ncbi.raw.out" ]]; then
  ncbi_raw_diff_lines=$(
    awk -F'\t' '
      $1 == "ncbi_vs_native_raw_diff_lines" ||
      $1 == "ncbi_vs_wasm_simd_raw_diff_lines" ||
      $1 == "ncbi_vs_wasm_scalar_raw_diff_lines" { total += $2 }
      END { print total + 0 }
    ' "$SCRATCH_DIR/README.tsv"
  )
  if [[ "$ncbi_raw_diff_lines" -ne 0 ]]; then
    printf 'NCBI oracle output differs; inspect %s/*.raw.diff\n' "$SCRATCH_DIR" >&2
    exit 1
  fi
fi

printf 'Wrote TBLASTX Wasm parity artifacts to %s\n' "$SCRATCH_DIR"
