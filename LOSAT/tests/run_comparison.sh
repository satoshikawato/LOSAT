#!/bin/bash
# Simple, one-run comparison using the historical .out/.log filenames.
# NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
# const string kArgQuery("query");
# const string kArgSubject("subject");
# const string kArgOutput("out");
# const string kArgNumThreads("num_threads");
# This test harness invokes NCBI only as a comparison oracle.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [[ "${1:-}" == --help ]]; then
    cat <<'HELP'
Usage: ./run_comparison.sh   (NCBI + native + serial/threaded Wasm)
       ./run_wasm_comparison.sh   (Wasm only, same cases and filenames)

Program selections:
  # LOSATN/BLASTN (both nucleotide tasks) and LOSATP/BLASTP, excluding TBLASTX
  BENCHMARK_PROGRAMS=megablast,blastn,blastp ./run_comparison.sh
  # LOSATN/BLASTN only
  BENCHMARK_PROGRAMS=megablast,blastn ./run_comparison.sh
  # LOSATP/BLASTP only
  BENCHMARK_PROGRAMS=blastp ./run_comparison.sh

Environment (also read by the three plot_*.py scripts):
  BENCHMARK_DIR       Results root containing blast_out/, losat_out/, plots/
                     (default: this tests directory)
  BENCHMARK_PROGRAMS Comma-separated tblastx,megablast,blastn,blastp (default: all)
  BENCHMARK_CASE     Substring of the LOSAT output stem (default: all pairs)
  LOSAT_THREADS     Requested multithread count (default: 8)

  BENCHMARK_TIMEOUT Maximum seconds per search (default: 600; GNU timeout)

Runner settings (0 disables, 1 enables):
  RUN_NATIVE=1 RUN_NCBI=1 RUN_LOSAT_WASM=1 RUN_LOSAT_WASM_THREADED=1
  BUILD_LOSAT_WASM=0 BUILD_LOSAT_WASM_THREADED=$BUILD_LOSAT_WASM
  LOSAT_BIN, LOSAT_WASM_BIN, LOSAT_WASM_THREADED_BIN override artifacts.
  NODE_BIN=node selects the Node executable for both Wasm targets.
  BLASTN_BIN, BLASTP_BIN, TBLASTX_BIN, MAKEBLASTDB_BIN override NCBI tools.

Each selected condition runs once. Wall time includes process startup and Wasm
compilation. Missing artifacts or failed searches stop the run. See README.md.
HELP
    exit 0
fi
[[ $# == 0 ]] || { echo "Unexpected argument: $1 (use --help)" >&2; exit 2; }

LOSAT_THREADS="${LOSAT_THREADS:-${LOSATP_THREADS:-${LOSAT_BLASTP_THREADS:-8}}}"
[[ "$LOSAT_THREADS" =~ ^[1-9][0-9]*$ ]] || { echo 'LOSAT_THREADS must be positive' >&2; exit 2; }
BENCHMARK_TIMEOUT="${BENCHMARK_TIMEOUT:-600}"
[[ "$BENCHMARK_TIMEOUT" =~ ^[1-9][0-9]*$ ]] || { echo 'BENCHMARK_TIMEOUT must be positive seconds' >&2; exit 2; }
BENCHMARK_PROGRAMS="${BENCHMARK_PROGRAMS:-tblastx,megablast,blastn,blastp}"
IFS=',' read -ra programs <<< "$BENCHMARK_PROGRAMS"
for program in "${programs[@]}"; do
    case "$program" in tblastx|megablast|blastn|blastp) ;; *) echo "Unknown program: $program" >&2; exit 2;; esac
done
BENCHMARK_DIR="${BENCHMARK_DIR:-$SCRIPT_DIR}"
mkdir -p "$BENCHMARK_DIR"
BENCHMARK_DIR="$(cd "$BENCHMARK_DIR" && pwd)"
LOSAT_OUT_DIR="$BENCHMARK_DIR/losat_out"
BLAST_OUT_DIR="$BENCHMARK_DIR/blast_out"
FASTA_DIR="$SCRIPT_DIR/fasta"
LOSAT_BIN="${LOSAT_BIN:-$SCRIPT_DIR/../target/release/LOSAT}"
LOSAT_WASM_BIN="${LOSAT_WASM_BIN:-$SCRIPT_DIR/../target/serial-command/wasm32-wasip1/release/LOSAT.wasm}"
LOSAT_WASM_THREADED_BIN="${LOSAT_WASM_THREADED_BIN:-$SCRIPT_DIR/../target/threaded-command/wasm32-wasip1-threads/release/LOSAT.wasm}"
NODE_BIN="${NODE_BIN:-node}"
RUN_NATIVE="${RUN_NATIVE:-1}"
RUN_NCBI="${RUN_NCBI:-1}"
RUN_LOSAT_WASM="${RUN_LOSAT_WASM:-1}"
RUN_LOSAT_WASM_THREADED="${RUN_LOSAT_WASM_THREADED:-1}"
BUILD_LOSAT_WASM="${BUILD_LOSAT_WASM:-0}"
BUILD_LOSAT_WASM_THREADED="${BUILD_LOSAT_WASM_THREADED:-$BUILD_LOSAT_WASM}"
for setting in RUN_NATIVE RUN_NCBI RUN_LOSAT_WASM RUN_LOSAT_WASM_THREADED BUILD_LOSAT_WASM BUILD_LOSAT_WASM_THREADED; do
    [[ "${!setting}" == 0 || "${!setting}" == 1 ]] || { echo "$setting must be 0 or 1" >&2; exit 2; }
done
[[ "$RUN_NATIVE$RUN_NCBI$RUN_LOSAT_WASM" != 000 ]] || { echo 'No runners selected' >&2; exit 2; }

selected() {
    [[ ",$BENCHMARK_PROGRAMS," == *",$task,"* && "$losat_stem" == *"${BENCHMARK_CASE:-}"* ]]
}
require_command() {
    command -v "$1" >/dev/null || { echo "Required executable unavailable: $1" >&2; exit 1; }
}
# Check inputs before running or building anything.
count=0
while IFS=$'\t' read -r task query subject name losat_stem ncbi_stem query_gencode db_gencode; do
    selected || continue
    for input in "$query" "$subject"; do
        [[ -s "$FASTA_DIR/$input" ]] || { echo "Missing FASTA: $input" >&2; exit 1; }
    done
    if [[ "$RUN_NCBI" == 1 ]]; then
        case "$task" in
            tblastx) require_command "${TBLASTX_BIN:-tblastx}"; require_command "${MAKEBLASTDB_BIN:-makeblastdb}" ;;
            blastp) require_command "${BLASTP_BIN:-blastp}" ;;
            *) require_command "${BLASTN_BIN:-blastn}" ;;
        esac
    fi
    count=$((count + 1))
done < "$SCRIPT_DIR/comparison_cases.tsv"
[[ "$count" -gt 0 ]] || { echo 'No matching comparison cases' >&2; exit 2; }
require_command timeout
if [[ "$RUN_NATIVE" == 1 ]]; then require_command "$LOSAT_BIN"; fi
if [[ "$RUN_LOSAT_WASM" == 1 ]]; then
    require_command "$NODE_BIN"
    NODE_BIN="$(command -v "$NODE_BIN")"
    "$NODE_BIN" -p '`Node: ${process.version}; V8: ${process.versions.v8}; executable: ${process.execPath}`'
    if [[ "$BUILD_LOSAT_WASM" == 1 || ( "$RUN_LOSAT_WASM_THREADED" == 1 && "$BUILD_LOSAT_WASM_THREADED" == 1 ) ]]; then
        printf 'Wasm compiler: '
        rustc --version
        printf 'Wasm build tool: '
        cargo --version
    fi
    if [[ "$BUILD_LOSAT_WASM" == 1 ]]; then
        (cd .. && cargo build --release --bin LOSAT --target wasm32-wasip1 --no-default-features --target-dir target/serial-command)
    fi
    if [[ "$RUN_LOSAT_WASM_THREADED" == 1 && "$BUILD_LOSAT_WASM_THREADED" == 1 ]]; then
        (cd .. && cargo build --release --bin LOSAT --target wasm32-wasip1-threads --features wasm-threads --target-dir target/threaded-command)
    fi
    # Require the exact requested command artifacts; never pick a stale deps/ file.
    "$NODE_BIN" - "$LOSAT_WASM_BIN" "$LOSAT_WASM_THREADED_BIN" "$RUN_LOSAT_WASM_THREADED" "$SCRIPT_DIR/wasi_artifact.js" <<'JS'
const fs = require('fs');
const [serial, threaded, enabled, inspector] = process.argv.slice(2);
const { inspectArtifact } = require(inspector);
for (const [path, threads] of [[serial, false], ...(enabled === '1' ? [[threaded, true]] : [])]) {
    inspectArtifact(fs.readFileSync(path), threads ? 'threaded-command' : 'serial-command');
    console.log(`Wasm ${threads ? 'threaded' : 'serial'}: ${path}; modified: ${fs.statSync(path).mtime.toISOString()}`);
}
JS
fi
mkdir -p "$LOSAT_OUT_DIR" "$BLAST_OUT_DIR"
export LC_ALL=C
TIMEFORMAT=$'real\t%3R\nuser\t%3U\nsys\t%3S'

# NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:48
# const string kArgOutput("out");
# Keep exit status beside timing: a failed search must never become a fast bar.
run_timed() {
    local stem="$1" status=0 elapsed="" field value
    shift
    echo "Running $(basename "$stem")"
    echo 'simple_benchmark=1' > "$stem.log"
    # Capture Bash's timer separately from command output so the displayed
    # seconds and the plotted wall time come from exactly the same measurement.
    # NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:48
    # const string kArgOutput("out");
    if { time timeout --kill-after=5s "$BENCHMARK_TIMEOUT" "$@" -out "$stem.out" >>"$stem.log" 2>&1; } 2>"$stem.time"; then
        status=0
    else
        status=$?
    fi
    cat "$stem.time" >> "$stem.log"
    while read -r field value; do
        if [[ "$field" == real ]]; then elapsed="$value"; fi
    done < "$stem.time"
    rm "$stem.time"
    echo "exit_status=$status" >> "$stem.log"
    [[ -n "$elapsed" ]] || { echo "Missing wall time; see $stem.log" >&2; exit 1; }
    if [[ "$status" == 0 ]]; then
        printf 'Finished %s: %s s\n' "${stem##*/}" "$elapsed"
    else
        printf 'Failed %s: %s s (exit %s); see %s.log\n' "${stem##*/}" "$elapsed" "$status" "$stem" >&2
        exit "$status"
    fi
}

declare -A prepared_dbs=()
while IFS=$'\t' read -r task query subject name losat_stem ncbi_stem query_gencode db_gencode; do
    selected || continue
    program="$task"
    args=(-query "$FASTA_DIR/$query" -subject "$FASTA_DIR/$subject" -outfmt 6)
    case "$task" in
        megablast|blastn) program=blastn; args+=(-task "$task") ;;
        tblastx) args+=(-query_gencode "$query_gencode" -db_gencode "$db_gencode") ;;
    esac
    single_stem="$losat_stem"
    [[ "$task" != tblastx ]] || single_stem="$losat_stem.n1"
    if [[ "$RUN_NATIVE" == 1 ]]; then
        run_timed "$LOSAT_OUT_DIR/$single_stem" "$LOSAT_BIN" "$program" "${args[@]}" -num_threads 1
        if [[ "$LOSAT_THREADS" != 1 ]]; then
            run_timed "$LOSAT_OUT_DIR/$losat_stem.n$LOSAT_THREADS" "$LOSAT_BIN" "$program" "${args[@]}" -num_threads "$LOSAT_THREADS"
        fi
    fi
    if [[ "$RUN_LOSAT_WASM" == 1 ]]; then
        run_timed "$LOSAT_OUT_DIR/$losat_stem.wasm" env NODE_NO_WARNINGS=1 "$NODE_BIN" "$SCRIPT_DIR/run_losat_wasi.js" "$LOSAT_WASM_BIN" "$program" "${args[@]}" -num_threads 1
        if [[ "$RUN_LOSAT_WASM_THREADED" == 1 ]]; then
            run_timed "$LOSAT_OUT_DIR/$losat_stem.wasm.n$LOSAT_THREADS" env NODE_NO_WARNINGS=1 "$NODE_BIN" "$SCRIPT_DIR/run_losat_wasi_threads.js" "$LOSAT_WASM_THREADED_BIN" "$program" "${args[@]}" -num_threads "$LOSAT_THREADS"
        fi
    fi
    if [[ "$RUN_NCBI" == 1 ]]; then
        ncbi_threads=1
        case "$task" in
            tblastx)
                ncbi_bin="${TBLASTX_BIN:-tblastx}"
                ncbi_threads="$LOSAT_THREADS"
                ncbi_stem="$ncbi_stem.n$LOSAT_THREADS"
                # NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1052-1054
                # if (m_Target == eDatabase && args[kArgDbGeneticCode] &&
                #     args.Exist(kArgDb) && args[kArgDb]) {
                #     opt.SetDbGeneticCode(args[kArgDbGeneticCode].AsInteger());
                # Preserve the historical database oracle so code 4 is applied.
                db="$BLAST_OUT_DIR/db/${subject%.*}"
                if [[ -z "${prepared_dbs[$subject]:-}" ]]; then
                    mkdir -p "$BLAST_OUT_DIR/db"
                    if ! "${MAKEBLASTDB_BIN:-makeblastdb}" -in "$FASTA_DIR/$subject" -dbtype nucl -parse_seqids -out "$db" >"$db.makeblastdb.log" 2>&1; then
                        echo "Database preparation failed; see $db.makeblastdb.log" >&2
                        exit 1
                    fi
                    prepared_dbs[$subject]=1
                fi
                args=(-query "$FASTA_DIR/$query" -db "$db" -outfmt 6 -query_gencode "$query_gencode" -db_gencode "$db_gencode")
                ;;
            blastp) ncbi_bin="${BLASTP_BIN:-blastp}" ;;
            *) ncbi_bin="${BLASTN_BIN:-blastn}" ;;
        esac
        run_timed "$BLAST_OUT_DIR/$ncbi_stem" "$ncbi_bin" "${args[@]}" -num_threads "$ncbi_threads"
    fi
done < "$SCRIPT_DIR/comparison_cases.tsv"
echo "Finished $count pairs. Results: $BENCHMARK_DIR"
