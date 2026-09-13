# Simple native / Wasm comparison

Use `run_comparison.sh` and the three `plot_*.py` scripts for everyday checks.
They share the 45 historical pairs in `comparison_cases.tsv` (TBLASTX, Megablast,
BLASTN and BLASTP), and write the existing `.out` / Bash `time` `.log` files.
No performance manifest, profiling phase, warm worker or certification run is
needed. `wasm_performance.py` remains available for detailed performance work.

From the repository root:

```bash
cd LOSAT
cargo build --release --bin LOSAT
cargo build --release --bin LOSAT --target wasm32-wasip1 --no-default-features --target-dir target/serial-command
cargo build --release --bin LOSAT --target wasm32-wasip1-threads --features wasm-threads --target-dir target/threaded-command
cd tests

./run_comparison.sh
./plot_overall_trend.py
./plot_comparison.py
./plot_execution_time.py
```

Requirements: Bash 4+, GNU `timeout`, Node.js with WASI support (use a supported
LTS version, such as Node 24), NCBI BLAST+ (`blastn`, `blastp`,
`tblastx`, `makeblastdb`), Python 3 with `matplotlib`, `pandas` and `seaborn>=0.12`.
`--bin LOSAT` builds the command artifact with `_start`; library Wasm artifacts
cannot run these CLI comparisons. Missing or incorrect artifacts fail explicitly.

To add only Wasm results to existing native/NCBI results:

```bash
./run_wasm_comparison.sh
./plot_overall_trend.py
./plot_comparison.py
./plot_execution_time.py
```

To run LOSATN/BLASTN and LOSATP/BLASTP without TBLASTX, select both nucleotide
tasks (`megablast` and `blastn`) and the protein task (`blastp`):

```bash
export BENCHMARK_PROGRAMS=megablast,blastn,blastp
./run_comparison.sh
```

This keeps the native, serial Wasm, threaded Wasm and NCBI runners enabled.
The exported selection also applies to the plotting scripts. For nucleotide
comparisons only, use `BENCHMARK_PROGRAMS=megablast,blastn`; for protein
comparisons only, use `BENCHMARK_PROGRAMS=blastp`. Use
`unset BENCHMARK_PROGRAMS` to restore the default selection of all programs.

To run a small selection and keep previous results intact, use the same exported
settings for execution and plotting:

```bash
export BENCHMARK_DIR="$(mktemp -d /tmp/losat-comparison.XXXXXX)"
export BENCHMARK_PROGRAMS=blastp
export BENCHMARK_CASE=WSSV.PajaWSV
export LOSAT_THREADS=4
./run_comparison.sh
./plot_overall_trend.py
./plot_comparison.py
./plot_execution_time.py
```

`BENCHMARK_PROGRAMS` accepts comma-separated `tblastx,megablast,blastn,blastp`;
`BENCHMARK_CASE` matches a substring of `losat_stem` in the table. An unmatched
selection is an error. By default all pairs run with `LOSAT_THREADS=8`.
The complete set includes long genomes; use these filters for quick checks.
Each search has a 600-second limit; set `BENCHMARK_TIMEOUT` to a larger positive
number of seconds for slower cases. A timeout stops the run and is excluded from
plots. Caught Wasm worker exceptions print their diagnostic and terminate the
command immediately (SIGTERM, shell exit 143), including while its main thread
waits inside Wasm. Both Wasm runners preserve explicit WASI exit codes. The outer
timeout still covers failures that cannot reach a JavaScript error handler.
Relative `BENCHMARK_DIR` paths resolve from this `tests` directory, even when
the scripts are launched elsewhere. Without `BENCHMARK_DIR`, selected results
in `tests/losat_out`, `tests/blast_out` and `tests/plots` are overwritten.

`RUN_LOSAT_WASM=0 ./run_comparison.sh` runs native and NCBI only.
`RUN_LOSAT_WASM_THREADED=0` disables threaded Wasm. `RUN_NATIVE=0` and
`RUN_NCBI=0` disable those runners. Plotting discovers available result files;
these runner switches do not hide previously collected series. Use a fresh
`BENCHMARK_DIR` to isolate a new run. See `./run_comparison.sh --help` for binary
overrides and optional Wasm build switches.

The script uses `node` from your current `PATH`; it does not install or pin Node.
Set `NODE_BIN=/path/to/node` to select a different executable for both serial and
threaded Wasm. At startup it prints the actual Node/V8 versions, executable path,
and Wasm file paths/modification times. When building Wasm it also prints the
Rust/Cargo versions. Rust compiles the release Wasm artifacts; Node executes them.

The threaded runner applies a shared-memory compatibility guard for all search
modes. Node 24.21.0 and 26.8.2 exposed false out-of-bounds traps after another
worker grew memory: NZ_CP006932 self BLASTN failed at `memory.fill`, and
AP027131/NZ_CP006932 BLASTP failed at `memory.copy`. V8's
[bulk-operation helpers](https://github.com/nodejs/node/blob/v24.21.0/deps/v8/src/wasm/wasm-external-refs.cc#L778-L810)
check per-instance bounds that can lag behind the actual shared memory.

Before compiling the threaded module, `wasi_shared_memory.js` guards every copy
and fill. If a range exceeds the instance's cached memory length, it executes
`memory.grow(0)` to refresh that instance, then performs the original operation.
This requests no additional pages; real invalid accesses still trap. The
on-disk Wasm artifact is unchanged, and normal optimizing compilation and
threading remain enabled. No extra flags or rebuild are needed. Adapter parsing
and compilation are included in the measured startup time. Keep `wasi_shared_memory.js`, `wasi_artifact.js`, and `wasi_thread_host.js`
beside `run_losat_wasi_threads.js` when copying the runner; set
`LOSAT_WASI_THREADS_DEBUG=1` to log the number of guarded instructions.
Unsupported module encodings fail explicitly. This guard addresses the observed
copy/fill bounds race; it does not suppress other memory faults or fix OOM.

Run the host regression checks with
`"$NODE_BIN" --test test_wasi_runners.js test_wasi_shared_memory.js`
(or substitute `node`). They exercise real Wasm commands, exit codes, worker
exceptions, shared-memory growth, overlapping copies, and invalid ranges.

Each completed search immediately prints its measured wall time, for example
`Finished WSSV.PajaWSV.losatp.wasm: 7.446 s`. Failed searches and timeouts also
print elapsed seconds and their exit status. The displayed seconds and the plot
logs use the same Bash timer.

The figures show BLAST+, native n1/nN, serial Wasm n1, and threaded Wasm nN.
With `LOSAT_THREADS=1`, native n1 runs once and the two Wasm builds stay distinct.
File naming remains compatible with previous results:

| Condition | LOSAT output/log stem |
| --- | --- |
| Native n1 | `<losat_stem>` (TBLASTX: `<losat_stem>.n1`) |
| Native nN | `<losat_stem>.nN` |
| Serial Wasm n1 | `<losat_stem>.wasm` |
| Threaded Wasm nN | `<losat_stem>.wasm.nN` |

NCBI TBLASTX retains `.nN`; other NCBI filenames have no thread suffix and run
with n1. TBLASTX uses a freshly prepared NCBI database of the same subject FASTA
(outside the timer), preserving the historical database search and honoring
`db_gencode=4`. LOSAT uses local `-subject`. This distinction follows NCBI
`blast_args.cpp:1052–1054`, which applies `db_gencode` with `-db`.

Each condition has **one wall-time sample, including process startup and Wasm
compilation**, with no warmup. Thread counts are requested counts; they are not
measurements of active workers. These plots are quick diagnostics, not release
certification or repeated-sample speed claims. `plots/execution_times.tsv`
contains the plotted seconds and source log paths. New logs record exit status;
failed or interrupted runs are excluded. Old logs without exit status can still be plotted but
cannot retrospectively prove successful execution.

Per-pair plots accept any available LOSAT series with NCBI. Overall distributions
use only the intersection of pairs available across the series present in each
program group, and report excluded pairs. Empty files count as valid zero-hit
results. A group whose paired results are all empty has no distribution to plot.
Timing plots may show partial series; absent timings are not zero-time bars.
Distribution similarity and hit counts do not establish byte parity: compare
raw `.out` files separately with `cmp`/`diff` when investigating a discrepancy.


## Threading contract gates

Use Rust 1.92.0 and Node 24.21.0 for the reproducible gates. Native, serial
command, threaded command and threaded reactor are separate artifacts. Build
all four Wasm kinds (including the internal serial reactor) in distinct target
directories; `--reverse-order` verifies the opposite build order:

```bash
python build_wasi_artifacts.py --target-dir ../target --output-dir ../../.tmp/wasi-artifacts
python check_wasm_threading.py --native ../target/release/LOSAT \
  --native-serial ../target/native-serial/release/LOSAT \
  --serial ../../.tmp/wasi-artifacts/losat-serial-command.wasm \
  --threaded ../../.tmp/wasi-artifacts/losat-threaded-command.wasm \
  --reactor ../../.tmp/wasi-artifacts/losat-threaded-reactor.wasm \
  --serial-reactor ../../.tmp/wasi-artifacts/losat-serial-reactor.wasm \
  --output-dir ../../.tmp/threading-gates
```

Build native without parallel support using `cargo build --release --bin LOSAT
--no-default-features --target-dir target/native-serial` from the crate directory.
The comparison gate requires NCBI BLAST+ 2.17.0 only as a test oracle. It retains
raw bytes, argv, stderr, artifact identity, and worker lifecycle records.

Every supported explicit N > 1 creates exactly N dedicated compute workers.
n1 creates none; serial targets reject n2/n4. Invalid or excessive requests and
invalid `LOSAT_WASI_THREAD_CAP` settings fail before spawning workers. Repeated
reactor calls join all guest workers and await host exit events. Recoverable
spawn rejection returns the underlying cause and clears old result bytes;
traps and worker failures terminate the host process and require a fresh one.

`LOSAT_WASI_THREADS_DEBUG=1` writes actual pool size and per-stage job/path
records to stderr, plus host spawn-attempt, spawned, ready and exit events.
Pool size is distinct from measured worker activity (null when unmeasured).
The detailed performance parser rejects absent or inconsistent evidence.
`LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS=1` describes sequential diagnostic chunking;
the old `LOSAT_TBLASTX_PARALLEL_SCAN_CHUNKS` name has been removed.

The CI workflow also runs `check_wasm_threading_regressions.py --jobs 3` against
all frozen manifest rows and `check_tblastx_thread_thresholds.py` for the
LC738874/LC738875 E-value 10/100/10000 sweep. `--jobs` overlaps independent
comparison processes; it does not change a search's requested thread count.
The frozen runner uses a 3600-second per-search deadline, configurable with a
positive `--timeout-seconds`; long genome cases can exceed a 900-second limit.
A deadline failure remains a failed attempt, even if a later retry succeeds.
These runs are correctness checks and must not be used for speed comparisons.
Native Gate A remains the frozen LOSAT output; the separate registered-platform
Gate B is unchanged. The local Linux oracle check is explicitly diagnostic.

`benchmark_wasm_threading.py --candidate-dir <build-root> --artifacts <wasi-dir>
--baseline-dir <baseline-build-root> --baseline-runners <saved-test-dir>
--oracle-dir <ncbi-bin-dir> --output-dir <evidence-dir>` records five workload
shapes, at least one warmup and five measured repetitions, and rejects incorrect
baseline output from timing comparisons. Cold processes, explicitly reused
compiled modules, and repeated reactor instances have separate results. Reuse
RSS is the process lifetime maximum, not an individual search's peak. Stage
selection never substitutes for measured worker activity.

BLASTP format checks cover long titles, empty query result sets, formats 0/7,
and custom fields through both streaming and rendered-alignment paths. `qacc`
and `qaccver` remain distinct columns. As in NCBI local FASTA searches, `stitle`
is `N/A` when no ASN BLAST defline object exists; the FASTA description is still
used by pairwise output.
