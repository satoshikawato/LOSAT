# Simple native / Wasm comparison

Use `run_comparison.sh` and the three `plot_*.py` scripts for everyday checks.
They share the 45 historical pairs in `comparison_cases.tsv` (TBLASTX, Megablast,
BLASTN and BLASTP), and retain the familiar `.out` / Bash `time` `.log` filenames inside a fresh run
directory. The standard Wasm comparison uses the same threaded artifact for
n1 and nN. A small `run.json` and per-output `.run.json` bind status, argv and
file hashes to that run. `wasm_performance.py` remains available for detailed
performance work; these quick checks do not require profiling or certification.

From the repository root:

```bash
cd LOSAT
cargo build --release --bin LOSAT
cargo build --release --bin LOSAT --target wasm32-wasip1-threads --features wasm-threads --target-dir target/threaded-command
cd tests

./run_comparison.sh
./plot_overall_trend.py
./plot_comparison.py
./plot_execution_time.py
```

Requirements: Bash 4+, Node.js with WASI support (use a supported
LTS version, such as Node 24), NCBI BLAST+ (`blastn`, `blastp`,
`tblastx`, `makeblastdb`), Python 3 with `matplotlib`, `pandas` and `seaborn>=0.12`.
`--bin LOSAT` builds the command artifact with `_start`; library Wasm artifacts
cannot run these CLI comparisons. Missing or incorrect artifacts fail explicitly.

To compare Wasm with a fresh NCBI oracle without running native LOSAT:

```bash
RUN_NATIVE=0 ./run_comparison.sh
./plot_overall_trend.py
./plot_comparison.py
./plot_execution_time.py
```

`run_wasm_comparison.sh` collects Wasm-only diagnostics in a new directory.
Without a matching fresh oracle, those results are not eligible for comparison
plots. Separate invocations cannot append results to an existing run.

To run LOSATN/BLASTN and LOSATP/BLASTP without TBLASTX, select both nucleotide
tasks (`megablast` and `blastn`) and the protein task (`blastp`):

```bash
BENCHMARK_PROGRAMS=megablast,blastn,blastp ./run_comparison.sh
```

This keeps the native, threaded Wasm n1/nN and NCBI runners enabled.
The plotting scripts automatically reuse the recorded selection. For nucleotide
comparisons only, use `BENCHMARK_PROGRAMS=megablast,blastn`; for protein
comparisons only, use `BENCHMARK_PROGRAMS=blastp`. Use
`unset BENCHMARK_PROGRAMS` to restore the default selection of all programs.

To run a small selection and keep previous results intact:

```bash
BENCHMARK_PROGRAMS=blastp BENCHMARK_CASE=WSSV.PajaWSV LOSAT_THREADS=4 ./run_comparison.sh
./plot_overall_trend.py
./plot_comparison.py
./plot_execution_time.py
```

`BENCHMARK_PROGRAMS` accepts comma-separated `tblastx,megablast,blastn,blastp`;
`BENCHMARK_CASE` matches a substring of `losat_stem` in the table. An unmatched
selection is an error. By default all pairs run with `LOSAT_THREADS=8`.
The complete set includes long genomes; use these filters for quick checks.
Searches run without a time limit. Failed searches stop the run and are excluded
from plots. Caught Wasm worker exceptions print their diagnostic and terminate the
command immediately (SIGTERM, shell exit 143), including while its main thread
waits inside Wasm. Both Wasm runners preserve explicit WASI exit codes.
No directory settings are needed for normal use. Each comparison saves into a
new `tests/benchmark-runs/<timestamp>-<pid>/` directory containing `blast_out/`
and `losat_out/`; the three plotting scripts save into that run's `plots/`.
`tests/benchmark-runs/latest.json` records the latest default attempt. Plotting
automatically reads that run and its program selection, case filter, and thread
count. If it failed or is still running, plotting stops with an error instead of
silently reading an older run. Existing run files are preserved.

`BENCHMARK_DIR` is an optional override for a custom destination or for plotting
an older run. Relative paths resolve from this `tests` directory, even when the
scripts are launched elsewhere. A comparison destination must not already exist.
Explicit destinations do not change the default latest-run selection. If you
previously exported `BENCHMARK_DIR`, use `unset BENCHMARK_DIR` to return to the
automatic location. Explicit program/case/thread settings override the recorded
plotting defaults.

`RUN_LOSAT_WASM_THREADED=0 ./run_comparison.sh` runs native and NCBI only.
`RUN_LOSAT_WASM=1` adds serial compatibility checks; it is independent of
`RUN_LOSAT_WASM_THREADED`. Build that optional artifact with
`cargo build-wasi-command-serial` or set `BUILD_LOSAT_WASM=1` together with
`RUN_LOSAT_WASM=1`. `BUILD_LOSAT_WASM_THREADED=1` builds the standard threaded
artifact independently. `RUN_NATIVE=0` and
`RUN_NCBI=0` disable those runners. Plotting admits only successful outputs with matching run metadata, argv,
log and output hashes. Runner switches do not import results from older runs. See `./run_comparison.sh --help` for binary
overrides and optional Wasm build switches.

The script uses `node` from your current `PATH`; it does not install or pin Node.
Set `NODE_BIN=/path/to/node` to select a different executable for both serial and
threaded Wasm. `run_comparison.sh` and `run_wasm_comparison.sh` now add
`--no-liftoff --no-wasm-tier-up` **only for TBLASTX**. This TurboFan compilation
profile was verified on Node 24.21.0 / V8 13.6.233.17-node.53. Five measured cold
runs per condition gave 6.56 → 4.44 seconds on the primary n8 fixture; two
preserved controls gave about 2.8× speedups. Peak process RSS increased by about
28, 47 and 37 MiB respectively. The TBLASTX memory tradeoff was explicitly accepted;
BLASTN and BLASTP retain their common Node settings. See the
[follow-up evidence](../../docs/evidence/wasm_performance_20260913/run-20260913-02/REPORT.md)
for the exact fixtures, historical measurement limits and acceptance policy.

Set `NODE_TBLASTX_ARGS_JSON='[]'` to disable the TBLASTX additions, including when
collecting an ordinary-compilation baseline. `NODE_ARGS_JSON='[]'` supplies common
Node arguments for all programs; TBLASTX-specific arguments follow them. Both
variables require JSON arrays of strings. The run manifest records the effective
argv by program, and each result records the actual ordered command.

These defaults belong to the comparison scripts. Direct Node invocations must
supply the flags before the runner script, for example:

```bash
"$NODE_BIN" --no-liftoff --no-wasm-tier-up run_losat_wasi_threads.js \
  "$LOSAT_WASM_THREADED_BIN" tblastx -query fasta/MelaMJNV.fasta \
  -subject fasta/PemoMJNVA.fasta -outfmt 6 -num_threads 8 -out /tmp/tblastx.out
```

The measured improvement applies to cold Node/WASI launches. It does not establish
a browser speedup or the same gain for reused modules/reactors. Unsupported flags
fail explicitly; there is no silent fallback to a different compilation mode.
At startup the comparison script prints the actual Node/V8 versions, executable
path, and Wasm file paths/modification times. When building Wasm it also prints the
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
`Finished WSSV.PajaWSV.losatp.wasm: 7.446 s`. Failed searches also
print elapsed seconds and their exit status. The displayed seconds and the plot
logs use the same Bash timer.

The figures show BLAST+ requested n1/nN, native n1/nN and threaded Wasm n1/nN results,
plus serial Wasm n1 when compatibility comparisons are enabled.
With `LOSAT_THREADS=1`, BLAST+ and native each run once and the two Wasm builds stay distinct.
LOSAT file naming remains compatible with previous results:

| Condition | LOSAT output/log stem |
| --- | --- |
| Native n1 | `<losat_stem>` (TBLASTX: `<losat_stem>.n1`) |
| Native nN | `<losat_stem>.nN` |
| Serial Wasm n1 | `<losat_stem>.wasm` |
| Threaded Wasm nN | `<losat_stem>.wasm.nN` |

Every NCBI task writes `<ncbi_stem>.n1` and `<ncbi_stem>.nN` output/log files.
BLASTP, BLASTN and megablast use `-subject`, matching LOSAT's local-subject
search conditions. NCBI `blast_args.cpp:3223–3239` reduces these searches to
one thread even when nN is requested. Both requested runs are retained, and
the execution-time caption explicitly identifies their effective n1 behavior.

TBLASTX retains `-db` to apply non-default subject genetic codes such as
`db_gencode=4`, following NCBI `blast_args.cpp:1052–1054`. Each run prepares
one `nucl` database per TBLASTX subject FASTA with `-parse_seqids`.
Database preparation is outside the search timer; separate `.makeblastdb.log`
and `.makeblastdb.json` files record its wall time, command, input hash, and
tool identity. DBs are reused only within that run. Subject-only selections
do not require `makeblastdb`. Task and genetic-code options remain unchanged.

The TBLASTX, BLASTN and megablast fixtures each contain one subject sequence.
Database searches distribute subject OIDs to search workers
(`blast_engine.c:1411–1475`); a single subject therefore does not supply eight
independent search jobs. In TBLASTX, frames and long-subject chunks are handled
sequentially within that worker (`blast_engine.c:478–593,804–841`). The nN
label records the requested count, not a claim that all workers are busy.

Runs using the previous all-DB or single-oracle protocols cannot be plotted
with the current target/thread labels. The plot scripts report that a fresh run is required;
previous records and rendered figures remain historical evidence.

Each condition has **one wall-time sample, including process startup and Wasm
compilation**, with no warmup. Thread counts are requested counts; they are not
measurements of active workers. These plots are quick diagnostics, not release
certification or repeated-sample speed claims. `plots/execution_times.tsv`
contains the plotted seconds and source log paths. New logs record exit status;
failed or interrupted runs are excluded. Old logs without explicit successful
status and matching provenance are historical evidence only and cannot be
plotted as current results. Timing bars additionally require exact raw equality
with the same run's NCBI n1 oracle for the selected target, including lexical numeric formatting.
This check also applies to NCBI nN timing bars.

Per-pair plots accept any available LOSAT series with NCBI. Overall distributions
use only the intersection of pairs available across the series present in each
program group, and report excluded pairs. Empty files count as valid zero-hit
results. A group whose paired results are all empty has no distribution to plot.
Timing plots may show partial series; absent timings are not zero-time bars.
Distribution similarity and hit counts do not establish byte parity: compare
raw `.out` files separately with `cmp`/`diff` when investigating a discrepancy.


## Threading contract gates

Use Rust 1.92.0 and Node 24.21.0 for the reproducible gates. Native, threaded
command and threaded reactor are separate artifacts. The default helper builds
the two threaded Wasm kinds in distinct target directories; `--reverse-order`
verifies the opposite build order:

```bash
python build_wasi_artifacts.py --target-dir ../target --output-dir ../../.tmp/wasi-artifacts
python check_wasm_threading.py --native ../target/release/LOSAT \
  --threaded ../../.tmp/wasi-artifacts/losat-threaded-command.wasm \
  --reactor ../../.tmp/wasi-artifacts/losat-threaded-reactor.wasm \
  --output-dir ../../.tmp/threading-gates
```

For compatibility checks, add `--include-serial` to the builder and pass
`--serial` and `--serial-reactor` with those explicit artifacts to the checker.
To also check native without parallel support, build it using
`cargo build --release --bin LOSAT --no-default-features --target-dir target/native-serial`
from the crate directory and pass `--native-serial` to the checker.
The comparison gate requires NCBI BLAST+ 2.17.0 only as a test oracle. It retains
raw bytes, argv, stderr, artifact identity, and worker lifecycle records.

Every supported explicit N > 1 uses N total compute threads: the caller plus
N-1 child workers. Diagnostics report pool_threads=N and caller_participates=true;
host spawn/ready/exit counts are N-1. n1 creates no child workers; serial targets reject n2/n4. Invalid or excessive requests and
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

Standard artifact/CI defaults are threaded command and reactor. The reusable or
manually dispatched `wasm-threading.yml` workflow accepts
`include_serial_compatibility: true` to build and check serial compatibility.
`build_wasi_artifacts.py --include-serial` and checker `--serial` arguments are
explicit opt-ins. `benchmark_wasm_threading.py` defaults to native/threaded;
include `serial` in `--kinds` when comparing the compatibility build. Frozen
v0.1.0 certification and release contracts retain their historical target scope.
