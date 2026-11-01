# LOSAT debug and diagnostics environment variables

Moved from `AGENTS.md` on 2026-10-08 so that it is read only when needed.

## Floating-point contract

The reports match NCBI's bytes only while every floating-point value is computed as NCBI's C
code computes it on the same machine. Every build keeps these rules (written down 2026-10-18,
plan v4 D6 (c)); `LOSAT/tests/check_fp_contract.py` checks them in CI: `source` in the `rust`
job and in the release workflow's quality job, `binary` on the release binary of the
`fast output regressions` job, on the Linux x64 release artifact, and on the release binary of
the Linux arm64 cell of `platform-checks.yml` (since 2026-10-18, plan v4 D6 (a)(b)).

- **No contraction.** rustc emits each `*`, `+`, `-` and `/` as its own IEEE operation, rounded
  once, and never fuses `a * b + c` into a fused multiply-add: it sets no `contract` or fast-math
  flag, whatever the target CPU or `#[target_feature]`. An FMA-capable target (x86-64-v3,
  aarch64, the `avx512f` helper of `x_newton_exact.rs`) changes only how `mul_add` itself is
  computed. Do not add fast-math LLVM options (`-C llvm-args=-fp-contract=fast`,
  `-enable-unsafe-fp-math`, `-enable-no-nans-fp-math` and the like), the `*_fast` and
  `*_algebraic` operations, FMA intrinsics (`_mm*_fmadd_*` and the other `_mm*_fm*` /
  `_mm*_fnm*`, `vfmaq_*`, `vfmsq_*`), or the Wasm `relaxed-simd` feature (its `relaxed_madd`
  fuses or not depending on the engine). The Wasm builds enable `+simd128` only.
- **One fused multiply-add**, `f64::mul_add` in `stats/spouge.rs::x_product_error`: the exact
  rounding error of a product in the double-double emulation of `expl` in the long-double line
  `result *= expl(-z * z) / z` of NCBI's `boost_erf.c`. The other call is its test
  `x_dekker_product_error_matches_fused_multiply_add`. `mul_add` rounds once (IEEE 754
  fusedMultiplyAdd) on every target, so an instruction and a library routine give the same bits.
  In the release binary of dev 79fab14e it sits in `stats::spouge::erf_impl_ncbi`:
  - on baseline x86-64, as calls (13 call sites, through one GOT slot) of the `fma` of Rust's
    compiler-builtins, not glibc's. That `fma` takes `vfmadd213sd` when the CPU has FMA, the FMA4
    form when it has only FMA4, and a software fused multiply-add otherwise; these two
    instructions are the binary's only FMA instructions;
  - with `-C target-cpu=x86-64-v3`, as 13 inline `vfmsub213sd`, and no other function holds an
    FMA instruction; the same with PGO (`-Cprofile-use`, dev 0b928fff);
  - on aarch64, as 13 inline `fnmsub` (a cross build for `aarch64-unknown-linux-musl` of dev
    0b928fff, counted with `llvm-objdump`; that static build also holds the 7 FMA instructions
    of musl's `log`, `log2` and `pow`, which the glibc build of the Linux arm64 cell calls in
    the shared libm instead).

  PGO and other inlining changes may move the site into a caller outside `stats::spouge`; the
  check then fails and names the function.

  A new `mul_add` needs the same argument (the fused result is the value NCBI computes) and an
  entry in the checker's `MUL_ADD_SITES` (and, in another module, `MUL_ADD_FUNCTION_PREFIXES`).
- **The C library's math functions.** `ln`, `exp`, `powf` and the other transcendental functions
  call the platform's libm (glibc on Linux, Apple's libm on macOS, the UCRT on Windows, wasi-libc
  in Wasm), as NCBI's C code does, so on one machine both get the same values. No build option
  changes which function is called: glibc picks its variant (for example `__log_fma` on an
  FMA-capable CPU) when it is loaded, from the CPU, not from LOSAT's target CPU. Do not
  replace these calls by another implementation (the `libm` crate in `LOSAT/Cargo.toml` is not
  called; the release binary holds none of its functions).
- **Release builds.** `cargo build --release`: fat LTO, one codegen unit, `panic = "abort"`
  (`LOSAT/Cargo.toml`) and the x86-64 loop alignment of `LOSAT/.cargo/config.toml`, for the
  baseline CPU of the target; the AVX2 kernels are chosen at run time
  (`is_x86_feature_detected!`). Extra flags for one build go in `CARGO_TARGET_<TRIPLE>_RUSTFLAGS`,
  which cargo joins with the configuration's `[target.'cfg(target_arch = "x86_64")']` flags;
  `RUSTFLAGS` replaces them. Measured 2026-10-18 (plan v4 D6 (a)(b); ws, 11 small bench5 jobs
  of all five programs at 1 and 8 threads, 8 rounds; every output identical): a build for
  `x86-64-v3` takes ×1.00 of the baseline's time at 1 thread and ×1.01 at 8 (TBLASTX ×1.01 to
  ×1.035: 2.9 % fewer instructions, 3 % more cycles), and PGO (`-Cprofile-use`, profile of 8 of
  the jobs) ×0.98 (best ×0.93, megablast at 1 thread; a 0.06 s BLASTN job ×1.03 to ×1.07).
  The release builds stay at the baseline CPU, without PGO.

## Debug and trace variables

### Diagnostic build

The trace, debug, dump and profile variables that no test, script or CI step sets are read only by
a build with the `diagnostics` cargo feature (`cargo build --release --features diagnostics`; the WASI
command adds `--target wasm32-wasip1 --no-default-features`). A build without it ignores them and
writes one stderr line for each that is set (`utils/xdiag.rs`). They take the values their entries below
give.

The variables of the diagnostic build (the entries marked "(diagnostic build)" below):

```text
LOSAT_DEBUG_BLASTP_TRACE_HSP
LOSAT_DEBUG_BLASTP_WSSV
LOSAT_DEBUG_CHAINING
LOSAT_DEBUG_CUTOFFS_ALL
LOSAT_DEBUG_EXTENSION
LOSAT_DEBUG_HSP_SAVING
LOSAT_DEBUG_SCAN_SOFF
LOSAT_DEBUG_SCAN_WINDOW
LOSAT_DUMP_TBLASTX_STAGE
LOSAT_LINK_WIDE_KEYS
LOSAT_STARTUP_TRACE
LOSAT_TRACE_CHAIN_HSP
LOSAT_TRACE_HSP
LOSAT_TRACE_HSP_MASKS
LOSAT_TRACE_LINK_SELECTIONS
LOSAT_X_ADJMEMO_STATS
LOSAT_X_DPCAPTURE
LOSAT_X_DPCAPTURE_EVERY
LOSAT_X_NEWTONEXACT_PROF
```

### All debug and trace variables

- `LOSAT_TRACE_HSP="qstart,qend,sstart,send"` (diagnostic build) trace a specific TBLASTX HSP.
- `LOSAT_TRACE_HSP_MASKS=1` (diagnostic build) print mask coverage for the traced TBLASTX HSP.
- `LOSAT_TRACE_CHAIN_HSP="qstart,qend,sstart,send"` (diagnostic build) trace TBLASTX chain
  selection for a specific HSP.
- `LOSAT_TRACE_LINK_SELECTIONS=1` (diagnostic build) print TBLASTX link-selection details.
- `LOSAT_LINKING_LEGACY=1` TBLASTX sum-statistics linking by the literal port
  `link_hsp_group_ncbi` instead of the default kernel (the incremental kernel of
  `tblastx/sum_stats_linking/linking_incr.rs`; same result). The linking traces
  (`LOSAT_TRACE_HSP`, `LOSAT_TRACE_CHAIN_HSP`, `LOSAT_TRACE_LINK_SELECTIONS`,
  `LOSAT_DEBUG_CHAINING`; all diagnostic build) also run the literal port, which prints them.
- Checks of the default TBLASTX linking kernel: `LOSAT_LINK_FAST_VERIFY=1` (every
  pass against a pass computed from scratch by NCBI's scans; slow on large
  groups), `LOSAT_LINK_FAST_SHADOW=1` (every group also linked by the literal
  port, all fields compared), `LOSAT_LINK_WIDE_KEYS=1` (diagnostic build; u128 keys for every
  group, to measure their cost), `LOSAT_LINK_STATS=1` (per-group counters,
  `[LINK_INCR_STATS]` lines).
- `LOSAT_DUMP_TBLASTX_STAGE=<dir>` (diagnostic build) append TBLASTX stage snapshots as TSV files.
- `LOSAT_TRACE_BLASTN_HSP="qstart,qend,sstart,send"` trace a specific BLASTN
  HSP.
- `LOSAT_TRACE_BLASTN_SEED="q,s"` trace a specific BLASTN seed.
- `LOSAT_TRACE_BLASTN_CONTEXT=<context_idx>` restrict BLASTN tracing by context.
- `LOSAT_TRACE_BLASTN_SUBJECT=<subject_id_or_index>` restrict BLASTN tracing by
  subject.
- `LOSAT_TRACE_BLASTN_STAGE=<seed|ungapped|prelim|traceback|purge|hitlist|all>`
  restrict BLASTN tracing by stage.
- `LOSAT_DEBUG_CUTOFFS=1` cutoff calculations (tblastx + blastn).
- `LOSAT_DEBUG_CUTOFFS_ALL=1` (diagnostic build) verbose TBLASTX cutoff diagnostics.
- `LOSAT_DEBUG_CHAINING=1` (diagnostic build) chaining debug (legacy; tblastx).
- `LOSAT_DEBUG_EXTENSION=1` (diagnostic build) tblastx extension debug.
- `LOSAT_DEBUG_HSP_SAVING=1` (diagnostic build) TBLASTX HSP-save diagnostics.
- `LOSAT_DEBUG_OUTPUT_FILTER=1` TBLASTX output filter diagnostics.
- `LOSAT_DEBUG_BLASTN=1` blastn hit loss diagnostics.
- `LOSAT_DEBUG_COORDS=1` blastn coordinate transforms.
- `LOSAT_DEBUG_COORDS_START=<int>` narrow selected BLASTN coordinate diagnostics.
- `LOSAT_DEBUG_SCAN_SOFF=<int>` (diagnostic build) tblastx scan debug center subject offset.
- `LOSAT_DEBUG_SCAN_WINDOW=<int>` (diagnostic build) tblastx scan debug window size.
- `LOSAT_DEBUG_BLASTP_TRACE_HSP="[qidx,sidx,]qstart,qend,sstart,send"` (diagnostic build) trace one
  BLASTP HSP through the composition-based redo (4 or 6 comma-separated fields).
- `LOSAT_DEBUG_BLASTP_WSSV=1` (diagnostic build) print the redo's `[BLASTP_WSSV]` line for the one HSP
  q 78-789 x s 0-665.
- `LOSAT_TIMING=1` timing breakdown. Every search program prints one line last, after its own lines
  (`utils/xtiming.rs`):
  `[TIMING] common: since_start .. outside_search .. read_fasta .. encode .. translate .. seg .. build_lookup .. search .. output .. (search_sum ..)`.
  `since_start` is the wall time from the entry of `main` to the line; `search` is the wall time during
  which at least one of the program's search intervals ran (batches searched side by side count once),
  `search_sum` adds those intervals up, and `outside_search = since_start - search` is the run's fixed cost.
  The other fields add up the intervals the program's driver marks (on several threads at once their wall
  times add up): reading the FASTA query and subject files, encoding, translating the query (BLASTX,
  TBLASTX), the query's SEG (DUST for BLASTN), building the lookup tables, writing the reports. Work done
  inside the search on its threads (per-subject translation, the lookup tables of BLASTX's parallel chunks and
  of TBLASTN's subject side) is in `search`; BLASTX's `read_fasta` is the `input` of its `blastx_run` line
  (`loop - search - render`), and BLASTN's streamed reports (`LOSAT_X_RPTSTREAM`) are in `search` too.
  A field a program does not have stays 0.
  TBLASTX also prints `[TIMING] x_unit: unit_sum .. seed .. extend .. (units=N hits=N extends=N seed_ns_per_hit .. extend_ns_per_call .. clock_overhead .. table_kib_mean .. table_kib_max ..)`
  after `x_pairpar`: the time of the default unit path's units (`unit_sum`, added over the threads that ran
  them) split into the seed stage (table reset, scan, two-hit test, saving init HSPs) and the two-hit
  extensions, with two reads of the time-stamp counter around each extension (x86-64; `Instant` elsewhere)
  whose measured cost is taken out of both parts (`clock_overhead`); an estimate
  (`tblastx/x_unit_probe.rs`). Without the variable the units run without the probe.
  TBLASTX prints, per query batch, after its other lines:
  `[TIMING] x_search: setup .. subjects .. collect ..` (three wall times that add up to
  `search_total`: before the first subject, the subject loop, after it) and
  `[TIMING] x_subject: enter .. translate .. prepare .. units .. units_release .. merge .. frame_loop .. seed_release .. link1 .. reap .. reevaluate .. link2 .. convert .. emit .. end .. sum .. (subjects=N)`
  (the stages of the per-subject work in the order they run; one clock per subject is read at
  every stage boundary, so the stages add up to `sum`, and `end` includes the release of the
  subject's translation; with several subjects the times are sums over the subjects and their
  threads). `link1` is the first linking call, on every HSP the seed stage saved; `link2` is the
  call `sum_stats_linking` times (`tblastx/x_search_timing.rs`). Two more lines follow, one per
  linking call:
  `[TIMING] x_link1: reset .. sort .. groups .. kernels .. assemble .. final_sorts .. replay .. sum .. (calls=N hsps=N groups=N largest_group=N kernel_sum ..s kernel_max ..s)`
  and `[TIMING] x_link2:` in the same form (the steps of the default linking call,
  `x_hsp_inplace::link_in_place`, in the order they run; one clock per call, so the steps add up
  to `sum`. `kernels` is the wall time of the per-group linking, `kernel_sum` the per-group
  times added over the threads that ran them, `kernel_max` the longest group, which `kernels`
  cannot go below; `largest_group` is the number of HSPs in the largest group; with several subjects
  `largest_group` and `kernel_max` are the largest over the calls, the other values sums. A linking call
  that takes the reference body (`LOSAT_DIAGNOSTICS`, `LOSAT_DEBUG_CUTOFFS`) is not timed there:
  `calls=0`;
  `tblastx/x_link_timing.rs`). The clock is not read without the variable.
  BLASTX prints six kinds of lines at the end of a run (`blastx/x_timing.rs`):
  `[TIMING] blastx_run: startup .. prolog .. read_stores .. loop .. epilog .. total .. (search .. render .. titles .. input ..)`
  (wall times of the thread that runs the program; inside `loop`, the time it waits for searches, writes reports (of
  which `titles` reads the reported sequences' titles with `-db`), and the rest, `input`);
  `[TIMING] blastx_waves: waves=N batches=N wall .. busy .. threads=N idle .. tail=N tail_wall ..`
  (at 2+ threads with several query batches: the waves of query batches searched side by side, the searches' own times added over the
  threads, `idle = wall x threads - busy`, and the batches searched one at a time with the whole pool: the end of the
  input, or a long query);
  `[TIMING] blastx_subject_waves: waves=N tasks=N wall .. busy .. longest .. latest_start .. last_end .. threads=N idle ..`
  (the parallel subject waves of a query chunk, the `scan_par` stage: a task is the subjects one worker searches in a
  wave; `busy` adds the tasks' times, `longest` each wave's longest task, `latest_start` / `last_end` each wave's
  latest task start / end after the wave began: `longest` near `wall` with `busy` well below `wall x threads` is uneven
  work);
  `[TIMING] blastx_subject_flow: loops=N tasks=N ahead=N taken=N wall ..` (`LOSAT_X_BXSUBJFLOW`, below: the
  subject loops run as a flow, their tasks, those searched ahead of the loop (by a helper, or by the loop's thread while
  it waited for a helper's task) and those of them the loop took (it searched the others itself at their turn), and
  the wall time, which is the batch line's `scan_par`, the hand-over to the sink included);
  `[TIMING] blastx_batch_side:` / `blastx_batch_alone: batches=N query_setup .. setup .. chunks_par .. lookup .. scan_par .. scan_inline .. replay .. chunk_end .. batch_end .. sum ..`
  (the stages of a query batch in the order they run; `side` = batches searched on one thread beside others, times
  added over the threads; `alone` = a batch searched by itself, with the pool if there is one, wall times; `replay`
  hands the subjects to the sink in subject order and includes `scan_inline`) and
  `[TIMING] blastx_batch_end_side:` / `blastx_batch_end_alone: init .. read .. kappa_spec .. kappa_serial .. loop .. finish .. (matches=N reused=N serial=N)`
  (inside the end of a batch: the composition-based redo made ahead on the pool, the redo done in the in-order loop,
  and the hit lists after it).
- `LOSAT_DIAGNOSTICS=1` general diagnostics counters.
- `LOSAT_X_STATS=1` with TBLASTX's default unit path: `[X_STATS] TBXSEED` lines at exit
  (`tblastx/x_unit_probe.rs`; with `LOSAT_TIMING` too, these take the place of the `x_unit` split). Per cell
  of a unit's diagonal table (one diagonal of the concatenated query): the outcomes of the two-hit test
  (`masked flag_reset too_far overlap context_boundary extended saved`, summing to `hits`); the hits per
  cell over the whole unit and over each `scansub` batch (`hbar` = hits / cells hit, histograms 1 / 2 / 3 /
  4-7 / 8+ by cells and by hits); the isolated hits (no other hit on the cell less than the window before
  or after); the runs of consecutive `masked` and of consecutive `too_far` outcomes on a cell (lengths 1 /
  2 / 3 / 4-7 / 8-15 / 16-63 / 64+); the extensions' cells (word, left, right; counted by repeating each
  extension, `extension/two_hit.rs::x_two_hit_cells`) and `rescan` = cells / distinct subject positions
  per cell (positions past the furthest one an earlier extension of the cell reached). The statistics
  repeat every extension, so the run is slower; reports do not change.
- `LOSAT_X_REACH=<directory>` records, at the decision points of the conditional switches of the strict set
  (`utils/xreach.rs`, `x_reach!`), whether each took its own path or the reference path while it was on, and
  writes `reach-<pid>-<ns>.tsv` (`switch new ref program`) into the directory at exit; output unchanged.
  `LOSAT/tests/switch_reach_check.py` compares the folded files with `LOSAT/tests/switch_reach_expected.tsv`
  (one row per switch of the set and per graduated one, `GRADUATED`: class `always`, `platform` or
  `conditional`, never `always` for a graduated one; for a conditional switch each
  side `required` or `gap` with the input that would reach it). CI runs the fast regressions with it and the
  check after them. A new switch of the set adds its row; a conditional one also marks its decision points
  and, where cheap, adds a fixture for each side, or writes `gap` with the reason.
- `LOSAT_STARTUP_TRACE=1` (diagnostic build) startup trace.
- `LOSAT_WASI_THREADS_DEBUG=1` threaded-WASI scheduling diagnostics.
- `LOSAT_WASI_PRESPAWN=auto|N` (Node threaded runner; `auto` when unset since
  2026-10-18, `0` boots none) boot N thread workers (`auto`: `-num_threads` - 1)
  before `_start`; with `auto`, a spawn that finds none left boots the next
  `-num_threads` workers together (`standby_refill` event); same output, extra
  `standby_*` debug events. See
  `LOSAT/tests/README.md`, "Threading contract gates".
- `LOSAT_WASI_RECYCLE=0|1` (Node threaded runner, reactors; `1` when unset since
  2026-10-18) keep the thread workers between calls; a returned thread's worker
  instantiates the module again and waits for the next thread. Same output.
- `LOSAT_WASI_GUARD_CACHE=<dir>|0` (Node threaded runner; on by default since
  2026-10-19 in `losat/wasi-guard` of the user's cache directory, `0` turns it
  off) reuse the shared-memory guard's output, checked by SHA-256 and rebuilt
  when it does not match; keeps the 8 entries used last. See
  `LOSAT/tests/README.md`, "Threading contract gates".
- `LOSAT_TBLASTX_PARALLEL_CHUNKS=1` force TBLASTX subject-chunk parallel path
  for diagnostics.
- `LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS=1` diagnostic-only sequential TBLASTX
  scan-interior chunking; this does not enable parallel scan work.

## Performance switches (on by default)

Ported from the round 3-5 performance work (2026-10, `perf/strict-set` and the packages after it).
Each switch of the strict set gives the reference result by construction (the same integer
results, or the same floating-point operations in the same order, or only a change of order,
placement or reuse), so none of them changes the output. Until main 62102f7a every switch was
off unless its variable was set. Since `perf/strict-default` (2026-10-12) the switches of the
strict set are on by default (43 at first; the block below lists the set as it is now), and
their variables mean
(`LOSAT/src/utils/xswitch.rs`):

| `LOSAT_X_<NAME>` | the switch is |
| --- | --- |
| `0` | off |
| not set, empty, or any other value | on (a value switch takes the number if it is in range, otherwise its built-in value) |

- `LOSAT_X_DEFAULT`, which turned every switch of the set off at once (the previous, reference
  paths), was removed on 2026-10-23 (switch-tidy item 3); a run ignores it. CI runs the default:
  the tests of a switch compare it with its own `=0`, and `LOSAT_X_REACH` (above) checks that the
  default reaches the fallback sides that `LOSAT/tests/switch_reach_expected.tsv` requires.
- `LOSAT_X_<NAME>=0` turns one switch off, to find which path a problem comes from. Single
  switches off are a diagnostic aid: `TBNKEEPFEW` only with four or more query batches;
  `TBNBATCHRULE` only where query batches are searched side by side (2 or more threads, four or
  more batches); `TBNOWNPOOLS` only with `TBNBATCHRULE`, natively; `TBNREDOLPT` only where the
  composition-based redo runs on the pool (2 or more threads, 2 or more matches); `REDOCONTAINIDX` only where
  one query's list of kept alignments in a composition-based redo reaches 32; `TBNXLATESIMD` only
  on a packed subject (`TBNSTORE`: `-subject`, the plans of `run_local`'s letters; `-db`) and with AVX2; `BPSUBJRUN` acts only in the parallel `-subject` loop of blastp (2 or
  more threads, several subjects); `BPDIAG` batches only with AVX2 (else the per-hit loop); `BPDEFER` acts only
  in `BPDIAG`'s AVX2 batch; `BPMASKED` only with `BPDEFER`'s deferred calls, on tables of at least `BPMASKEDMIN` hits per used cell; `BPGAP` runs the
  vector alignment only with AVX2, BMI1 and LZCNT (else the reference function); `BPGAP16` acts only
  in `BPGAP`'s vector alignment, where its costs fit 16 bits; `BPQBLOCK` takes the block path only for one output of outfmt 0, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15 or 16, or a
  `run_local` call with several outputs, the hit records or an observer, with more than one block (see its entry);
  `BNSTREAM` only in blastn's `-subject` search at 2 or more threads;
  `TBXNEXT` only in TBLASTX's parallel subject loop (2 or more threads, several subjects; native
  builds and the threaded Wasm module); `TBXBOARD` in native builds with 2 or more threads (every batch: one
  subject, or several in `TBXNEXT`'s loop); `BPCLUT` acts on `blastp -task blastp-fast` only, and `BPLKPFLAT` only where the search's pool
  builds the flat table of `BPCLUT` (`BPCLUT` on, 2 or more threads), so with `BPCLUT=0` it changes
  nothing.
- White space around a value is ignored (`set LOSAT_X_BPSUBJRUN=0 && ...` in cmd.exe stores `0 `).
  Only `0` means off; `00`, `false` or `off` are "any other value".
- The variables are read once per process (once per Wasm instance). The Web
  reactors get an empty WASI environment, so they always run the default.
- The `...SHADOW` / `...VERIFY` switches run the new and the reference computation side by side
  and stop with an assertion on the first difference. They, the counters, the timing aids and the
  tuning variables are not part of the set: they are off unless set to any value, and a shadow
  switch turns its switch on whatever its variable says. A
  shadow compares the kernel the default selects.
- The names are provisional (`x_` prefixes, `EXPERIMENT` comments). A switch graduates once it has
  been on `dev` for 7 days and in two full-set platform-checks runs (push to `main`, weekly or
  dispatch) with equal outputs (on every platform, at 1 and 4 threads) and no failed test (2026-10-22, replacing "a few releases"): its
  variable and its `...SHADOW` mode go, an always-on switch loses its previous code, and a
  conditional or platform switch keeps its rule and its fallback path and moves to `GRADUATED` in
  `xswitch.rs` (its `x_reach!` marks and its row of `switch_reach_expected.tsv` stay).
  `LOSAT_X_<NAME>=0` does not reach a graduated switch's previous path. In the descriptions below `LOSAT_X_<NAME>=1` shows the switch, not a
  setting that is needed, and "strict-set candidate" marks the round-5 switches, which are in
  the set.

The strict set (the variables that are on; `STRICT_SET` in `xswitch.rs`), written as the
environment that turned it on before the change, one variable per line in the order of
`STRICT_SET` (by name; a unit test in `xswitch.rs` compares the two). The numbers are the
built-in values. A switch that joins the set adds its line here, in order, and says since when
and on what evidence in its own entry below; it adds no paragraph after this block and no
count of the set anywhere, so that two pull requests that each add a switch do not conflict:

```bash
export LOSAT_X_BNCSR=1
export LOSAT_X_BNCSRPACK=1
export LOSAT_X_BNDBCUT=1
export LOSAT_X_BNDBSTREAM=1
export LOSAT_X_BNDBWIDE=1
export LOSAT_X_BNDCONE=1
export LOSAT_X_BNFEWSUBJ=1
export LOSAT_X_BNHITGROUP=1
export LOSAT_X_BNHLGROW=1
export LOSAT_X_BNMAPPED=1
export LOSAT_X_BNPVONCE=1
export LOSAT_X_BNQPACK=1
export LOSAT_X_BNSEEDLISTS=1
export LOSAT_X_BNSTORE=1
export LOSAT_X_BNSTREAM=1
export LOSAT_X_BNTABPAR=1
export LOSAT_X_BNTBRANGE=1
export LOSAT_X_BPCLUT=1
export LOSAT_X_BPCLUTSCANBF=1
export LOSAT_X_BPDBCUT=1
export LOSAT_X_BPDBKEEP=1
export LOSAT_X_BPDBWIDE=1
export LOSAT_X_BPDEFER=1
export LOSAT_X_BPDIAG=1
export LOSAT_X_BPFFILL=1
export LOSAT_X_BPFLEAN=1
export LOSAT_X_BPFSEED=1
export LOSAT_X_BPFSTAGE=1
export LOSAT_X_BPGAP=1
export LOSAT_X_BPGAP16=1
export LOSAT_X_BPGAPBAND=1
export LOSAT_X_BPHITBF=1
export LOSAT_X_BPHITSFLAT=1
export LOSAT_X_BPIHPOOL=1
export LOSAT_X_BPLISTEDLEAN=1
export LOSAT_X_BPLK16=1
export LOSAT_X_BPLKDROP=1
export LOSAT_X_BPLKPAHEAD=1
export LOSAT_X_BPLKPFLAT=1
export LOSAT_X_BPLKREUSE=1
export LOSAT_X_BPMASKED=1
export LOSAT_X_BPQBLOCK=50000
export LOSAT_X_BPREDOLPT=1
export LOSAT_X_BPREDOPAR=1
export LOSAT_X_BPSEGWIN=4096
export LOSAT_X_BPSTRIPEFILL=1
export LOSAT_X_BPSUBJRUN=1
export LOSAT_X_BPTREECELLS=1
export LOSAT_X_BPUNGAPSENT=1
export LOSAT_X_BXBUILDPAR=1
export LOSAT_X_BXCAPFILL=4
export LOSAT_X_BXCHUNKCROSS=1
export LOSAT_X_BXCHUNKFLOW=1
export LOSAT_X_BXCLUTCAP=3
export LOSAT_X_BXCLUTPOOL=1
export LOSAT_X_BXCROSS=1
export LOSAT_X_BXCROSSDEFER=1
export LOSAT_X_BXCTXONCE=1
export LOSAT_X_BXDIAG=1
export LOSAT_X_BXDIAGKEEP=1
export LOSAT_X_BXENDMEM=1
export LOSAT_X_BXFLOW=1
export LOSAT_X_BXHITBF=1
export LOSAT_X_BXKAONCE=1
export LOSAT_X_BXKAPPAFIT=1
export LOSAT_X_BXPOOLN1=1
export LOSAT_X_BXSCANBYTE=1
export LOSAT_X_BXSUBJFLOW=1
export LOSAT_X_BXTASKRES=1
export LOSAT_X_BXUNGAPSENT=1
export LOSAT_X_BXUNITCAP=8
export LOSAT_X_DBIDXREAD=1
export LOSAT_X_DCCSR=1
export LOSAT_X_DCSEEDLEAN=1
export LOSAT_X_DIAGLV=1
export LOSAT_X_DIAGLV2=1
export LOSAT_X_DISCLUTPAR=1
export LOSAT_X_DPEMASK=1
export LOSAT_X_DPNEON2=1
export LOSAT_X_DPNIBBLE=1
export LOSAT_X_DPTRACEPOOL=1
export LOSAT_X_DPWINDOW=1
export LOSAT_X_FOLDLEAN=1
export LOSAT_X_GREEDYROWFUSE=1
export LOSAT_X_MB24=1
export LOSAT_X_MBCSR=1
export LOSAT_X_MBFILL=1
export LOSAT_X_MBHITPF=1
export LOSAT_X_MBLUTKEEP=1
export LOSAT_X_MBLUTONCE=1
export LOSAT_X_MBLUTPREFAULT=1
export LOSAT_X_NEWTONBATCH=1
export LOSAT_X_NEWTONDROP=1
export LOSAT_X_PAIRPARDISC=1
export LOSAT_X_PAIRPARDISCARR=1
export LOSAT_X_PAIRPARFOLD=1
export LOSAT_X_PAIRPARGATE=1
export LOSAT_X_PAIRPARHOME=1
export LOSAT_X_PAIRPARMERGE=1
export LOSAT_X_PAIRPARN1=1
export LOSAT_X_PAIRPARNA=1
export LOSAT_X_PAIRWISEPAR=1
export LOSAT_X_QGROUPPAR=1
export LOSAT_X_QPREPPAR=1
export LOSAT_X_QSTREAM=1
export LOSAT_X_REDOCONTAINIDX=1
export LOSAT_X_RPTSTREAM=1
export LOSAT_X_SEEDCOMPACT=1
export LOSAT_X_SEEDUNMAP=1
export LOSAT_X_TBNBATCHRULE=1
export LOSAT_X_TBNDIAGPACK=1
export LOSAT_X_TBNFLOWHELP=1
export LOSAT_X_TBNKEEPFEW=1
export LOSAT_X_TBNLINKCTX=1
export LOSAT_X_TBNLINKSETS=1
export LOSAT_X_TBNLUTLEAN=1
export LOSAT_X_TBNOWNPOOLS=1
export LOSAT_X_TBNPAIRCAP=1
export LOSAT_X_TBNREAPIDX=1
export LOSAT_X_TBNREDOLPT=1
export LOSAT_X_TBNSCANBYTE=1
export LOSAT_X_TBNSCANRESUME=1
export LOSAT_X_TBNSTAGEDROP=1
export LOSAT_X_TBNSTORE=1
export LOSAT_X_TBNUNGAP=1
export LOSAT_X_TBNUNITXLATE=1
export LOSAT_X_TBNXLATESIMD=1
export LOSAT_X_TBXBAND=1
export LOSAT_X_TBXBOARD=1
export LOSAT_X_TBXEXTBATCH=256
export LOSAT_X_TBXFEWGROUPS=1
export LOSAT_X_TBXHITBF=1
export LOSAT_X_TBXKMERGE=1
export LOSAT_X_TBXLUT=1
export LOSAT_X_TBXNEXT=1
export LOSAT_X_TBXPIPE=1
export LOSAT_X_TBXPREFETCH=1536
export LOSAT_X_TBXQFRAMES=1
export LOSAT_X_TBXQSLICE=1
export LOSAT_X_TBXRESOLVE=1
export LOSAT_X_TBXRPTPAR=1
export LOSAT_X_TBXSCANBF=1
export LOSAT_X_TBXSTORE=1
export LOSAT_X_TBXUNGAPSENT=1
export LOSAT_X_TBXUNITMASK=256
```

Value switches: `LOSAT_X_TBXPREFETCH` (2 to
4,194,304, built-in 1,536: KiB of a subject's largest unit table), `LOSAT_X_TBXEXTBATCH` (2 to
4,194,304, built-in 256: Ki letters of a subject's largest unit query-group span), `LOSAT_X_TBXQSLICE`
(2 to 1,073,741,824: letters of a query slice, on every group longer than that; `1` or unset: the
built-in rule, slices of 131,072 letters where a group's span times the threads reaches 6 Mi
letters), `LOSAT_X_TBXUNITMASK` (2 to 4,194,304, built-in 256 for `1` or unset: KiB of the largest
masked unit diagonal table a unit takes for its size alone at eight threads or more),
`LOSAT_X_BPQBLOCK` (1 or
more, built-in 50,000: query letters per block; a value that is not a positive number is the
built-in value).

`LOSAT_X_SEEDBUCKET` was in the block of the round 3-4 set and was left out when the set became
the default: it numbered the lookup hits of one chunk scan in a 32-bit counter, exact only up to
2^32 hits per scan, which a query and a subject of a few megabases pass. It now numbers the hits
per flush, which is exact for any number of hits, and joined the set on 2026-10-12. It graduated
on 2026-10-22 (see its entry below): its rule, now without a variable, acts on the TBLASTX
searches that the PAIRPAR unit path leaves to the scan loop and whose
diagonal table has at least 2^21 cells (a query of 524 kb or more): a query of 0.9-1.05,
1.8-2.1, 3.6-4.2 or 7.2-8.4 Mb against a subject longer than the table's spare cells
(`x_pair_par.rs` `guard`), and on threaded Wasm also a subject of more than one chunk (more than
15 Mb) at more than one thread. Measured there on 2026-10-12 (alternating pairs with and without
the switch, the same output in every run; Linux x86_64 with a 36 MiB L3, node 26 for Wasm), as
time with the switch over time without. On a quiet machine: a 1.0 Mb query against 1.0 Mb (the
smallest table it acts on, 8 MiB) at 8 threads 0.84-1.12 natively (median 0.88 of five pairs) and
0.77-0.98 on threaded Wasm, at one thread 0.85-1.09 natively (median 0.99) and 0.79-0.96 on
serial Wasm; a 1.0 Mb query against 5.0 Mb 0.99 and 1.08 (native, 8 threads); a 2.0 Mb query
against 5.0 Mb (16 MiB table) 0.62 and 0.67. With other work on the machine the scan order slows
down more than the bucketed one: most pairs at 1.0 Mb were then between 0.40 and 0.89, single
ones up to 1.44. A 0.6 Mb query against a 15.1 Mb subject on threaded Wasm with a 4 GiB memory:
0.76 and 0.85 (loaded machine). The Wasm memory reaches the same size with and without it.

`LOSAT_X_PAIRLUT` and `LOSAT_X_HSPINPLACE` were added after the default changed, as switches
that were off unless set, and joined the set on 2026-10-12. Both lower the memory of TBLASTX:
the unit path (PAIRPAR) then builds only the tables it scans with, and the HSP list
of a subject is sorted and linked without a second whole copy. PAIRLUT acts only
where the unit path runs. Both graduated on 2026-10-22 (see their entries below): `LOSAT_X_HSPINPLACE`
has no variable and no `LOSAT_X_HSPINPLACESHADOW`; PAIRLUT has no variable and no
`LOSAT_X_PAIRLUTSHADOW`.

`LOSAT_X_TBNOIDCHUNK` was added after the default changed, as a switch that was off unless set,
and joined the set on 2026-10-13. It changes only which subjects one entry into the tblastn
search pool takes (2 or more threads, several subjects) and needs no other switch. Since
2026-10-14 its chunks are bounded by residues, as NCBI's database reader cuts them for its
search threads. It graduated on 2026-10-22 together with `LOSAT_X_TBNOIDNEXT`,
`LOSAT_X_TBNOIDUNITS`, `LOSAT_X_TBNSSIDE` and `LOSAT_X_TBNBATCH` (no variable; the rules stay,
see their entries below).

`LOSAT_X_SEEDRELEASE` was added after the default changed, as a switch that was off unless set,
and joined the set on 2026-10-13. It lowers the peak memory of TBLASTX further: the tables that
only the seed stage reads are freed before the HSP list of a query batch's last subject is
merged and linked. It changes no value and no order of computation, acts on every TBLASTX path
(the unit path and the scan loop) and needs no other switch. It graduated on
2026-10-22 (see its entry below).

`LOSAT_X_GREEDYFUSE` was added after the default changed, as a switch that was off unless set,
and joined the set on 2026-10-14. It ran the preliminary score-only greedy alignment of megablast
as one fused loop over distances. It changed no value and no order of computation and acted on
the blastn searches whose preliminary stage takes the score-only greedy alignment on an unpacked
subject. It graduated on 2026-10-22 (see its entry below).

`LOSAT_X_LINKSORT` was added after the default changed, as a switch that was off unless set, and
joined the set on 2026-10-14. It speeds up TBLASTX at two or more threads: the linking call sorts
its output order once where it sorted twice, and a long index sort runs as a stable parallel sort
on the search pool. It changes only how an index order is computed, not the order. It graduated on
2026-10-22 (see its entry below), with its shadow `LOSAT_X_LINKSORTSHADOW`.

`LOSAT_X_NEWTONLANES8` was added after the default changed, as a switch that was off unless set,
and joined the set on 2026-10-15. On aarch64 the Newton lanes (blastp and tblastn) solve eight
problems per round instead of four; it changes only which independent problems share a vector, not
the operations of a problem. It graduated on 2026-10-22 together with `LOSAT_X_NEWTONLANES`: no
variable; the eight lanes run on aarch64 only, the lanes stay four wide on x86-64 and Wasm.

`LOSAT_X_TBNSTATE` and `LOSAT_X_TBNXLATE` were added after the default changed, as switches that
were off unless set, and joined the set together on 2026-10-15. They cut the work tblastn does per
subject on the unit path: a thread keeps the scan state of the units from one
subject to the next instead of making and clearing it per subject and per unit, and the six frames
of a subject are read from codon tables instead of through letters. Neither changes a value; they
were timed together because `TBNSTATE` alone was 3 % slower on the fly genome at 8 threads. Both
graduated on 2026-10-22 (no variable, no `LOSAT_X_TBNXLATESHADOW`; the clearing per unit and the
letter-wise subject translation of tblastn are gone).

`LOSAT_X_TBNOIDNEXT` and `LOSAT_X_TBNOIDUNITS` were added after the default changed, as switches
that were off unless set, and joined the set together on 2026-10-15. They change how tblastn hands
the chunks of `LOSAT_X_TBNOIDCHUNK` to its threads (2 or more threads, several subjects): the
threads take the next chunk when they are free instead of waiting at the end of each entry, and
the units of a subject longer than a chunk's budget (a chromosome) run on the pool instead of on
the thread that took it. Neither changes a value or the order in which the results are entered;
they are timed together because `TBNOIDUNITS` alone was slower on the human genome at 8 and 16
threads. On main's path before the kept scan state joined the set they cost 5-11 % more CPU time on genomes;
with it the CPU time is the same. Both graduated on 2026-10-22 (no variable; the pool takes the
next chunk when a thread is free, and a subject beyond a chunk's budget runs its units on the
pool; the entries of chunks and of one subject per thread are gone).

`LOSAT_X_LINKSLICES` was added after the default changed, as a switch that was off unless set, and
joined the set on 2026-10-15. It lowers the memory of TBLASTX and, at two or more threads, its time:
the linking call links the groups of a subject's HSP list where they lie in the one sorted list,
where it cut the list into one Vec per group and copied the linked groups into a new list, so the
list no longer exists twice. It changes where the groups are linked, not a value or an order. It
graduated on 2026-10-22 (see its entry below), with its shadow `LOSAT_X_LINKSLICESSHADOW`.

`LOSAT_X_BPSUBJRUN` was added after the default changed, as a switch that was off unless set, and
joined the set on 2026-10-16. With two or more threads blastp's `-subject` loop searches the
subjects in runs of consecutive subjects, each run with one diagonal table that is cleared once and
one search scratch per thread, instead of a new scratch for every piece rayon splits the subjects
into and a clearing of the table for every subject. It changes no value; a run's later subjects
start at the diagonal offset the one-thread loop has. On the workstation (16 physical cores, medians
of alternating pairs, every output hash equal, max RSS equal) `blastp-fast` of the human proteome
against the yeast proteome took 0.89 of the time at 16 threads (6 pairs; 0.88 with
`LOSAT_X_BPQBLOCK=20000`), 0.90 and 0.92 at 8 and 0.92 at 4, exact `blastp` 0.998 at 16 (8 pairs)
and 0.995 at 8 (4 pairs). One thread and `-db` are not affected.

`LOSAT_X_TBNSTORE` was added after the default changed, as a switch that was off unless set, and
joined the set on 2026-10-16. tblastn `-subject` holds each subject packed (about a quarter of a
byte per base) instead of its letters, translates a subject's frames from the packed bases and
makes only the windows of a subject that the redo and the report read; it changes where the bases
are held, not a value. Timed on the workstation (main 20cfc555 with the switch's commits, rounds
alternating with `LOSAT_X_TBNSTORE=0`, the same reports): at 8 threads 0.94-1.00 of the time
without it (one protein against the zebrafish transcriptome 0.94 and genome 0.95, 55 proteins
against 16S 1.00, searches of 2 to 4 batches against the fly, zebrafish and human genomes and the
transcriptomes 0.97-1.00), at one thread 0.98 (zebrafish genome) and 0.99 (zebrafish
transcriptome); peak memory fell in every case, for example 2,556-2,627 to 1,261-1,263 MiB (one
protein against the zebrafish genome, 8 threads) and 1,863 to 975-979 MiB (the human
transcriptome).

`LOSAT_X_DPNEON` was added after the default changed, as a switch that was off unless set, and
joined the set on 2026-10-15. On aarch64 it let the vector X-drop kernels run, built over NEON, the
kernel bodies x86_64 and wasm32 run. It changed no value and no order of computation and acted on
aarch64 only. It graduated on 2026-10-22 (see its entry below).

`LOSAT_X_TBNBATCHRULE` was added after the default changed, as a switch that was off unless set,
and joined the set on 2026-10-16. Where query batches are searched side by side (the rule that
was `LOSAT_X_TBNBATCH`) it chooses which tblastn query batches
are searched side by side (the full waves of a small subject set) and which one after another
with the pool inside each; every batch is searched as before and the reports keep the batch
order, so no value changes. Timings are in its entry below.

The eight blastp switches `LOSAT_X_BPSEED`, `LOSAT_X_BPUNGAP`, `LOSAT_X_BPGAP`, `LOSAT_X_BPLUT`,
`LOSAT_X_BPPFREE`, `LOSAT_X_BPPCOMPACT`, `LOSAT_X_BPREDOFIT` and `LOSAT_X_BPQBLOCK` were added
after the default changed (the OrthoFinder work, 2026-10-13/14), as switches that were off unless
set, and joined the set on 2026-10-14, after a 49-job full-scale run of 7-species protein
searches with all of them on (the blocks at 50,000 letters): every output equal to NCBI's
`-db` output (or the known `-subject` / `-db` difference), and the max RSS below NCBI's on
every job. They change no value and no order of computation. `LOSAT_X_BPQBLOCK` is a value
switch whose built-in value is 50,000 query letters: a `-outfmt 6` blastp search of more than
50,000 query letters in more than one query batch now runs in query blocks.
`LOSAT_X_BPSEED`, `LOSAT_X_BPUNGAP`, `LOSAT_X_BPLUT`, `LOSAT_X_BPPFREE`, `LOSAT_X_BPPCOMPACT` and
`LOSAT_X_BPREDOFIT` graduated on 2026-10-26 with `LOSAT_X_BPCOMPOREUSE` and `LOSAT_X_BPSCANBF` (see the
entries; the tuning variable `LOSAT_X_BPSEED_PREFETCH` became the constant 16).

`LOSAT_X_BNSTREAM` was added after the default changed, as a switch that was off unless set, and
joined the set on 2026-10-16. With two or more threads blastn hands each subject's HSPs of a
subject set (`-subject`) to the hit lists in OID order from a window of at most 64 subjects per
thread while the threads search the next ones, instead of holding the HSPs of every subject of a
query batch until all are searched (one thread always hands them over at once). It changes no
value and not the order in which the hit lists take the subjects; it lowers the peak memory and,
on the workstation, the time (0.94 at 8 threads and 0.88 at 16 on 100 16S queries against 28,036
16S sequences). A database search (`-db`) has its own drivers and is not affected.

`LOSAT_X_TBXNEXT` and `LOSAT_X_TBXFRAMES` (graduated 2026-10-26: no variable) were added after the default changed, as switches that
were off unless set, and joined the set together on 2026-10-17. With two or more threads TBLASTX's
subject loop gives each thread the next subject when it is free, instead of rayon's pieces of the
subject range cut before the search (a genome's chromosome arms, listed first, fell into one piece
and were searched one after another on one thread); and a subject's frames are made in windows when
the seed stage reads them and freed after it, the frames of its letters only when HSPs reach the
re-evaluation. Neither changes a value or the order of the reports. They joined together because
the hand-out alone holds the frames of several long subjects at once. On the workstation (node 0,
medians of alternating pairs, every report equal to NCBI's `-num_threads 1` report): 3 fruit fly
transcripts against the fruit fly genome database (1,870 sequences) 11.59 -> 9.75 s at one thread
(NCBI BLAST+ 2.17.0 `-db` 10.05 s), max RSS 257 -> 107 MiB (NCBI 162); 11.48 -> 2.58 s at 8 threads
(NCBI 2.73 s), 257 -> 316 MiB (NCBI 359); 3.52 s at 4 threads (NCBI 3.50) and 2.58 s at 16 (NCBI
2.71); 10 transcripts against the zebrafish cDNA database x0.98 at one thread and x0.94 at 8. The
paths without them (`LOSAT_X_DEFAULT=0`) ran at x0.991-1.000 of main's time (V-PERF).

`LOSAT_X_BPDIAG` (the two-hit diagonal test of a lookup cell eight hits at a time) was added after
the eight, off unless set, and joined the set on 2026-10-17 after a second run of the 49 jobs with
all nine on: every output equal to the run with the eight (so to NCBI's `-db` output or the known
difference), the max RSS below NCBI's on every job, and the sum of the job times 0.457 of NCBI's
(0.597 with the eight; NCBI's times from an earlier run of the same jobs). It changes no value and
no order of computation; without AVX2 it keeps the per-hit loop. Its diagnostic
`LOSAT_X_BPDIAG_SIMD=0` (the scalar batch) stays outside the set.

`LOSAT_X_SEEDCOMPACT` was added after the default changed, as a switch that was off unless set, and
joined the set on 2026-10-17. A unit of TBLASTX's parallel seed stage keeps its saved ungapped hits
as 28-byte records, the shape NCBI BLAST+ keeps them in, copied out at their exact length, and the
merge makes each chunk's list at its exact length and moves the first one into the frame's list
instead of copying it. It changes no value and no order. On the workstation (8 rounds, every report
equal to NCBI BLAST+ 2.17.0) max RSS went from 109 to 87 MiB at one thread and from 157 to 148 MiB at
8 threads on the 0.66 Mb by 0.61 Mb pair (AP027131 x AP027133, genetic code 4), and from 130 to 121
MiB and from 215 to 209 MiB on a genome against itself (NZ_CP006932); the whole search took
x0.97-0.99 of the time.

`LOSAT_X_BPCLUT` and `LOSAT_X_BPLKPFLAT` were added after the default changed, as switches that were
off unless set, and joined the set together on 2026-10-17 (Owner decision 2026-10-16, "if there is a
merit, put it in"). `blastp -task blastp-fast` builds each query block's compressed lookup table flat
(count, layout, fill) and, with two or more threads, on the search's pool; `LOSAT_X_BPLKPFLAT`
schedules that pool build with one synchronization per pass, reuses the previous block's table
vectors and makes the contexts' Karlin blocks tasks of the build. `LOSAT_X_BPLKPFLAT` acts only where
the pool builds the flat table of `LOSAT_X_BPCLUT`, so both joined in one change. Neither changes a
value: the table and the report are the reference's. Timings (workstation, 2x Xeon 6338; every output hash equal
to the expected one) are in the entries below: with both on, `blastp-fast` of the human proteome
against the yeast proteome with `-subject` took 0.784 of the time at 16 threads with
`LOSAT_X_BPSUBJRUN` also on (0.741 with `LOSAT_X_BPQBLOCK=20000`), and `LOSAT_X_BPLKPFLAT` alone,
with `LOSAT_X_BPCLUT` on in both arms, 0.893 (0.858) at 16 threads and 0.938 (0.928) at 8; exact
`blastp` builds no compressed table and is not affected.

`LOSAT_X_BPDEFER` (the candidates of a lookup cell's AVX2 batches called after the cell's last batch
of eight) was added after `LOSAT_X_BPDIAG`, off unless set, and joined the set on 2026-10-17 after a
third run of the 49 jobs with all ten on: every output equal to the run with the nine, the max RSS
below NCBI's on every job (at most 0.923 of it), and the sum of the job times 0.407 of NCBI's (0.457
with the nine; NCBI's times from the same earlier run). On the same jobs, both arms in one run, its
summed CPU cycles were 0.938 of the run without it. It changes no value and keeps the order of the
calls; it acts only where `LOSAT_X_BPDIAG`'s AVX2 batch runs.

`LOSAT_X_TBNOWNPOOLS` was added after the default changed, as a switch that was off unless set,
and joined the set on 2026-10-17. With `LOSAT_X_TBNBATCHRULE` it searches
the tblastn query batches the rule leaves one after another at the same time, each on a search
pool of its own with its share of the threads; every batch is the same search and the reports
keep the batch order, so no value changes. Timings are in its entry below.

What runs per target:

- x86_64: the vector X-drop DP (SSE4.1 or AVX 128-bit; the 256-bit kernels with AVX2;
  the `DPW128` kernels), the AVX2 builds of the greedy rows and of the Newton kernels where the CPU has AVX2,
  the delayed megablast lookup fill (`MBDELAY`) with the prefetch instruction.
- wasm32 (built with simd128): the simd128 X-drop DP (the `DPW128` kernels replace the first kernel),
  the batched megablast lookup fill (`MBBATCH`) instead of `MBDELAY`, the Dekker product in
  the ErfC emulation (`DEKKER`, graduated); everything else is the same Rust as on x86_64
  without AVX2.
- aarch64: the NEON X-drop DP (the `DPW128` kernels where the 128-bit kernel would run; no
  256-bit AVX2 kernel), the batched lookup fill (`MBBATCH`) instead of `MBDELAY`; the Dekker product is
  compiled for wasm32 only; the Newton lanes take eight problems per round;
  everything else as on wasm32.
- Any other target: no vector X-drop DP is compiled (the scalar port runs), and the Newton lanes stay four wide; otherwise as
  on aarch64.
- One thread (`-num_threads 1`, a build without the `parallel` feature, wasm32 without
  `wasm-threads`): the scheduling switches (`AHEAD`, `GREEDYSPEC`, `PAIRPAR` for
  blastn/megablast, `BXPAR`,
  `BXPOOL`, `BXBATCH`, `BXCHUNK`, `BNSTREAM`, and the rules of the graduated `SEGSHARE`,
  `TBNBATCH`, `TBNOIDCHUNK`, `TBNOIDNEXT`, `TBNOIDUNITS`) leave the reference code running. The unit paths have no thread-count gate: the tblastn units
  (only its parallel redo needs a pool) and the TBLASTX `PAIRPAR`
  units run in order on the calling thread.

Same values, faster computation:

- DPFAST (graduated 2026-10-22; no variable, no `LOSAT_X_DPSHADOW`; always on) X-drop gapped DP
  (score-only and traceback) in 128-bit vectors (16 x i8 or 8 x i16 lanes; AVX2/SSE on x86_64,
  simd128 on Wasm, NEON on aarch64), blastn task blastn and the protein programs
  (`utils/xdrop_simd.rs`). The scalar ports of `Blast_SemiGappedAlign`, `s_BlastAlignPackedNucl` and
  `ALIGN_EX` stay for the calls that do not fit 8/16 bits and for a target without a vector kernel
  (`KERNEL_COMPILED`); the side-by-side comparison at the call sites went with `LOSAT_X_DPSHADOW`.
  The random unit test `simd_kernels_match_scalar_loops_on_random_problems` still compares the
  kernels with the scalar loops. Timing aids (separate variables, they stay): `LOSAT_X_DPNO8=1` (no
  8-bit kernel), `LOSAT_X_DPNOAVX=1` (no AVX2).
- DPAVX2 (graduated 2026-10-22; no variable; class platform, the name stays in `GRADUATED` for the
  reach table) the same X-drop DP in 256-bit AVX2 vectors (32 x i8 or 16 x i16 lanes, x86_64 with
  AVX2 only; the rule is `x_xdrop_avx2::available()`): the same values, lane width and band steps
  as the 128-bit kernels, so the same result and edit script; a call the 256-bit kernel does not
  take runs on the 128-bit kernel (`utils/x_xdrop_avx2.rs`). `LOSAT_X_DPNOAVX=1` turns it off.
- `LOSAT_X_DPEMASK=1` (with the 256-bit AVX2 kernel, speed; in the strict set from the pull request that
  added it, 2026-10-19; `=0` returns to the previous kernel) the 8-bit
  256-bit X-drop kernel (`align8`: blastp's and the translated programs' score-only DP, blastn's
  DP) stores the E (`score_gap_col`) of a dropped cell as minus infinity and no longer clears E
  right of the last survivor at the end of each row; it also keeps the row buffers' high-water
  mark as a max instead of a branch. A dropped cell's E is below the X-drop threshold in both
  forms, so no surviving cell's score, drop decision or script bit the traceback walk reads changes
  (the argument is at the head of `utils/x_xdrop_avx2.rs`); the fuzz test against the scalar NCBI loops and the capture replay
  run both instances. Without the clearing stores the next row's E loads are forwarded instead of
  waiting for those stores to be written (line R, R1). On the workstation (one core, captured
  problems, cycles per row): blastp score-only x0.88-0.93 (OrthoFinder Blast0_6, WSSV x
  PajaWSV), blastn x0.89-0.91 (AP027202 x LC738875), blastn traceback x0.97; store-forward blocks
  0.5-1.0 -> 0.03-0.12 per row, 7-9 fewer instructions per row. The 8-bit kernel also takes the
  traceback calls whose scores fit 8 bits (blastn; protein traceback with an unscaled matrix, for
  example `-comp_based_stats 0`). The 16-bit kernel (the traceback of the composition-based
  redo, whose matrix is scaled) is unchanged: the same change there measured x0.98-1.01.
  The rows of the vector kernels (`XdropScratch`: H and E on 16- and 8-bit lanes and the column letters `sb`; the
  traceback rows stay a `Vec`) start on 64-byte boundaries (`Rows` in `utils/xdrop_simd.rs`, no switch: placement
  only, the same values and lengths; perf/dc-fold-lean). The 256-bit kernels store 32 bytes at `PAD + first` plus
  multiples of 32 and load them in the next row, so whether those stores split a cache line now depends on the
  band start only. From glibc's 16-byte alignment it also depended on the heap's history: one more allocation
  elsewhere (that branch's) made blastn of 31,168 dmel cDNAs against 3 Mb at 16 threads x1.023 slower, with the
  same instructions and blocked store forwards x3.1, all in this kernel; with the rows aligned x1.002 against dev
  [0.999-1.004], dev with the same rows x1.003.
- `LOSAT_X_DPWINDOW=1` (every program's vector X-drop DP, memory; in the strict set from the pull request
  that added it, 2026-10-23; a value switch: `=<n>`, 128 or more, is the band start at which the rows move,
  built-in 1024; `=0` or `LOSAT_X_DEFAULT=0` returns to rows that hold every column) the vector kernels
  (`utils/xdrop_simd.rs`, `x_xdrop_avx2.rs`, `x_xdrop_w128.rs`) kept their DP rows (`H` and `E` on 8- and
  16-bit lanes and the column letters, 10 bytes per column) indexed by column, as NCBI's `score_array` is,
  so a scratch grew, zero-filled, to the last column its longest band reached: 524,288 lanes (5 MiB) on
  Sakai x EDL933, and a search grows one scratch per worker and stage (blastn `-subject` 16 threads: 92).
  Now a kernel whose band start reaches the bound moves the lanes from `first - 64`, rounded down to a
  multiple of 64, to lane 0 and lowers every column it keeps by the same amount (`XdropScratch::x_shift`:
  the band, the buffers' write ends, the filled column letters, the fence and bad-letter columns); the best
  cell and the traceback rows keep absolute columns. The lanes keep their values and the freed lanes are
  zero as lanes never written are, so every score, drop decision, script byte and cell count is the same.
  The fuzz tests against the scalar NCBI loops run each kernel with the bound at 128, 320 and off, and
  every family and lane width moves rows there (20,000 cases); `band_window_keeps_every_result` compares
  long problems on and off. All eight kernels move their rows: the six kernels of these three
  families and the two NEON2 kernels (`x_xdrop_neon2.rs`, `LOSAT_X_DPNEON2`, on by default on
  aarch64; their fuzz and window tests run in the platform checks). Memory on the workstation
  (dev cd8e0f89 / this switch / `=0`, max RSS MiB, `v4-B/results/b16c/mem1/`, same output): Sakai x
  EDL933 `-subject` blastn 16 threads 592-597 / 526-550 / 575-622, dc-megablast 16 threads 689-736 /
  614-630 / 651-684, 4 threads 289-292 / 268-272 / 296-303 (blastn 4 threads and 1 thread, the D.
  melanogaster `-db` cells and blastp unchanged within 15 MiB); NCBI `-subject` (one thread) 164-177. The
  rest of the genome-pair excess is the traceback rows of the long alignments, one per worker. Time
  (workstation node 1, against dev cd8e0f89, per-pair median wall, `v4-B/results/b16c/summary_ab1.txt`):
  Sakai x EDL933 blastn 16 threads x0.954 (8 pairs), 4 threads x0.986 (6), dc-megablast 16 threads x0.886
  (8; user x0.838), 4 threads x0.985 (6); D. melanogaster cDNAs x genome `-db` dc-megablast 16 threads
  x1.011 (4 pairs, 0.993-1.015; user x1.000, instructions x0.9995); blastp sp4 `-subject` x1.000 (6),
  yeast `-db` x1.001 (4); user-space instructions x0.970-0.9995 in every cell.
- `LOSAT_X_DPNIBBLE=1` (every program's vector X-drop traceback, memory; in the strict set from the
  pull request that added it, 2026-10-24; `=0` returns to one byte per cell) the
  traceback rows of the vector kernels (`XdropScratch::trace`, NCBI's `edit_script`) held one byte of
  `SCRIPT_*` bits per DP cell of the band, so a long gapped alignment kept one byte per cell: up to 36 MiB
  for one Sakai x EDL933 traceback, and each worker keeps its largest. Now the kernels store a code of
  four bits per cell (bits 0-1 the operation: 0 `SCRIPT_GAP_IN_A`, 1 `SCRIPT_SUB`, 2 `SCRIPT_GAP_IN_B`;
  bit 2 extend A, bit 3 extend B), two adjacent cells per byte, packed in registers before a store of half
  the size (x86 `pmaddubsw` + `packuswb`, wasm32 a shift, an or and a byte shuffle, aarch64 `usra` +
  `xtn` / `uzp1`), and the walk (`kernel::finish`) follows the codes with the same state machine and
  pushes each operation's NCBI value. A code comes from the same comparisons as the byte, so every script,
  edit operation, score and offset is the same. All eight kernels (256-bit AVX2, 128-bit DPW128 and base,
  NEON2); the fuzz tests against the scalar NCBI loops run every family in both layouts (20,000 cases),
  and `nibble_trace_keeps_every_result` compares long problems in both. The rule (conditional): a
  search's scratch takes the codes when it is made on a thread of a Rayon pool, that is in a search on
  two threads or more; a one-thread search
  (no pool) keeps one byte per cell, because there the packing cost time and the saving is one trace.
  Workstation, node 1, against dev a806eaa5 (per-pair median wall, max RSS MiB, `v4-B/results/b16e/`,
  same output): Sakai x EDL933 `-subject` dc-megablast 16 threads x0.962 (8 pairs), 613 -> 414; blastn 16
  threads x0.980 (8), 529 -> 373; dc-megablast 4 threads x0.987 (6), 277 -> 208; blastn 4 threads x1.002
  (6, 0.993-1.004), 244 -> 197; with the codes at one thread both tasks took x1.009 (cycles x1.010-1.012)
  for 140 -> 123 and 137 -> 121 MiB, hence the rule. D. melanogaster cDNAs x genome `-db` dc-megablast 16
  threads x1.014 (6 pairs, 0.985-1.045; user x0.996, cycles x0.994: steps of `time %e`), blastp sp4
  `-subject` 16 threads x1.002, 1,000 proteins 1 thread x1.002, yeast `-db` x0.997, blastx 50 kb x dmel
  proteome `-db` x1.000. NCBI `-subject` (one thread) 164-177 MiB. Against dev aa8281a7 (the merged form:
  the kernels' script codes in `const` blocks keep a one-byte kernel's register allocation, one thread at
  dev's instructions, x0.998 / x0.999): dc-megablast 16 threads x0.984 (8 pairs), 633 -> 476; blastn 16
  threads x0.986 (8), 535 -> 379; 4 threads x0.993 / x1.000 (6), 217 -> 189 / 257 -> 197.
- `LOSAT_X_DPTRACEPOOL=1` (every program's vector X-drop traceback, memory; in the strict set from the pull
  request that added it, 2026-10-26; `=0` returns to one trace per scratch) each scratch of a search
  (`XdropScratch`, one per worker and look-ahead helper) kept the largest traceback rows it ever held
  (`trace` and `rows`, NCBI's `edit_script` and `edit_start_offset`), so on 16 threads every worker ended
  the traceback stage of a genome pair holding a long trace (Sakai x EDL933: 130 dc-megablast calls of 1 MiB
  or more, spread over the whole stage). Now a trace that passes 1 MiB while it grows takes the largest idle
  pair of a pool the scratches share (the rows so far are copied over; the thread's own small buffers wait
  in a thread-local) and gives it back after the traceback call; the pool keeps at most two idle pairs and
  frees the smallest (two left less memory than four at the same time), and the last scratch to go (the end
  of a search's traceback stage) frees them, as the scratches' own traces went with them before. The kernels
  write and the walk reads the same rows, so every script, edit operation, score and offset is the same; the
  1 MiB test runs at the kernels' row reserve, which takes only the two buffers (`x_trace_reserve`), so the
  kernels run the same instructions (one thread x1.0000 against dev). The rule (conditional): a trace that
  grew past 1 MiB by itself goes to the pool on a thread of a Rayon pool (searches on two threads or more);
  a one-thread search keeps its one trace, as it has nothing to share it with.
  `trace_pool_keeps_every_result` runs long problems through the pool on a two-thread pool. Workstation,
  node 1 (per-pair median wall, max RSS MiB, `v4-B/results/b16e2/`, same output) against dev 3be64b7d (12
  pairs at 16 threads, 6 at 4): dc-megablast 16 threads x0.945 (0.933-0.959), 473 -> 316; blastn 16 threads
  x0.971 (0.931-1.023), 360 -> 262; 4 threads x0.993 / x1.000, 190 -> 155 / 193 -> 182; D. melanogaster
  cDNAs x genome `-db` dc-megablast 16 threads x1.000 (8 pairs), 105 -> 105; blastp 3,000 proteins x
  Species6 `-subject` 16 threads x0.998 (6), 55 -> 54. One thread (6 rounds, one core): instructions x1.0000
  / x1.0001, cycles within the A/A. Probe 4 chose two idle pairs over four (dc-megablast 16 threads 307 vs
  327 MiB at x0.941 vs x0.934) and 1 MiB over 2 MiB. NCBI `-subject` (one thread) 164-177 MiB: what stays
  above it at 16 threads is the stage start, where 16 workers run long tracebacks at once.
- DPW128 (graduated 2026-10-22; no variable; class platform, the name stays in `GRADUATED`) the
  128-bit X-drop DP kernels with the per-row savings of the 256-bit AVX2 kernels (first and last
  survivor kept while the row runs,
  row gap as a prefix maximum with a two-operation carry, the first lane's missing diagonal as a
  lane mask, column residues copied 16 at a time, row tables of `Rows28` / `Flat` matrices built
  with vector instructions, byte shifts in the form V8 compiles to `palignr`): the same values,
  lane width, blocks and result as the 128-bit kernels; a call it does not take (8-bit lanes with
  `16 * gap_extend - (gap_open + gap_extend) > 128`) runs on the 128-bit kernel
  (`utils/x_xdrop_w128.rs`). The rule: where the 128-bit kernel would run (`x_xdrop_w128::available()`:
  SSE4.1 on x86_64 after the 256-bit kernels, which keep priority on an AVX2 machine; always on
  wasm32 simd128 and aarch64). On wasm32 it replaces the simd128 128-bit kernel (the only DP kernel there).
- DPNEON (graduated 2026-10-22; no variable; class platform, the name stays in `GRADUATED`; in the
  strict set from 2026-10-15) the vector X-drop kernels on aarch64: the 128-bit kernels and the DPW128
  kernels built over a NEON `arch` / `ext` layer (`mod arch` in `utils/xdrop_simd.rs`, `mod ext` in
  `utils/x_xdrop_w128.rs`; the kernel bodies are the ones x86_64 and wasm32 run). They always run
  on aarch64; there is no switch back to the scalar port, which stays for the calls the kernels
  decline. Platform checks no longer have the mode `noneon`; they hash the default (NEON) and, on
  the Linux and macOS arm64 runners, `nodpneon2` (`LOSAT_X_DPNEON2=0`, the 128-bit NEON kernels).
- `LOSAT_X_DPNEON2=1` (aarch64 only, speed; in the strict set since 2026-10-24: on unless
  `LOSAT_X_DPNEON2=0`, which runs the DPW128 / 128-bit NEON kernels)
  the X-drop DP kernels of the 256-bit AVX2 kernels on NEON, each
  256-bit vector in two registers (32 x i8 or 16 x i16 lanes per block, the same band steps, the
  row gap in the form of the DPW128 kernels): the same values and result as the 128-bit kernels; a
  call they do not take (8-bit lanes with `32 * gap_extend - (gap_open + gap_extend) > 128`, 16-bit
  lanes with `16 * gap_extend - (gap_open + gap_extend) > 32768`) runs on the DPW128 /
  128-bit kernels (`utils/x_xdrop_neon2.rs`). It takes the place the 256-bit AVX2 kernels have on
  x86_64 (first in the dispatch). They move their DP rows down with the band like the other
  kernels (`LOSAT_X_DPWINDOW`). Evidence: #208, bench run 38590719673 on ubuntu-24.04-arm, 1 thread
  x0.867-0.974 on the cases with X-drop DP, x0.994-1.004 without, 4 threads x0.902-0.979, equal
  output hashes. Platform checks hash `=0` on the Linux and macOS arm64 runners (mode `nodpneon2`);
  the mode `dpneon2` was equal to the default on both arm64 runners in every run since 2026-10-15.
- GREEDYFAST (graduated 2026-10-22; no variable, no `LOSAT_X_GREEDYSHADOW`, no `X_GREEDY_TEST_MODE`
  hook; always on) megablast greedy alignment, one distance in three passes of 256 cells (two on x86_64 CPUs with
  AVX2 since `LOSAT_X_GREEDYROWFUSE`), match extension 8 bases at a time, on an unpacked subject (`x_greedy_row` in
  `blastn/alignment/greedy.rs`).
  `x_greedy_row` runs every traceback call and the score-only calls the fused loop (GREEDYFUSE) declines;
  the ported distance loop stays only for lengths that do not fit the slices (no production call passes
  them). The greedy functions take an unpacked subject only; the packed-subject scans (`rem != 4`) were
  deleted at graduation.
- GREEDYPRELIM (graduated 2026-10-22; no variable, no `LOSAT_X_GREEDYPRELIMSHADOW`, no counters;
  always on) the preliminary (score-only) megablast greedy alignment reads a one-byte-per-base copy
  of the packed ncbi2na subject (the NCBI-resolved codes) with `compressed_subject = false` instead
  of the packed bytes: the same comparisons, bounds and results
  (`algorithm/blastn/alignment/x_greedy_prelim.rs`), with the fast rows. Memory: one byte per subject
  base while the subject is searched (made on the first chunk with ungapped hits, shared by its chunks
  and threads). A score-only greedy call always runs on the copy: the packed score-only call
  (`greedy_gapped_alignment_score_only`, `compressed_subject = true`), its re-export and the unit test that
  compared the two forms were deleted at graduation (D2 review; nothing the default runs reached them).
- GREEDYFUSE (graduated 2026-10-22, in the set since 2026-10-14; no variable, no
  `LOSAT_X_GREEDYFUSESHADOW`; always on) the score-only megablast greedy alignment on an unpacked
  subject (the preliminary stage): the whole
  loop over distances of `BLAST_GreedyAlign` in one function, one pass over the diagonals of a
  distance, eight bases per slide step, the same cells of the two rolling rows and of `max_score`
  written in the C order (`algorithm/blastn/alignment/x_greedy_fused.rs`, exactness argument in
  its module doc). It replaces the per-distance `x_greedy_row` there; a traceback call (an edit block)
  and a call it declines run `x_greedy_row`, and the affine alignment keeps its own path. x86_64 with
  AVX2/BMI runs a build with `tzcnt`/`lzcnt`, other CPUs and targets the plain build of the same code. The
  unit tests compare it with `x_greedy_row` (a call with an edit block) and the two builds with each other.
- COMPFAST, NEWTONEXACT (graduated 2026-10-22; no variables, no `LOSAT_X_COMPSHADOW`, no
  `LOSAT_X_NEWTONEXACTSHADOW`) composition adjustment (Newton). COMPFAST: the loops of the Newton
  iteration reordered with the same element operations in the same order (`step_bound`,
  `ScaledSymmetricProductA`, `MultiplyByA`, `MultiplyByAtranspose`, the divisions of `Dinv`, the
  column-major Cholesky factor and forward solve). NEWTONEXACT: the iteration with fixed-size storage, a
  structured right-looking column Cholesky and AVX2, same element operations
  (`core/composition_adjustment/x_newton_exact.rs`). Both always run. The line-by-line ports stay where
  the default still takes them: `optimize_target_frequencies_reference` for a problem without the
  relative-entropy constraint (`desired_re <= 0`: the unconstrained rule, or a relative entropy that could not be computed). `optimize_target_frequencies` asserts
  that the alphabet size is 20, as NCBI's only caller passes `COMPO_NUM_TRUE_AA`, so the ports for a larger
  vector or matrix (`step_bound_reference`, `scaled_symmetric_product_a_reference`, `factor_ltriang_pos_def`,
  `solve_ltriang_pos_def`) were deleted. The original-order
  `MultiplyByA` / `MultiplyByAtranspose` loops and the bitwise shadow comparisons were deleted. The unit tests
  compare the exact module and the lanes with the port bitwise. `LOSAT_X_NEWTONEXACT_LEVEL=0|2` (scalar | try
  AVX-512) and `LOSAT_X_NEWTONEXACT_PROF=1` (diagnostic build; cycle counters per phase) are separate variables and stay.
- NEWTONLANES, NEWTONLANES8 (graduated 2026-10-22; no variables, no `LOSAT_X_NEWTONLANESSHADOW`; class
  conditional / platform, the names stay in `GRADUATED` as the keys of `x_reach!`) Newton problems of the
  composition adjustment solved ahead, four problems side by side (one per AVX2 lane, each lane with the
  operations of `x_newton_exact`; `core/composition_adjustment/x_newton_lanes.rs`): blastp the first
  adjustment of the next 16 matches of a query (serial and query-parallel redo; the single-query
  match-parallel redo solves none), tblastn the adjustments of the coming windows of a match (at least 32
  problems; the window ranges, translation + SEG, computed ahead are handed to the redo loop). A call uses a
  stored result only when its input (q, row and column probabilities, relative entropy) has the same bits;
  otherwise the exact module runs. NEWTONLANES8 (aarch64 only; in the strict set from 2026-10-15 to
  2026-10-22): batches of five or more problems use eight lanes (`F8`, four 128-bit NEON operations per
  value); once the queue is empty and at most four lanes hold a problem, those problems move to four lanes
  (bit copies of the values kept between rounds) and finish there. Each problem keeps its operation
  sequence, so the results are the same bits (`x_newton_lanes.rs`, module header). On other targets the
  lanes stay four wide: on x86-64 the two 256-bit operations per value gained nothing (0 to 3 % slower on
  blastp / tblastn), and Wasm is not measured. Measured 2026-10-15 on the GitHub runners, bench5's
  blastp_wssv / tblastn_small at one thread, wall relative to four lanes: Linux arm64 (Neoverse-N2)
  0.981 / 0.959, macOS arm64 0.974 / 0.893 (macOS rounds are noisy); on Neoverse-N2 no lanes gave
  0.967 / 0.941. With `LOSAT_X_STATS=1` it prints
  `[X_STATS] NEWTONLANES prefetched= hits= misses= unused= ...` and, on aarch64,
  `[X_STATS] NEWTONLANES8 batches8= rounds8= tail_rounds=` (rounds on eight lanes and after the move; no
  `xstats` feature needed). Strict-set candidate (S-F, round 5).
- SEGFAST, SEGMEMO (graduated 2026-10-22; no variable, no `LOSAT_X_SEGSHADOW`) SEG window as a residue
  histogram updated in O(1), entropy memoised per state vector. Always run (`SegMasker::mask_sequence`)
  for a window of at most 15 letters and the 20-letter alphabet (`x_fast_applies`); other settings run
  `mask_sequence_reference`, the port of `s_SegSeq`, and a segment of more than 127 letters is trimmed by
  the port of `s_Trim`. The memo falls back to the port of `s_Entropy` for a state vector that does not
  pack (more than 15 letters). The unit tests compare the histogram path with the port on random sequences.
- SEGTRIMMEMO (graduated 2026-10-26, in the set from 2026-10-19; no variable; always on for the maskers that
  ask for it; the trim without the memo went with the variable) (speed) SEG's trim of a candidate segment (the
  port of `s_Trim`: about L^2 / 2 evaluations of `s_GetProb` for a segment of L letters up to
  maxtrim + 1 = 51, about 1,275 for longer ones) is kept per thread by the segment's letters and
  `maxtrim`, and a segment with the same letters takes the stored offsets instead of being trimmed
  again (`utils/x_seg_trim_memo.rs`; the histogram SEG path). The offsets are a function of those letters,
  `maxtrim` and the fixed 20-letter alphabet; the whole slice is compared, the hash only picks the
  slot. Only the maskers that repeat segments use it (`SegMasker::x_with_trim_memo`): TBLASTN's
  redo windows (a window around every alignment, so a low-complexity stretch of the subject is
  trimmed once per query whose window covers it: AP027131 x AP027133 118,765 trims of 12,584
  distinct segments, 89 % of the trim work repeated) and BLASTX's query SEG (the chunks of a split
  query repeat the whole query's segments: 51 % of the trim work). TBLASTX d04 would repeat 2 %, and
  the other maskers keep no table. A thread keeps at most 32,768 segments and 1 MiB of letters
  (table up to 768 KiB), then clears the table; segments over 127 letters are not kept (line C,
  v4-C DECISIONS 31-34).
  Workstation, one core or eight, 8 rounds (2026-10-19): TBLASTN AP027131 x AP027133 n1 0.924 of
  the time (pairs 0.921-0.926), n8 0.971; BLASTX AP027131 x AP027133 n1 0.927, n8 0.985 (memo on
  every masker, e91c8864); with the two maskers of this version (cc7306b2, load 12-45) n1 0.932 /
  0.935. WSSV x PajaWSV, OrthoFinder q0 x Species6 n8 and TBLASTX d04 do not use it. Max RSS +1 to
  +4 MiB (BLASTX n8 68 MiB; NCBI n1 88 MiB). At eight threads the repeats fall on other threads'
  tables, hence the smaller gain.
- SUMPMEMO (graduated 2026-10-27, in the set from 2026-10-23; no variable; always on) (BLASTX, TBLASTN and
  TBLASTX, the sum P-value of linked HSPs; speed) `blast_sum_p` (the port of
  `s_BlastSumP`, `stats/sum_statistics.rs`) keeps per thread the arguments (`r` and the bits of `s`) and values of its
  last two calls of `blast_sum_p_calc` (the port of `s_BlastSumPCalc`: for r > 4, a Romberg integral of a Romberg
  integral, hundreds of `exp` and `log` per call) and returns the stored value for the same arguments. The value is a
  function of `r` and `s` alone. TBLASTX's linking asks for the sums of the best chain of each ordering method in every
  round (`link_hsps.c:901-953`), and the chain that was not chosen is asked for again in the next round: AP027131 x
  AP027133 (genetic code 4) at one thread calls `blast_sum_p_calc` 59,452 times, 56,780 of them with the arguments of
  one of the two calls before (the same chain up to 296 times in a row), MelaMJNV x PemoMJNVA 16,121 and 14,796. BLASTX
  and TBLASTN of AP027131 x AP027133 do not reach `blast_sum_p_calc`. Workstation, one core, the same binary with
  `=0` and without, 8 rounds in alternating order (2026-10-23, outputs identical in every run): TBLASTX AP027131 x
  AP027133 x0.942 of the time (cycles x0.946, instructions x0.943), MelaMJNV x PemoMJNVA x0.950 (cycles x0.945);
  eight threads (eight cores, siblings idle) x0.923 and x0.965 (cycles x0.944 and x0.946). LOSAT Web (Node, serial
  reactor, 8 rounds in one process against dev 32a63700): the whole search x0.945 and x0.946 of the cycles, the sum
  P-value stage x0.09 and x0.15.
- DISPTABLE (graduated 2026-10-27, in the set from the pull request that added it, 2026-10-24; no variable; always
  on; the scan of each codon stays as the reference of the unit test) (TBLASTX and TBLASTN, the displayed residues;
  speed) `display_residue` (`tblastx/report.rs`) reads a
  codon's displayed residue from a table of the 4,096 codons of 4-bit codes, made once per genetic code at its first
  use (`display_table`, each entry `display_codon` of its codon; `GeneticCode::display_slot`), and the bases' codes
  from two 256-byte tables (plus strand and complement). NCBI's display does the same (`CTrans_table`: `sm_BaseToIdx`
  and `m_AminoAcid`, `Genetic_code_table.cpp`). Before, each base's code was found among the IUPAC letters and each
  codon translated by the nested loops of `display_codon`. The residues are the same (unit test
  `display_table_gives_the_residues_of_the_scan` over every byte, codon and code). TBLASTX translates the displayed
  rows again for the identities of every HSP, the outfmt 0 rows and LOSAT Web's hit records; TBLASTN for its subject
  rows. The table of a code is filled once per process (4,096 calls of `display_codon`). Workstation, the same binary
  with `=0` and without, 8 rounds in alternating order (2026-10-24, outputs identical in every run): TBLASTX
  AP027131 x AP027133 (genetic code 4) with outfmt 0 at one thread x0.979 of the time (cycles x0.980, instructions
  x0.987), MelaMJNV x PemoMJNVA x0.957 (cycles x0.960); with outfmt 6 (the identities only) x0.991 and x0.982;
  eight threads (eight cores, siblings idle) with outfmt 0 x0.994 and x0.971 (cycles x0.982 and x0.965). LOSAT Web
  (Node, serial reactor, 8 rounds in one process, the same reactor with `=0`): the whole search x0.984 and x0.967
  of the cycles, the report stage x0.63 and x0.46. TBLASTN was not timed.
- DUSTFAST, DUSTRING (graduated 2026-10-22; no variables, no `LOSAT_X_DUSTSHADOW`) DUST perfect-interval
  list built in one merge, window kept in a ring buffer. Always run: `TripletWindow::find_perfect` is the
  merge (the literal `std::list` port was deleted), and `DustMasker::mask_subsequence` uses the ring buffer
  for a window of at most 64; `mask_subsequence_reference`, the port of `CSymDustMasker::operator()`, runs
  for a larger window. The unit test compares the ring buffer with that port on random sequences.
- KARLINFAST (graduated 2026-10-22; no variable) Karlin K convolution in eight lanes (same per-lane
  order). Always run where `first == 0` and eight positions remain; the scalar loop runs for `first != 0`
  and for the remaining positions. The unit test that compared the lanes with the scalar-only run went
  with the variable.
- ERFMEMO (graduated 2026-10-22; no variable) Spouge ErfC memoised per thread by the argument's bits;
  always on, `erf_impl_ncbi` computes the missing values.
- ERFABSORB (graduated 2026-10-26, in the set from 2026-10-19; no variable; always on) (speed) Spouge ErfC returns 2.0 for an
  argument below -6 without evaluating ErfImpl. The reference computes `2.0 - ErfImpl(-z, TRUE)`
  there, and ErfImpl of 6 or more is below 2.2e-17 (a fifth of 2^-53; the bound is argued at
  `erfc_ncbi` in `stats/spouge.rs`), so the difference rounds to 2.0. The unit test
  `x_erfc_below_minus_six_is_two_and_the_switched_function_has_the_reference_bits` checks the
  reference on a grid, on random arguments and on the doubles next to -6.
  ErfImpl's long-double `expl` (double-double, 40 series terms and nine squarings) costs about
  5,000 cycles a call; on WSSV x PajaWSV (blastp, one thread) 13,046 of the 18,964 calls that
  reached it had an argument below -6 (line C, v4-C DECISIONS 25). Workstation, one core, 8
  rounds (2026-10-19): WSSV x PajaWSV n1 0.969 of the time with this switch alone (0.943 with
  BPCOMPOREUSE as well, n8 0.958), tblastn AP027131 x AP027133 0.994, blastx
  AP027131 x AP027133 0.987, OrthoFinder q0 x Species6 n8 0.995 and n1 1.00 (both switches).
- DEKKER (graduated 2026-10-22; no variable; platform: wasm32 only) ErfC product error by Dekker's
  product instead of `fma` (`stats/spouge.rs` `x_product_error`), where the product and the
  factors are in range (2^-900 < |product| < 2^900, |a|, |b| < 2^900); other targets and values
  out of range keep `mul_add`.
- CODONFAST (graduated 2026-10-22; no variable) an unambiguous codon translated by one table read; a
  codon with an ambiguity code still runs the port of `s_CodonToAA`'s loops (`GeneticCode::get`).
- PURGEFAST (graduated 2026-10-22; no variable; always on) blastn common-endpoint purge without
  shifting the array per removal (`filtering/purge_endpoints.rs`). The NCBI-literal shifting loops
  of both passes went; the survivors and the removed entries land in the slots the C shifts would
  leave them in.
- ITREEFAST (graduated 2026-10-22; no variable, no `LOSAT_X_ITREESHADOW`, no `X_ITREE_TEST_MODE`
  hook; always on) blastn interval tree with side indexes (grid, endpoint bitmaps) that
  answer "none contained" / "no such endpoint" without walking the tree. The tree walk stays: it
  answers while tracing is on (`LOSAT_DEBUG_*`, the index may name another container) and where
  the index was dropped (an HSP the index does not describe, or a size limit).
- CTXFAST (graduated 2026-10-22; no variable, no `LOSAT_X_CTXSHADOW`; always on) blastn
  query-offset-to-context lookup by binary search over the runs of the contexts instead of a
  per-base table. The per-base table (`direct_map`) stays for overlapping runs.
- CTXO1 (graduated 2026-10-27, in the set from 2026-10-19; no variable; always on, and the binary search
  over the runs still answers where the layout does not allow the form) (blastn, speed) the
  context of a lookup hit's query offset in constant time (`blastn/blast_engine/x_ctx_o1.rs`, argument
  in its header): where every context has the same length `L`, `offset / (L + 1)` by a multiplication
  (NCBI's `BSearchContextInfo` shortcut), else one `u32` per block of `2^k` offsets (the context at the
  block's start and the place of the block's sentinel), with `k` the largest value at which at most one
  block in 256 holds two sentinels; those few blocks keep the binary search. The answer is the per-base
  table's at every offset. The binary search reads up to 16 entries of 24 bytes per context for each
  lookup hit (R8: 52 % of the cycles of blastn-short with 30,000 queries, 318 cycles a hit). The block
  table takes 4 bytes per `2^k` query letters (375 KB for 30,000 queries of 18 to 30 letters, 2.8 MB for
  both strands of 31,168 cDNAs of 30 to 71,382 letters); it is not built above 64 MiB, nor where the
  contexts are not laid out one after another with one sentinel, and the binary search answers.
  `LOSAT_X_CTXO1FORM=blocks` (tuning, not in the strict set) builds the block table where the lengths
  are equal. The binary search of the crowded blocks is covered by a unit test.
  Measured (workstation, Xeon 6338, booked cores of one node, alternating order, medians, every report
  identical; base = the binary before): blastn-short of 30,000 queries of 24 letters against 3 Mb of
  D. melanogaster 2L at one thread 57.15 -> 29.98 s (x0.525; NCBI `-db` 28.87 s; the block table instead
  of the division x0.571), of lengths 18 to 30 57.36 -> 33.74 s (x0.588; NCBI 48.81 s); dc-megablast of
  31,168 cDNAs against the same subject x0.822 at one thread, x0.939 at 16; blastn of the same x0.949
  / x0.960; max RSS unchanged. `=0` x1.000-1.005 of the binary before (1.6 % more instructions on
  blastn-short: the binary search moved out of line behind the check for this form).
- `LOSAT_X_BNPVONCE=1` (BLASTN, the serial scan and the strip scan of a contiguous table laid out as packed CSR
  of 4^11 cells or fewer, speed; in the strict set from the pull request that added it, 2026-10-31
  (perf/abr-pv-once); `=0` tests the presence vector first) tests each word once, in the head's presence bits
  (`x_bn_pv_once`). The head has one bit per cell; the presence vector's bit stands for four cells of 4^11, so at the
  15 % of filled cells of ABRicate's form it is set for about half of the words, which then test again. A larger head
  (megablast's 4^12 cells, 4 MiB of bits) keeps the presence vector in front. Evidence (ws, one core of node 0, eight rounds, against the
  previous `dev` in the same run): ABRicate's ncbi DB x its assembly x0.950, Sakai draft 200 contigs x the same DB x0.973,
  ECTyper x0.979, Easyfig ex2a x0.996; megablast instructions x1.00000; the same reports, the same memory.
- BNUNGAPCOUNT (graduated 2026-10-27, in the set from 2026-10-22; no variable; always on; the ported loops
  stay as the unit test's reference) (blastn, speed) the exact ungapped extension
  (`s_NuclUngappedExtendExact`: blastn-short and other word sizes below 11, and the exact pass after the
  approximate one reaches the reduced cutoff) with the step count of each direction worked out first:
  `min(q_off, s_off)` to the left, `min(q_avail, s_avail)` to the right, the C loop conditions on
  `(s, base)` written as those counts (`blastn/x_ungap_count.rs`, argument in its header). The same
  steps, tests and results; the ported loops also tested the query offset against both ends (never true
  within the counts) and bounds-checked the matrix cell, and their two-part loop condition on `(s, base)`
  cost 0.8 branch misses a call more than the counts (premise probe on bs24: NCBI's loop shape without the
  extra tests kept the misses, the counted loops dropped them).
  R8 (`LOSAT_X_CTXFAST` / `CTXO1` era): 167 cycles, 574 instructions and 3.10 branch misses a call
  against NCBI's 135 / 497 / 2.56.
  Measured 2026-10-22 on the workstation (node 1, one pinned core, HT sibling idle, warmed 3 s; 3 rounds,
  the arm order reversed every other round; `perf stat` cycles against the binary before, `dev`
  249c7454; outputs identical in every arm): blastn-short of 30,000 queries of 24 letters against
  dmel_dna_3M x0.916 (30.00 -> 27.50 s), of 18 to 30 letters x0.920 (33.63 -> 30.96 s): 9 % fewer
  instructions and 18 % fewer branch misses, about 36 cycles and 96 instructions less a call over R8's
  220 M calls; blastn and dc-megablast of 31,168 dmel cDNAs against the same subject x0.995 / x1.000,
  megablast Sakai x MG1655 x0.997; max RSS unchanged. `=0` x1.008-1.010 on blastn-short (8 more
  instructions a call: the mode check now inlined in front of the ported loops, which were already a
  call of their own), x0.999-1.003 on the others.
- LUTARENA (graduated 2026-10-22; no variable; always on) protein lookup-table chains carved from one
  arena per table (`tblastx/lookup/backbone.rs`, all four protein programs); the `Vec` chain went.
- MBDELAY, MBBATCH (graduated 2026-10-22; no variables; class platform, the names stay in
  `GRADUATED` for the reach table) megablast lookup construction touches (prefetches) hash cells
  ahead; same update order. The rule is the CPU: x86_64 prefetches each word's cell and applies the
  update 16 words later (`MBDELAY`); every other target collects batches of 64 words, reads their
  cells in one loop, then applies the updates (`MBBATCH`, built-in 64, `X_MB_BATCH`). Both run in
  `x_fill_mb_queued` (`blastn/lookup.rs`); the unqueued loop of `build_mb_lookup` went, as nothing
  reaches it any more (PAIRPAR's pool fill, when it takes the build, comes first).
- ENVCACHE (graduated 2026-10-22; no variable) blastn debug variables (`LOSAT_DEBUG_COORDS`,
  `LOSAT_DEBUG_COORDS_START`) read once per process instead of per extension (`utils/xenv.rs`).
- BXLAZYCTX (graduated 2026-10-22; no variable; always on) blastx seed context lookup made after the diagonal
  test, as NCBI does (`blastx/seed.rs` `word_finder`, `blastx/x_seed_scan.rs`).
- BXLEAN (graduated 2026-10-22; no variable; always on) does not build records nobody reads (candidate copies
  when the observer is passive, `x_passive`); no tree or DP space for a (chunk, subject) without initial HSPs. The tblastn search never collects the
  test-only comparison events (the serial preliminary stage that collects them runs in unit tests
  only); this has no switch.
- BXLUT (graduated 2026-10-22; no variable, no `LOSAT_X_BXLUTSHADOW`; always on) blastx lookup table filled in
  place from per-cell counts (same append order as `s_AddNeighboringWords` and `BlastAaLookupFinalize`),
  neighbour-word lists memoised per thread (`tblastx/lookup/x_lut_direct.rs`). The chain builder
  (`build_lookup_from_prepared`) stays for diagnostics mode, a backbone index wider than 16 bits, and the
  tblastx and blastp tables that still choose between the two builders; the shadow and its
  `[X_SHADOW] BXLUT` line went (the unit tests still compare the two builders field by field).
- BXSCAN (graduated 2026-10-22; no variable, no `LOSAT_X_BXSCANSHADOW`; always on) blastx seed loop reads the
  scanner's per-thread offset-pair buffer in place instead of a fresh buffer and a copy per subject
  (`blastx/x_seed_scan.rs`); not in the diagnostic stages, where the reference `word_finder` runs (BLASTX has
  no bucketed order: its seed loop is the scan-order one). The shadow and its `[X_SHADOW] BXSCAN` and
  `[X_SHADOW] BXHITBF` lines went; the unit tests still compare the direct loop with `word_finder`.
- SEEDBUCKET (graduated 2026-10-22; no variable, no `LOSAT_X_SEEDBUCKETSHADOW`; class conditional,
  the name stays in `GRADUATED` for the reach table; it was in the strict set since 2026-10-12)
  TBLASTX two-hit stage: hits buffered per diagonal range and processed range by range, HSPs
  restored to scan order (`tblastx/x_seed_bucket.rs`). The rule stays: only for diagonal tables of
  at least `LOSAT_X_SEEDBUCKET_MIN_CELLS` cells (default 2^21) and never in a traced or debugged
  run; the fallback is the scan-order loop (`scan_chunk_plain`). It acts only where the PAIRPAR unit
  path leaves a batch to the scan loop. Tuning (variables that stay): `LOSAT_X_SEEDBUCKET_MIN_CELLS`,
  `LOSAT_X_SEEDBUCKET_CELLS` (cells per bucket), `LOSAT_X_SEEDBUCKET_BUDGET` (hits per flush;
  values above 2^31 are read as 2^31, and memory, 12 bytes per hit reserved up front, limits a
  real setting long before that). A test that needs the bucketed order on a small input uses one
  thread (the cost rule sends the batch to the scan loop), `LOSAT_X_SEEDBUCKET_MIN_CELLS=0` and
  checks `[X_STATS] PAIRPAR units=0`. With `LOSAT_X_STATS=1` the line
  `[X_STATS] SEEDBUCKET bucketed_scans=... flushes=...` tells how many scans ran in the bucketed
  order (`LOSAT/tests/tblastx_seed_bucket.rs` uses it). The shadow (HSP lists and the whole
  diagonal table compared with the scan-order loop on every chunk) went with the variable. BLASTX
  keeps the scan-order loop: an opt-in
  bucketed order for BLASTX (blastx_ap 8.4 -> 8.7 s, no gain measured) and its shadow were
  removed, and SEEDBUCKET acts on TBLASTX only.

Rebuilding less:

- SEGSHARE (graduated 2026-10-22; no variable, no `LOSAT_X_SEGSHARESHADOW`; conditional: the
  query-parallel redo and the windows of `x_redo_par`, 2+ threads) blastp: the subject SEG masks
  shared by the whole search, as in the one-thread path. The one-thread loop keeps one cache for
  the search; the match-parallel redo makes one per match (`BlastpKappaSubjectRangeCache::new`).
- BPCOMPOREUSE (graduated 2026-10-26, in the set since 2026-10-19; no variable; always on) (blastp, speed) an
  entry of the subject SEG
  cache also keeps the amino-acid composition of its masked bytes, read once when the entry is
  made, and the composition-based redo and its Newton prefetch take it instead of reading the
  same bytes again for every query (`read_aa_composition` on the same bytes gives the same
  bits). A range SEG does not mask (near-identical or not biased) is read as before; a masked
  range whose window the redo skips has its composition read once anyway. The unit test
  `x_masked_subject_range_carries_the_composition_of_its_bytes` checks the stored composition against a
  fresh read. On WSSV x PajaWSV (one
  thread) the subject composition was read 9,434 times from 86 different byte strings, all of them
  cache entries (line C, v4-C DECISIONS 25). Each entry holds a pointer and, when the switch is
  on, 232 bytes on the heap. Workstation, one core, 8 rounds (2026-10-19): WSSV x PajaWSV n1 0.978
  of the time with this switch alone (0.943 with ERFABSORB as well); max RSS unchanged.
- TBNQSIDE (graduated 2026-10-22; no variable; always on) tblastn: query SEG and lookup table
  once per query set (`XQuerySide`); the per-subject build is gone.
- TBNSSIDE (graduated 2026-10-22; no variable; conditional: its rule stays) tblastn: subject
  six-frame translations once per search of four or more
  query batches (20,000 query residues per batch), up to 256 MiB of translated residues
  (`XSubjectSide::KEEP_FROM_BATCHES`, `KEEP_LIMIT`; whether a plan is kept is decided from the
  subject's length by the first batch that reaches it, `XSubjectSide::keeps`, and a plan that is
  not kept is made per use as the plans of a search without a side are: since 2026-10-28 without
  residues when its first frame has 1 MiB or more, `TBNUNITXLATE`). The batches go to
  `x_search_batches` at any thread count. A
  search of one to three batches keeps none and translates a subject again for each batch, as
  NCBI does (`blast_engine.c:627-631`): one batch has no later reader, and two or three give
  the time the kept translations would save for their memory. Timed with two and three batches
  at 4 and 8 threads, with other work on the machine: on 28,036 16S subjects (40.6 Mb) the
  search without them has 115 MiB less (max RSS 218-226 -> 103-113 MiB) and takes a median of
  1.03-1.11 of the time, on one 606 kb genome the same memory and 1.05-1.08 of the time
  (0.01-0.06 s). Four batches at 8 threads gave 1.06 on the 16S set, so the saving does not
  start at four batches: the number is where the memory is preferred to it. Side by side
  (the rule of `TBNBATCH`) needed the kept translations, so two or three batches with at least
  one batch per thread are searched one after another instead of side by side: 1.02-1.07 of the
  time on the 16S set and 1.12-1.24 on the genome (0.04-0.09 s) with one batch per thread,
  0.75-0.86 with three batches on two threads. Since `TBNKEEPFEW` (2026-10-19) a search of four
  or more batches keeps the translations only when they hold no more than the threads have in
  flight; the batches are laid out as before (the layout reads the number of batches, not the
  kept translations).
- `LOSAT_X_TBNBATCHRULE=1` tblastn, where query batches are searched side by side (speed; in the
  strict set since 2026-10-16, on by default, `=0` turns it off; `tblastn/stage_d_pipeline.rs`
  `x_batches_side_by_side`): which query batches go side by side. With `=0`
  every batch goes side by side once there is one batch per thread. By default, only the batches of the
  full waves (a multiple of the threads) go side by side and the others one after another with the pool inside each, and
  only when the subjects hold at most 24 Mb or there are fewer than two subjects per thread
  (otherwise every batch one after another). The report does not change (every batch is searched
  as before, the reports keep the batch order). With `LOSAT_WASI_THREADS_DEBUG=1` a line
  `[losat-thread-stage] program=tblastn stage=query_batches batches= threads= side_by_side=
  one_after_another= own_pools=` gives the split (printed whenever the batches are dispatched on a
  pool, with or without the rule; `own_pools` is 0 with `LOSAT_X_TBNOWNPOOLS=0`, without the rule, in
  WebAssembly, and where its conditions below are not met). Timed on the workstation
  (perf-si-20261012 session 17, node 1, one hash per
  case; time of each arm / time of main's side by side): 16S (28,036 records, 40.6 Mb) x 10
  batches at 8 threads 0.68 (13.2 -> 8.9 s), 5 batches at 4 threads 0.69 with
  `LOSAT_X_TBNBATCH=0` (the variable of the rule that graduated 2026-10-22); one 0.6 Mb genome x
  10 batches at 8 threads 0.77 with the rule (8 side
  by side, 2 after), 1.07 with `TBNBATCH=0`; three E. coli genomes (15.7 Mb) x 16 batches at 12
  threads 0.80 with the rule; three fly chromosomes (76.9 Mb) x 16 batches at 16 threads 1.29
  with `TBNBATCH=0` (the rule keeps them side by side), the same bases in 79 pieces 0.95.
  Where all batches are a multiple of the threads and the set is small (0.3-15.7 Mb, one to
  5,000 subjects) side by side took 0.26-0.96 of the time of one after another.
- `LOSAT_X_TBNOWNPOOLS=1` tblastn, with `LOSAT_X_TBNBATCHRULE` (speed;
  in the strict set since 2026-10-17, on by default, `=0` turns it off; `tblastn/stage_d_pipeline.rs`
  `x_batches_on_own_pools`): the query batches `TBNBATCHRULE` leaves one after another (all of
  them when there are fewer batches than threads, else those after the full waves) are searched
  at the same time instead, each on a search pool of its own with its share of the threads
  (`threads / k`, the first `threads % k` one more), after the shared pool has ended. Only with
  fewer than two subjects per thread, and for two to `threads - 1` batches when the subjects
  hold at most 24 Mb (a pool may then have one thread), otherwise two to `threads / 2`. Native
  builds only: in WebAssembly only the main thread can start threads, and each pool is built on
  a thread of its own, so the batches stay one after another there. The report does not change
  (every batch is the same search, the reports keep the batch order; `tests/tblastn_batch_rule.rs`);
  the `query_batches` line above gives `own_pools=`, and one `[losat-thread-pool]` line is printed
  per pool. Timed on the workstation (perf-si-20261012 session 19, node 0, one hash per case,
  outputs equal; time with the switch / main's default): one 4.7 Mb genome (E. coli MG1655) x 16
  batches at 32 threads 0.20 (`-subject`, 3.86 -> 0.79 s) and 0.19 (`-db`), x 10 batches at 32
  threads 0.24, x 5 batches at 32 / 16 / 8 threads 0.27 / 0.45 / 0.97, x 2 batches at 8 threads
  0.89; one 0.6 Mb genome x 16 batches at 32 threads 0.15, x 5 batches at 16 / 8 threads 0.45 /
  0.61; three fly chromosomes (76.9 Mb) x 16 batches at 32 threads 0.81, x 5 batches at 16 threads
  0.89. Left out by the conditions, as measured with all of them on pools of their own: three fly
  chromosomes x 5 batches at 8 threads 1.44 (two pools of one thread), 5,000 16S records (7.7 Mb) x
  5 batches at 8 / 16 threads 1.38 / 1.09, 28,036 16S records x 5 batches at 8 / 16 threads 1.67 /
  0.98, the fly transcriptome x 5 batches at 16 threads 1.16. With `=0`, 16 batches on the 4.7 Mb
  genome at 32 threads took 1.60 (`-subject`) and 2.86 (`-db`) times NCBI BLAST+ 2.17.0 `-db` at 32
  threads (1.59 s).
- LUTSPLIT (graduated 2026-10-22; no variable; always on) tblastn: build only the half (contexts
  or table) a caller uses (`x_lookup_table_for_profile`, `x_lookup_contexts_for_profile`). The
  contexts-only caller still takes the full `build_ncbi_lookup_for_profile` in diagnostics mode.
- IDEALMEMO (graduated 2026-10-22; no variable) BLOSUM62 ideal Karlin block computed once per process
  (`x_compute_blosum62_ideal_karlin_params` is the memo's source).
- BXSEGMEMO (graduated 2026-10-22; no variable, no `LOSAT_X_BXSEGMEMOSHADOW`; always on) blastx redo stage: the
  subject SEG once per subject instead of once per window (per-thread memo of the last subject, byte-equal
  input; `blastx/x_seg_memo.rs`; the masker runs on a miss). The shadow and its `[X_SHADOW] BXSEGMEMO` line went.
- HSPINPLACE (graduated 2026-10-22, in the set since 2026-10-12; no variable, no `LOSAT_X_HSPINPLACESHADOW`;
  always on) (TBLASTX, memory) the HSP list of a subject is sorted, linked as sub-slices
  and replayed into NCBI's output order without a second whole copy alive (`tblastx/x_hsp_inplace.rs`,
  exactness argument in its header). The index sorts of `sort_ungapped_hits_by_index_replay` (score order,
  link order) move the records inside the list, cycle by cycle, after the index sort; the linking call links
  its groups as sub-slices of the list (LINKSLICES), with the same per-group functions (serially or on the
  pool), and applies the output replay as a permutation in place. A replay that leaves an HSP out or repeats
  one returns the reference's cloned list (`LOSAT_X_STATS=1` counts both:
  `[X_STATS] HSPINPLACE sorts= linkings= hits= groups= replay_inplace= replay_fallbacks= gate_reference=`).
  The per-group scratch of the literal port is not reserved unless that port can run (`LOSAT_LINKING_LEGACY`,
  a linking trace, `LOSAT_LINK_FAST_SHADOW`). A linking call with `LOSAT_DEBUG_CUTOFFS` or
  `LOSAT_DIAGNOSTICS` runs the reference linking body in `sum_stats_linking/linking.rs` (`gate_reference`;
  its cutoff debug output lives there); its first sort is the in-place sort. The reference sort copy of
  `sort_ungapped_hits_by_index_replay`, the shadow and its `[X_SHADOW] HSPINPLACE` line went.
- BPLUT (graduated 2026-10-26, in the set since 2026-10-14; no variable; always on) (BLASTP, memory) the lookup
  table of the default (non-compressed) BLASTP lookup is
  filled in place from per-cell counts by BLASTX's (BXLUT's) builder (same append order as
  `s_AddNeighboringWords` and `BlastAaLookupFinalize`; `tblastx/lookup/x_lut_direct.rs`, exactness
  argument in its header), so the per-cell chains of the reference builder, about 1.5 times the
  finished table, are never held. In diagnostics mode the reference builder runs (and a backbone index
  wider than 16 bits falls back to it). The unit test `bplut_table_and_contexts_equal_build_ncbi_lookup`
  still compares the table and the contexts with `build_ncbi_lookup` field by field.
- BPREDOFIT (graduated 2026-10-26, in the set since 2026-10-14; no variable; always on) (BLASTP, memory) the
  composition-based redo builds each redone HSP list with room for
  exactly its alignments, and a list from which HSPs were removed afterwards gives back their room
  (`blastp/x_redo_fit.rs`): the list's vector used to grow by `push`, whose first allocation has room for
  four HSPs of 168 bytes, while most lists hold one. The lists hold the same HSPs in the same order, so there
  was no shadow. With `LOSAT_TIMING=1` each search of a query set prints
  `[TIMING] x_bpredofit: lists=<n> hsps=<m>`. Measured once (human proteome block 1, 566 proteins, 50,404
  letters, x human proteome, one thread, with `LOSAT_X_BPQBLOCK=50000` and the other memory switches; max RSS):
  121,304 -> 82,180 kB (NCBI -db 151,476 kB in the same run).
- `LOSAT_X_BPSUBJRUN=1` (BLASTP, speed; in the strict set since 2026-10-16, on by default; acts only
  where a `-subject` search gives its subjects to the pool: more than one thread and more than one
  subject) the parallel preliminary loop searches the subjects in
  runs of consecutive subjects instead of one by one (`blastp/blast_engine/x_subject_runs.rs`, exactness
  argument in its header). Without it each piece rayon splits the subjects into makes its own scratch
  (with its own diagonal table and offset-pair array) and every subject starts by filling the whole
  diagonal table (`prepare_independent_subject`; 65,536 entries of 4 bytes for a query block of 50,000
  letters). With it the subjects are cut into runs as guided self-scheduling cuts work (each run the next
  1 / (4 x threads) of the letters left, so the runs shrink to single subjects at the end and end together;
  runs of equal size made exact blastp up to 2.5 % slower at 16 threads), one task per thread takes the runs
  in order from one counter (one `rayon::scope`) with one scratch, and a run fills the table once,
  at the one-thread order's diagonal offset of its first subject; its later subjects are searched as the
  one-thread loop searches them, at the offset the previous subject's `finish_subject`
  (`Blast_ExtendWordExit`) left, which is the one-thread order's, and with every other per-subject reset.
  An entry left by an earlier subject of the run lies below the later subject's offset by at least the
  window, so it reads as a cleared entry (the property the `-db` path uses for its blocks); the reset at
  INT4_MAX / 4 (offset back to the window, table cleared) happens inside a run as in the one-thread loop.
  Results are merged by subject index as without the switch. Memory: one scratch per thread instead of one
  per rayon piece. At one thread, and with `-db`, nothing changes.
- `LOSAT_X_BPDBCUT=1` (BLASTP, speed; in the strict set since 2026-10-17 (perf/perfmb-default-bundle, evidence at the end of this entry): on unless it is `0`; acts only where a `-db` search runs its preliminary stage on the pool:
  more than one thread; every BLASTP task (blastp, blastp-fast, blastp-short), the only program whose search
  reaches the `-db` block driver) the preliminary stage of `-db` reads the database in blocks of
  clamp(sequence bytes / (4 x threads), 4 KiB, 1 MiB) instead of 1 MiB, where the sequence bytes are the letters, one per sequence and one more (the size
  of a one-volume `.psq` file) (`blastp/blast_engine/database/x_db_cut.rs`, exactness argument in its
  header). Without it a database of a few MiB is a few blocks, and a round of the parallel stage has fewer
  tasks than threads: OrthoFinder's `BlastDBSpecies6` (2.92 MB) is 3 blocks, so 3 of 16 threads search
  (12.85 s of wall for 33.7 s of CPU on the workstation). With it the same database is about 64 blocks at 16
  threads; a database of at least 4 x threads MiB keeps blocks of 1 MiB, and one thread is unchanged. This
  changes dbsearch decision (6) ("blocks of 1 MiB, 8 per thread and round") at more than one thread while the
  switch is on; the rounds of 8 blocks per thread stay. The cut only moves where blocks end: every block
  starts from a cleared diagonal table at the one-thread order's offset of its first sequence, and the
  round's results are merged in ordinal order, as for blocks of 1 MiB. A search that fails (a damaged
  database) returns the one-thread order's error: the blocks of 1 MiB are read and prepared again in order up
  to the first error, because a smaller block can meet another error first or name another byte range.
  Memory: a worker holds one block of at most the cut (or one longer sequence), and a round covers at most
  the bytes it covers without the switch; more threads search at once, each with its search scratch. Measured on the
  workstation (alternating runs, node 0, the same outputs; NCBI `-db -num_threads 1` tables for the
  OrthoFinder jobs): blastp-fast human (23,879 proteins) x `BlastDBSpecies6` at 16 threads 19.92 -> 8.69 s
  (x0.44, 4 pairs; max RSS 95.1-96.8 -> 93.5-95.8 MiB), at 8 threads 21.0 -> 13.0 s, at 4 threads 26.0 -> 21.8 s;
  exact blastp of 3,000 human proteins at 16 threads 29.5 -> 9.3 s (x0.31, 2 pairs; max RSS 37.2 -> 43-44
  MiB); OrthoFinder jobs atha x scer 281.3 -> 87.7 s and cele x scer 217.1 -> 68.4 s at 16 threads (max RSS
  55 -> 57 and 53 -> 48 MiB). User time grows for exact blastp (82.8 -> 114 s) where all 16 cores now run.
  Default ON since 2026-10-17 (perf/perfmb-default-bundle; workstation (2x Xeon Gold 6338, 16 physical cores of one
  node with their siblings idle, alternating rounds, every output the expected one; perf-mb
  `logs/defaultset_ws_report.md` unless named); a reduction is 1 - with / without the switch, the other switches of
  the candidate set on): blastp-fast of 23,879 human proteins against `BlastDBSpecies6` (`-db`, 2.9 MB) at 16 threads
  -52.8 % wall (four rounds, -50.5..-60.4 %), exact blastp of 3,000 of them -69.7 % (-69.2..-70.1 %), at 4 threads
  -11.6 / -28.1 %; exact blastp of the proteome `Species1` against it at 16 threads 49.9 s with the candidate set against 174.6 s without
  it (NCBI `-db` 75.5 s); the yeast proteome all against all (6,066 proteins, 2.9 M aa; dbsearch `REPORT_db6.md`
  follow-up) at 16 threads 32.7-35.1 s against 99.2-99.3 s without it and 60.5-61.6 s for NCBI `-db`. Max RSS -5..+9
  MiB, below NCBI's in every cell (116 vs 582, 44 vs 131 MiB). It does not act at one thread or on a database of at
  least 4 MiB per thread (Swiss-Prot keeps its 1 MiB blocks).
- `LOSAT_X_BNDBCUT=1` (BLASTN, every task, speed; in the strict set from the pull request that added it (perf/bn-db-cut,
  evidence at the end of this entry): on unless it is `0`; acts only where a
  `-db` search runs its preliminary stage on the pool: more than one thread) the preliminary stage of `-db` reads the
  database in blocks of clamp(sequence bytes / (4 x threads), 4 KiB, 1 MiB) instead of 1 MiB, where the sequence bytes
  are a quarter of the letters and one per sequence (the sequence files without their ambiguity data)
  (`blastn/blast_engine/x_db_cut.rs`, exactness argument in its header; BLASTP's `LOSAT_X_BPDBCUT` cuts by the same
  rule). Without it a database of a few MiB is a few blocks, and a round of the parallel stage has fewer tasks than
  threads: 28,036 16S rRNA genes (10.3 MB of `.nsq`) are 10 blocks, so 10 of 16 threads search, and their first 2,000
  (0.72 MB) are one block, so one thread searches. With it both are about 64 blocks at 16 threads; a database of at
  least 4 x threads MiB keeps blocks of 1 MiB, and one thread is unchanged. The cut only moves where blocks end: a
  sequence's result does not depend on what the worker's scratch searched before it (which the rounds of 1 MiB blocks
  already rely on), and the round's results are merged in ordinal order. A search that fails (a damaged database)
  returns the one-thread order's error: the blocks of 1 MiB are read again in order up to the first error, because a
  smaller block can meet another error first or name another byte range. Memory: a worker holds one block of at most
  the cut (or one longer sequence); more threads search at once, each with its scratch, as in the parallel `-subject`
  path. Measured on the workstation (node 1; 16 or 8 booked cores with their siblings idle, or one core; a 3 s warm-up
  before each run, alternating order, medians of 3 timed pairs after an untimed one, against `LOSAT_X_BNDBCUT=0` with the same binary; every
  output equal to NCBI 2.17.0 `-db -num_threads 1`; line G `results/g3ab/ab.tsv`), 100 16S queries: megablast x the
  first 2,000 16S at 16 threads 2.60 -> 0.97 s (x0.373; NCBI `-db` 3.08 s), at 8 threads 3.08 -> 1.81 s (x0.588);
  megablast x 28,036 16S at 16 threads 6.42 -> 5.84 s (x0.910; NCBI 10.26 s), at 8 threads 12.81 -> 10.06 s
  (x0.785); dc-megablast x the first 2,000 16S at 16 threads 6.19 -> 1.39 s (x0.225; NCBI 6.48 s), at 8 threads 6.52
  -> 2.46 s (x0.377); megablast of the Sakai genome (one 5.5 Mb query) x 28,036 16S at 16 threads 0.60 -> 0.54 s
  (x0.900), at 8 threads 1.10 -> 0.78 s (x0.709). One thread unchanged (x0.999-1.005). User and system time fall
  (dc-megablast at 16 threads: user 31.9 -> 18.6 s, system 8.71 -> 0.37 s), and the median max RSS falls in every cell (139-143
  -> 119-122 MiB and 147-151 -> 117-121 MiB for megablast at 16 threads), below NCBI's except for the Sakai query at
  16 threads: 197-203 MiB with the switch and 200-203 MiB without it, against NCBI's 194 MiB.
- `LOSAT_X_BNDBSTREAM=1` (BLASTN, every task, memory; in the strict set from the pull request that added it
  (perf/bn-db-stream, evidence at the end of this entry): on unless it is `0`;
  acts only where a `-db` search runs its preliminary stage on the pool: more than one thread) the parallel preliminary
  stage of `-db` has no rounds: the threads claim the blocks once each in ordinal order, and one of them hands each
  block's results to the hit lists in that order while the others search the next blocks, at most 2 blocks per thread
  past the last one handed over (16 with `LOSAT_X_BNDBWIDE`; `blastn/blast_engine/x_db_stream.rs`, exactness argument in its header; the driver is
  `x_bnstream::search_in_order`, that of the `-subject` path's `LOSAT_X_BNSTREAM`, here with its `assist`: a thread
  that waits for room in the window, or the owner for the next block, runs the pool's pending jobs first, such as the
  chunks of a long sequence that another thread's search hands out). The blocks are listed 65,536 at a time. Without
  it the stage searches rounds of 8 blocks per thread and merges a round's results when the whole round is searched,
  so the preliminary HSP lists of all sequences of the round are in memory together; blocks hold at most 1 MiB, so
  every database below 8 x threads MiB of sequence bytes is one round (128 MiB at 16 threads). Megablast of 100 16S queries x the 28,036 16S rRNA
  genes at 16 threads: the round held 69 MiB of preliminary HSP lists, against 16 MiB for the lookup table, and the
  arenas kept the freed space to the end (103 MiB free in them at exit; probe of the live allocations by thread and
  region, v4-B session 10). The order of the merges is the one-thread order, every block is searched once, and a
  thread keeps one scratch for the blocks it claims in increasing order (a sequence's result does not depend on what
  the scratch searched before it, which the rounds' `map_init` already relies on); a failed search returns the first
  error in ordinal order after merging the blocks before it, as a round does, and blocks not yet claimed are not
  searched. Memory measured on the workstation (node 0, untimed, v4-MEM's `mem_run.py`: VmHWM and RssAnon every 50 ms;
  both arms and NCBI 2.17.0 `-db` in the same run; every output equal to NCBI's `-num_threads 1` output; v4-B
  `results/b14ab/mem.tsv`), megablast of 100 16S queries x 28,036 16S, max RSS with / without the switch (NCBI): 16
  threads outfmt 6 77-78 / 114-115 MiB (141), outfmt 8 81-82 / 117-120 (108), outfmt 9 78-79 / 118-119 (107), outfmt
  12 79-83 / 120-123 (109); 8 threads outfmt 6 70-74 / 105-108 (146), outfmt 8 77-78 / 115-117 (114), outfmt 9 68-75 /
  108-110 (115), outfmt 12 79-80 / 113-118 (115); one thread unchanged (33 MiB with outfmt 8, NCBI 112). Time measured
  on the workstation (node 1; 16, 8 or 2 booked cores with their siblings idle, a 3 s warm-up before each run, alternating
  order, against `LOSAT_X_BNDBSTREAM=0` with the same binary; every output identical; v4-B `results/b14ab/ab3.tsv`,
  `ab4.tsv`, `ab5.tsv`), 100 16S queries x 28,036 16S at 16 threads: megablast (outfmt 8) 5.885 -> 5.055 s (x0.859),
  dc-megablast 12.76 -> 10.17 s (x0.797), blastn 13.01 -> 10.56 s (x0.812); megablast at 8 threads 10.18 -> 9.41 s
  (x0.924), at 2 threads 36.47 -> 35.97 s (x0.986); megablast of the Sakai genome x 28,036 16S at 16 threads 0.575 ->
  0.490 s (x0.852; max RSS 199 -> 174 MiB). User-space instructions fall where the round held many results
  (dc-megablast 1,570 -> 1,500 G at 16 threads): one scratch per thread, and no merge after the round. A database of
  a few long sequences is unchanged (the D. melanogaster genome, 1,870 sequences with arms of 23-32 Mb: blocks of one
  sequence and the nested chunk stage, which the waiting threads join): 1,500 2L cDNAs, blastn per-pair median x1.009
  (9 pairs, 0.971-1.033) at 16 threads and x1.001 (8 pairs) at 8, dc-megablast x0.995 (7 pairs) at 16; max RSS equal
  with and without the switch (blastn at 16 threads 568-640 / 572-629 MiB, which is above NCBI's 382 MiB either way: a
  per-thread excess this switch does not touch).
- `LOSAT_X_BNDBWIDE=1` (BLASTN, every task, speed; in the strict set from the pull request that added it
  (perf/bn-db-wide2, evidence at the end of this entry): on unless it is `0`; acts only where `LOSAT_X_BNDBSTREAM`'s
  window runs: `-db`, more than one thread, and a list of more blocks than 2 per thread) the window of
  `LOSAT_X_BNDBSTREAM` is 16 blocks per thread (`X_BNDBWIDE_BLOCKS_PER_THREAD`), and it narrows to the old 2 per thread
  while more than 1 MiB per thread of results wait for the owner (`x_bnstream::WaitingBound`, the same 2 per thread
  keep the head moving; `blastn/blast_engine/x_db_stream.rs`). The old window was shorter than the blocks: on core_nt
  at 16 threads its 32 blocks were about 24 ms of the wall clock, and a block took 54 ms at p99, so the threads ran out
  of room behind every slow block. They waited 3,344 s of 10,860 s of thread time, evenly over the run (11.8 of the 16
  CPUs busy in every 30 s bin; 15.8 for a window of 16 per thread without the bytes, and for NCBI), and the owner's merges took only 11.9 s (probe
  perf/mb-wait-probe, `dropin-survey-20261029/BIGDB/STATE.md`). Which thread searches a block changes, and the order of
  the merges does not, so reports are unchanged for the reason given under `LOSAT_X_BNDBSTREAM`. A list of 2 blocks per
  thread or fewer behaves as without the switch. Reach: recorded with the window; the fixtures' databases do not reach
  the new side, `LOSAT/tests/blastn_db_stream.rs` checks both sides and the reports at 2, 4 and 8 threads, and
  `x_bnstream`'s `the_bytes_narrow_the_window_to_its_floor` the bound.
  Evidence (ws, node 1, booked CPUs with their siblings idle, one booking per round, three rounds with the arm order
  rotated, medians, every output equal per case and equal to NCBI's where it ran; dev caf4c166e against `=0` and on, both
  built from da1f0db42; `dropin-survey-20261029/BIGDB/s2/ab1.tsv`), `draft200` with the options of BlobToolKit's
  megablast, dev / `=0` / on (max RSS, MiB, dev / `=0` / on): core_nt (92 volumes, 298 GB, 575,156 blocks) at 16 threads
  651.9 / 657.6 / 451.0 s (844-960 / 861-924 / 841-852), against NCBI 2.17.0 516.8 s (x0.873 of NCBI; 257 GiB), with user
  7,580 -> 7,068 s and system 277 -> 94 s; ref_prok_rep_genomes (27 volumes) at 16 threads 45.95 / 47.35 / 42.29 s
  (549-559 / 549-577 / 452-459), against NCBI 48.59 s (x0.870; 25 GiB); at 8 threads 87.38 / 87.25 / 82.28 s (322-353 /
  331-351 / 283-297); at 32 threads 30.02 / 30.00 / 27.25 s (908-974 / 964-981 / 759-779); its first volume at 2 threads
  14.11 / 14.11 / 13.83 s (129-135 / 129-135 / 130-133). A small database, 100 16S queries x the 28,036 16S at 16 threads,
  is the same within the rounds' spread: megablast 4.43 / 4.42 / 4.44 s (68-70 / 68-69 / 68-72), dc-megablast 9.68 / 9.70 /
  9.61 s (62-63 / 61-63 / 61-67), blastn 9.94 / 9.95 / 9.97 s (48-49 / 47-49 / 48-51). RSS falls where the waits fall.
  Results pile up only behind a stall: untimed on node 0 (`perf/bn-db-wide-tune`, two rounds, 5 arms, 16 threads), the
  results waiting at their peak were 28-36 MiB on the 16S with the old window as well, and max RSS was flat across
  windows of 2, 8 and 16 blocks per thread with budgets of 1 MiB, 4 MiB and none. On databases made dense with hits
  (draft200 or 100 16S queries, 16 threads on node 1, two rounds, the same output; `dropin-survey-20261029/BIGDB/s3/dense/`,
  2026-10-31), bound / no bound / `=0` / NCBI: 1,600 mutated Sakai copies in 1 kb contigs with seven slow stretches of a repeated
  1 kb window 25.6 / 24.7 / 27.8 / 33.5 s, max RSS 114 / 133 / 100 / 3,090 MiB, results waiting at their peak 22 / 43 / 9.7 MiB. 50
  mutated copies of the 28,036 16S 190 / 189 / 192 / 307 s, 131 / 125-134 / 130 / 671 MiB, 62-70 / 49-83 / 67-81 MiB; with a long
  subject every 140,000 genes 209 / 209 / 210 / 422 s, 317 / 324-332 / 317-330 / 913-930 MiB, 54-62 / 56-64 / 59-62 MiB. The bytes
  act only behind a stall on a database of many hits, and there they cost 3 % of the time for 14 % of the memory; without
  them no run held more than 83 MiB. One thread does not run the window: one core of node 0, the first volume
  dev / `=0` / on 27.45 / 27.20 / 27.34 s, 110 MiB in each arm (three rounds; `dropin-survey-20261029/BIGDB/s2/n1_ab.tsv`).
- `LOSAT_X_BNDCONE=1` (BLASTN `-task dc-megablast`, speed; in the strict set from the pull request that added it
  (perf/dc-one-thread, evidence at the end of this entry): on unless it is `0`; acts only where a chunk in scope has
  no pool of two threads or more: a search on one thread; nothing while
  `LOSAT_X_PAIRPARDISC` / `LOSAT_X_PAIRPARDISCARR` is `0`) the seed stages of `LOSAT_X_PAIRPARDISC` (diagonal hash) and
  `LOSAT_X_PAIRPARDISCARR` (diagonal array) also run without a pool: a chunk in their scope (a subject this thread searches by itself, unmasked, `scan_range` 0, no
  debug or trace) is scanned in their strips on this thread and folded in strip order, as the runs of their shadows
  without a pool already were (`blastn/blast_engine/x_pair_seed_disc.rs`); megablast given `-template_type` (a
  discontiguous table) takes them too. Before, such a chunk took the reference fold (`process_offset_pairs`), which
  keeps neither the memo of `LOSAT_X_DIAGLV2` nor the lean folds of `LOSAT_X_FOLDLEAN` (it reads the CSR table of
  `LOSAT_X_DCCSR` as the stage does). Exact for the stages' reason: the owner folds every
  pair in the reference order with the reference statements whatever the number of threads (their entries). Many subjects on a pool keep the reference
  fold in every worker (outside the scope). Reach: recorded where a chunk in the scope has no pool. Tests:
  `LOSAT/tests/blastn_pairpar_disc.rs` and `blastn_pairpar_disc_arr.rs` (`the_stage_runs_at_one_thread`: the same
  bytes as the reference, the stage's chunks counted, none with the switch `0`), `blastn_diag_lv.rs` (the memo
  answers at one thread with the switch).
  Evidence (ws, node 1 timed, dev c9303778 against perf/dc-one-thread 84c92737, medians of the per-pair ratios, pairs
  in rotated order, every output dev's; v4-A15 pr4, megablast pr6 and mbctr_s8): one thread, 31,168 dmel cDNAs against 3 Mb of dmel genome with
  `-subject` 20.06 -> 14.50 s (x0.723, 6 of 6 pairs) and as `-db` 20.08 -> 14.47 s (x0.721); 1,775 cDNAs x0.730 (20
  pairs); 236 cDNAs x0.836 (60); 7 cDNAs of 6,993 letters (the diagonal array) x0.843 (80); one 16S against a 16S
  `-db` (the array, 28,036 chunks) x0.953 (30); 10 16S against it (the hash, 56,072 chunks of one subject each) x0.977
  (8); the switch at `0` x0.998 to x1.001. Where it does not act: 100 16S against the 16S `-db` at 16 threads x0.999,
  the 31,168 cDNAs at 16 threads x1.001, megablast of the 31,168 cDNAs at one thread x1.0015 with the switch
  at `0` x1.0017 (16 pairs each) and dev's instructions (-0.002 %: the code layout, the hot functions 320 to 768 bytes
  further on in their pages). Memory (peak RSS, the same runs): one thread on the 31,168
  cDNAs `-subject` 85.0 -> 86.3 MiB and `-db` 83.3 -> 85.6 MiB (`MALLOC_ARENA_MAX=1`: 85.1 -> 86.7 and 83.6 -> 85.3
  MiB; NCBI at one thread 313.4 and 151.0 MiB), the other cases within 1 MiB. After the review fixes (6b9281b6, v4-A15
  pr7, the same dev): one thread on the 31,168 cDNAs x0.723, 236 cDNAs x0.836, the 7 cDNAs x0.844, 236 cDNAs against
  15,000 200-bp pieces of the genome as `-subject` x0.994 (16 pairs), megablast with a coding 18 template (a
  discontiguous table) x1.000; the one-thread slots of `LOSAT_X_QBATCHPAR` (4 threads with `LOSAT_X_QBATCHPARSLOTS=4`; slots, removed on 2026-10-28)
  on the 31,168 / 1,775 cDNAs x0.726 / x0.732 (peak RSS +13 MiB on the 31,168, the four slots' strips), by default at
  4 threads x1.001 / x0.999; megablast at one thread x0.9997 with the switch at `0` x0.9996 (the layout difference
  above is gone in this binary), 16 threads unchanged.
- `LOSAT_X_BNFEWSUBJ=1` (BLASTN `-subject`, every task whose seed stage has strips (PAIRPAR,
  `LOSAT_X_PAIRPARDISC`, `LOSAT_X_PAIRPARDISCARR`), speed; in the strict set from the pull request that added it
  (perf/bn-few-subj, evidence at the end of this entry): on unless it is `0`;
  acts only on a pool of 2+ threads) a `-subject` search whose longest work unit is a large part of its letters searches
  its subjects one at a time on the whole pool, and a subject's batch of chunks the same way its chunks one at a time, so
  that the seed stage scans each chunk in strips (`blastn/blast_engine/x_few_subj.rs`, exactness argument in its
  header). Without it a search of two subjects or more searches each subject as one task of the pool, and a subject of
  two 5 Mb chunks or more (`MAX_DBSEQ_LEN`) its chunks as tasks, both without strips: the Sakai genome against MG1655 and
  its 92 kb plasmid scanned MG1655 on one thread. The rule: the longest unit's letters x S^(p) > all the letters, the
  unit a subject up to 5 Mb (subjects) or a chunk with its overlap (a batch), S^(p) per kind of table at the largest of
  2 / 4 / 8 / 16 threads not above p: dc-megablast 1.7 / 2.5 / 3.2 / 3.4, blastn 1.7 / 2.9 / 4.1 / 4.7, any other
  task 1.5 / 2.0 / 2.5 / 2.7 (the lowest one-thread / new-shape wall the probe below measured, floored to tenths;
  integers in tenths). The subjects keep their tasks with `-lcase_masking` (a masked chunk has no strips), and a
  batch's chunks keep theirs when the subject has lowercase letters under it; never for `-db` (plan v4 line D). The
  new shape is the one of a one-subject search (strips, `LOSAT_X_QGROUPPAR`, the speculative traceback), and the
  subjects reach `sink` in OID order in both shapes. It has no shadow of its own (it only chooses between
  two shapes with the same bytes). Memory: the shape holds one chunk's seed stage on
  the whole pool. The pull request that added the switch also cut the one-hit strips of `LOSAT_X_PAIRPARFOLD` to 64
  per thread with a look-ahead of 2 per thread (about 1/32 of a chunk's seeds, as dc-megablast's; before 8 / 2, a
  quarter), for every search and not under this switch (`=0` keeps them). That also lowers one-subject blastn (Sakai
  against MG1655 alone at 16 threads 282 -> 213 MiB median, at 4 threads 193 -> 151; NCBI 192 MiB) at x0.95-1.00 of
  the time. Test: `LOSAT/tests/blastn_few_subj.rs`.
  Measured on the workstation (node 1; 16 / 8 / 4 / 2 booked cores with their siblings idle, a 3 s warm-up before each
  run, rotated order; every output equal to NCBI 2.17.0 `-subject`; line G `results/g4probe/ab.tsv` (the shape forced
  for any 2+ subjects or chunks, exp/bn-few-subj-force 80fca36a, 3 timed pairs) and `results/g4conf2/ab.tsv` (this
  binary, 8 timed pairs)), the Sakai genome (5.5 Mb) as the query, with the switch against `LOSAT_X_BNFEWSUBJ=0`:
  against MG1655 + pO157 dc-megablast 1.88 -> 0.51 s at 16 threads (x0.271; NCBI 8.55 s), x0.330 / 0.485 / 0.716 at 8 /
  4 / 2; blastn 4.66 -> 1.40 s (x0.303; NCBI 22.17 s), x0.349 / 0.476 at 8 / 4; megablast 0.34 -> 0.25 s (x0.735; NCBI
  1.41 s), x0.911 at 2; against EDL933 (5.53 Mb, two chunks) dc-megablast x0.423 / 0.456 / 0.565 / 0.778 at 16 / 8 / 4 /
  2, blastn x0.275 / 0.313 / 0.438 / 0.679, megablast x0.711 / 0.725 / 0.800 at 16 / 8 / 4 and x0.907 at 2 (this
  binary); against MG1655 + EDL933 (longest
  unit 0.49) dc-megablast x0.493 / 0.557 / 0.734 at 16 / 8 / 4, blastn x0.416 / 0.501 / 0.704, megablast x0.881 / 0.889
  at 16 / 8 (this binary); against MG1655 in 3 records dc-megablast x0.614 / 0.706 and blastn x0.679 / 0.770 at 16 / 8; 75,000
  reads of 150 bases against MG1655 + pO157 megablast x0.667 / 0.841 at 16 / 8; dc-megablast of 2,000 pieces of 1 kb of
  Sakai against MG1655 + pO157 0.66 -> 0.13 s at 16 threads (x0.197; NCBI 3.07 s), x0.434 at 4 (this binary). The
  probe's blastn and megablast cells ran the strips before the 64 / 2 change, which takes x0.95-1.00 of their time. The
  cells it leaves
  as they were: MG1655 in 10 records (the new shape x1.01-1.48 slower), 2,000 16S genes (x5.2 slower), MG1655 + EDL933
  at 2 threads (x0.94-1.06) and with megablast at 4 (x0.96), MG1655 in 3 records at 2 / 4 threads (x0.86-0.96) and with megablast (x0.94-1.04). Max RSS
  (median, MiB, with / without; NCBI): blastn against MG1655 + pO157 203 / 195 at 16 threads (192), 180 / 151 at 8,
  150 / 125 at 4; megablast against MG1655 + EDL933 149 / 161 at 16; dc-megablast of the 2,000 pieces 66 / 64 at 16
  (125). On node 0 (untimed, 4 cores for 16 / 8 threads, two runs each, `results/g4rss3/rss.tsv`; these read the old
  shape about 20 % below the probe's 16 cores) the cells with EDL933 or 3 records, at 16 threads: blastn against
  MG1655 + EDL933 420-484 / 494-539 MiB, against EDL933 616-647 / 671-703, against MG1655 in 3 records 198-199 /
  220-236; dc-megablast against MG1655 + EDL933 658-711 / 639-768, against EDL933 671-672 / 665-689; at 8 threads the
  new shape's runs are within or above the old one's (blastn against MG1655 + EDL933 378-397 / 353-383,
  dc-megablast against EDL933 449-477 / 418-460). With EDL933, blastn and dc-megablast take 2.3-4.4 times NCBI's max
  RSS (164-286 MiB) in both shapes (outside this switch); megablast stays below NCBI's.
- `LOSAT_X_BNMAPPED=1` (BLASTN, every task, memory; in the strict set from the pull request that added it
  (perf/bn-mapped-pool, evidence at the end of this entry): on unless it is `0`;
  acts only with 64-bit glibc (Linux) and on buffers of 1 MiB or more) the large per-subject buffers of BLASTN are
  kept where glibc returns them to the system when they are dropped: the blocks and sequences read from a `-db`
  (`seqdb` `x_read_block_mapped`, `x_read_sequence_mapped`), their BLASTNA letters (`blastn/blast_engine/database.rs`
  `x_blastna`), the ungapped hits of a subject as the seed stages push them (`run.rs` `x_push_ungapped`) and the
  scratch of their sort by score (`utils/x_mapped.rs` `x_sort_by_mapped`). Such a buffer gets a capacity above glibc's
  largest dynamic mmap threshold (32 MiB and one page) once it needs 1 MiB, so glibc maps it; its pages beyond its
  length are not resident until they are touched. Without it glibc raises its threshold to the size of each mapped
  block that is freed, so after a few frees the chromosome-sized buffers of the worker threads come from their arenas,
  and the space they leave stays resident there for that thread only: the D. melanogaster genome `-db` (arms of 23-32
  Mb) at 16 threads had 605 MB free in the arenas at the peak (probe of the live allocations, v4-B session 11). A
  mapped block that is freed and mapped again faults on every page (x1.04 wall at 16 threads there), so the buffers
  that are done after every subject go to pools shared by all threads (`utils/x_mapped.rs` `XPool`), from which the
  next buffer that must grow past 1 MiB is taken with its pages resident: the database's blocks, sequences and BLASTNA
  letters to `X_BYTES` (idle bound `LOSAT_X_POOLIDLE`, below), the sort's scratch to `X_SORT` (`run.rs`, idle bound
  `LOSAT_X_SORTIDLE`, below). The hits are mapped but not pooled: a thread keeps its scratch, and with it its hit
  buffer, from one subject to the next (as NCBI keeps its `init_hitlist`), so idle hit buffers in a pool only added to
  those in use (Sakai x EDL933 `-subject` blastn at 4 threads: max RSS 298-327 MiB with a pool of 64 MiB, 245-264
  without, 238-261 with the switch off). The sort is a stable merge sort with the same comparator (runs of 64 KiB
  sorted with `sort_by`, merged pairwise into the scratch and back, the left run first among equal elements), so the
  order is that of `sort_by`; its scratch is as long as the hits (a half-size scratch, the left run copied out before
  each merge, held no less at the peak).
  Every pooled buffer is written before it is read (zeroed and read into, or BLASTNA letters written whole), so values
  and order never change; it only chooses where the bytes live and when they are freed. Not mapped: the merge of the
  folded strips (`x_pair_seed_part.rs`) and the per-query sorts of `LOSAT_X_QGROUPPAR`. With other allocators (macOS,
  Windows, musl, Wasm) and on 32 bits (glibc's ceiling is 512 KiB there) the capacities are those of `Vec` and nothing
  is pooled. No shadow (allocation only). Tests: the unit tests of `utils/x_mapped.rs` (capacities, the sort against
  `sort_by` for 0-200,003 elements and runs of 4, the pools' lending and idle bound); the gate runs every BLASTN
  fixture with `=0` as well.
  Memory measured on the workstation (node 0, untimed, v4-MEM's `mem_run.py`: VmHWM every 50 ms; with the switch and
  without it in the same run, rotated order; every output equal to NCBI 2.17.0's; v4-B `results/b16ab/mem8.tsv`,
  `mem9.tsv` with the max RSS of the timed runs `ab8.tsv`, `ab9.tsv`; 5-23 runs per arm), max RSS median with /
  without the switch (NCBI): 1,500 D. melanogaster 2L cDNAs x the D. melanogaster genome `-db`, blastn at 16 threads
  321 / 624 MiB (372), at 8 threads 272 / 343 (336); dc-megablast at 8 threads 186 / 222 (193), at 16 threads 268 / 284
  (191; above NCBI either way: the traceback decoded whole chromosome arms, NCBI fetches only the HSPs' ranges, as
  `LOSAT_X_BNTBRANGE` below now does); the
  Sakai genome x EDL933 `-subject` (blastn at 4 threads hit lists of 1.0-11.4 MiB, dc-megablast below 1 MiB), blastn at
  4 threads 264 / 264, at 16 threads 585 / 586, dc-megablast at 4 threads 294 / 292, at 16 threads 686 / 697 (equal
  within the run-to-run spread of 20-60 MiB; NCBI 177 / 164 MiB on its one `-subject` thread: the long tracebacks of
  several threads at once, outside this switch); megablast of 100 16S x 28,036 16S `-db` at 16 threads 85 / 87 (139).
  Time measured on the workstation (node 1; 16, 8 or 4 booked cores with their siblings idle, a 3 s warm-up before each
  run, alternating or rotated order, against `LOSAT_X_BNMAPPED=0` with the same binary; every output identical; v4-B
  `results/b16ab/ab8.tsv`, `ab9.tsv`), per-pair median: the D. melanogaster cells blastn at 16 threads x0.999 (12 pairs,
  0.980-1.046), at 8 x0.992 (10), dc-megablast at 8 x1.011 (12, 0.977-1.038); Sakai x EDL933 blastn at 16 threads
  x1.008 (16, 0.964-1.092), at 4 x0.998 (14), dc-megablast at 16 x0.996 (12), at 4 x1.008 (14, 0.988-1.016); megablast
  of the 16S cells at 16 threads x1.000 (4). User-space instructions are equal; system time rises where buffers are
  mapped (D. melanogaster blastn at 16 threads 2.5 -> 3.3 s, inside the same wall).
- `LOSAT_X_POOLIDLE=<MiB>` (tuning variable of `LOSAT_X_BNMAPPED`; not in the strict set, read once when the pool
  first takes a buffer back) the bound of the bytes that the idle buffers of `X_BYTES` (the database's blocks,
  sequences and BLASTNA letters) may have touched, 0 to 65,536. The oldest idle buffers beyond it are freed, and a
  buffer that alone touched more is freed at once. Unset or any other value gives 64, without a message, as the value
  switches of `utils/xswitch.rs` give their built-in value. The output does not depend on it. D. melanogaster blastn
  at 16 threads (above): max RSS 273 MiB with 0, 317-331 with 64, 346-367 with 256; with 0 every reuse faults again
  (x1.04 wall).
- `LOSAT_X_SORTIDLE=<MiB>` (tuning variable of `LOSAT_X_BNMAPPED`; not in the strict set, read once when the pool
  first takes a scratch back) the same bound for `X_SORT`, the scratch of the ungapped hits' sort, 0 to 65,536. Unset
  or any other value gives 16 (one scratch of the largest hit lists measured: 11.4 MiB, Sakai x EDL933 blastn),
  without a message. The output does not depend on it. Sakai x EDL933 blastn at 4 threads (above), per-pair against
  16: 8 x1.009 (6 of 8 pairs slower; the scratches of 9.6-11.4 MiB are mapped again at every sort), 0 x1.009-1.027,
  for 10 MiB less max RSS (median 251 against 261 MiB); elsewhere no change in memory.
- `LOSAT_X_BNTBRANGE=1` (BLASTN `-db`, every task, memory and time; in the strict set from the pull request that added
  it (perf/bn-tb-range, evidence at the end of this entry): on unless it is `0`)
  the traceback of a database search (`blastn/blast_engine/database.rs` `traceback_listed`) decodes a sequence of
  1 Mi bases or more to BLASTNA only in the ranges of the HSPs the preliminary stage kept of it, as NCBI's partial fetch
  does: each HSP widened by 1024 bases on both sides within the sequence, sorted by start, and merged into the range
  before unless it starts more than 1024 bases after its end (`x_partial_ranges`; `BLAST_SetupPartialFetching`,
  blast_traceback.c:1320-1370, blast_seqsrc.c:467-519), with `FENCE_SENTRY` (201) before and after each range inside
  the sequence (`seqdb` `x_nucleotide_blastna_ranges`; `GetAmbigPartialSeq`, seqdbvol.cpp:1541-1600). The buffer has
  the whole length and is zeroed by the allocator; with 64-bit glibc its capacity is above the mmap ceiling, so it is
  mapped and its pages outside the ranges are never written or resident (other allocators may zero it whole; Wasm's
  dlmalloc does, so there the switch saves the decoding only). The traceback kernels stop at a fence as NCBI's do; their
  fence branches count on the thread (`blastn/alignment/x_fence.rs`), every traceback call compares the count around
  it and marks the sequence (`BlastnSubjectRef::x_fenced`, also from the look-ahead helpers and the `LOSAT_X_QGROUPPAR`
  tasks), the traceback loop makes no HSP once it is marked, and `traceback_listed` decodes the whole sequence and traces
  it again from the same kept HSPs (NCBI refetches the whole sequence and traces the list again from its copy,
  blast_traceback.c:307-310, 516-520, 1651-1688). An alignment that met no fence read only decoded bytes, which are the
  whole decode's (the ambiguity runs applied in stored order, each to the ranges it overlaps), so values and order never
  change. NCBI turns partial fetching on per database (longest sequence at least 5000, mean length at least 2048,
  seqsrc_seqdb.cpp:234-260); LOSAT decides per sequence (a shorter one is decoded whole: its BLASTNA is small). Test
  variables (not in the strict set, read once): `LOSAT_X_BNTBRANGEMIN=<bases>` the shortest sequence decoded in ranges
  (0 = every sequence; default 1,048,576), `LOSAT_X_BNTBRANGEPAD=<bases>` NCBI's 1024 (16 to 1,048,576: the start of a
  gapped alignment reads up to 11 bases on each side of the seed outside the kernels' fence tests), so that the tests
  meet fences; `LOSAT_X_BNTBRANGELOG=1` prints one line at exit (sequences in ranges, traced again whole, their bases,
  the bases of their ranges). Tests: the unit tests of `database.rs` (NCBI's ranges; the decode in ranges against the
  whole decode on the `nucl_mv` fixture) and `seqdb/volume.rs` (the runs over ranges against the whole decode with
  shuffled runs, offsets of 2^31 and more and runs past the end); the gate runs every BLASTN fixture, the `-db` fixtures
  at 1 and 4 threads, `check_losat` and the capture with `LOSAT_X_BNTBRANGEMIN=0` (all sequences in ranges, no
  sequence traced again) and with `LOSAT_X_BNTBRANGEPAD=16` too (279 of the 371 sequences of the `-db` fixtures meet a
  fence and are traced again), all equal to NCBI 2.17.0's bytes; the shadow found no difference.
  Measured on the workstation (v4-B session 14, `results/b16b/`; every output equal to NCBI 2.17.0's, the same with
  and without the switch). 1,500 D. melanogaster 2L cDNAs x the D. melanogaster genome `-db` (62 of its sequences of
  1 Mi bases or more traced by blastn, 1,236.6 M bases, 79.9 M of them in ranges; dc-megablast 10, 216.7 M, 5.0 M; no
  sequence traced again). Max RSS (node 0, untimed, v4-MEM's `mem_run.py`, rotated order, 3 runs each, `mem1.tsv`),
  median with / without the switch (NCBI): dc-megablast at 16 threads 129 / 281 MiB (188), at 8 threads 111 / 187 (193);
  blastn at 16 threads 205 / 331 (388), at 8 threads 146 / 273 (361); megablast and dc-megablast of 100 16S x 28,036 16S
  `-db` at 16 threads (no sequence of 1 Mi bases) 84 / 81 (139) and 66 / 66 (120). Time (node 1; 16 or 8 booked cores
  with their siblings idle, a 3 s warm-up before each run, alternating order, against `LOSAT_X_BNTBRANGE=0` with the
  same binary, `ab1.tsv`), per-pair median: dc-megablast at 16 threads x0.918 (8 pairs, 0.903-0.953), at 8 x0.904 (8,
  0.868-0.957); blastn at 16 x0.854 (6, 0.834-0.865), at 8 x0.866 (6, 0.861-0.886); megablast of the 16S cells at 16
  x0.995 (4, 0.979-1.017); again on the pull request's code (`ab2.tsv`): megablast of the D. melanogaster cells at
  16 threads x0.815 (8, 0.788-0.844; max RSS 263 against 313 MiB), at 8 x0.845 (6, 0.806-0.886; 263 against 308),
  dc-megablast at 16 x0.928 (6, 0.865-0.946). User-space instructions x0.958-0.968 (megablast x0.73): the whole
  sequences' decoding is gone.
- `LOSAT_X_DBIDXREAD=1` (every program's `-db`, memory; in the strict set since 2026-10-19 (perf/db-index-pread, evidence
  at the end of this entry): on unless it is `0`) the index file of each volume
  (`.pin` / `.nin`) stays open after the volume is opened, and the reader reads the ranges of its offset arrays that it
  needs by position (`seqdb/index.rs` `IndexFile::open`): a block's sequence (and ambiguity) offsets with one read when
  the block is read, a header's two offsets when a title is read, the offsets after a block's first sequence in windows
  of 1,024 entries and more when the next block is planned (the end of the volume first, one entry), and the block's
  offsets for BLASTP's diagonal offset of the next block. Without it the whole index file is read into a buffer when the
  volume is opened and the arrays are held until the end of the run: 4 bytes per sequence for the header offsets and 4
  for the sequence offsets (and 4 for the ambiguity offsets of a nucleotide volume), 296 MiB for UniRef50 (38.8 M
  sequences, five volumes) with a 133 MB buffer of the largest volume while it is checked. The arrays are checked as
  before when the volume is opened (ascending offsets, the protein sentinel bytes, the nucleotide ambiguity offsets
  inside their sequence's data, the arrays inside the file), in one pass that reads 1 MiB at a time, with the same
  errors in the same order. The checks of every title that the block search of `LOSAT_X_BPQBLOCK` and the report
  streaming of `LOSAT_X_RPTSTREAM` (blastn) run before they start read the headers in ranges of 1,024 sequences (one read
  of the index file and one read of the header file per MiB, `SeqDb::visit_deflines`) instead of two reads per sequence,
  and blastn's ambiguity check reads the offsets of the same ranges at once (`SeqDb::visit_ambiguity`), the ambiguity
  data of the sequences that have any as before. Every value read is the one of the held arrays, so blocks, reports and
  errors do not change; a read of the index file that fails after the volume was opened is a LOSAT error of the read of
  the block (or title) that needs it. After the open, every range of the arrays read again must keep the order the open
  checked (the header and sequence offsets ascending, each ambiguity offset inside its sequence's data): an index file
  changed in place during a run is a database format error ("... differ from those checked at open (the file changed)")
  instead of a subtraction that underflows. Like NCBI (which keeps the index file mapped, `seqdbfile.hpp` `m_Lease`), a
  volume now holds three open files (index, header, sequence) instead of two. Evidence 2026-10-19 (perf/db-index-pread,
  ws node 0, untimed memory with `RssAnon` sampled every 0.2 s; `blastp -db` UniRef50 outfmt 6, every output equal to
  NCBI 2.17.0's): peak anonymous memory with the switch / dev 0797cc40 / NCBI: q1 n1 3.5 / 300 / 8.8 MiB, q1 n8
  18.5 / 315 / 72.6 MiB, q10 n8 44.4 / 340 / 99.6 MiB, q100 on volume `uniref50.04` n4 38.6 / 88.0 / 84.5 MiB, q100 on
  `uniref50.03` n8 104.5 / 222.8 / 133.5 MiB (max RSS q1 n8 24 / 379 MiB, NCBI 11.6 GiB with the mapped files); the
  same binary with `LOSAT_X_DBIDXREAD=0` 346 MiB on q1 n8. Timed (same binary on / off, 8 cores of node 0 with their
  siblings idle, warmed, order alternated): median ratio of wall time x1.005 (0.965-1.073, 8 pairs) OrthoFinder
  blastp-fast Species4 x BlastDBSpecies6, x1.000 (0.996-1.020, 8) Swiss-Prot q100, x1.004 (0.998-1.014, 8) megablast
  16S q100, x0.971 (0.940-1.003, 4) `uniref50.04` q100 (max RSS 53-56 against 104 MiB). Gates: ci_fast all programs
  with the switch, with `=0` and with `LOSAT_X_DEFAULT=0` (223 cases, 0 failures each); `-db` fixtures of every program
  at 1 and 4 threads under the shadow (0 differing, the exit line `identical` in every run that ends without an error).
- `LOSAT_X_BPLISTEDLEAN=1` (BLASTP `-db`, memory; in the strict set from the pull request that added it (perf/bp-db-tail,
  evidence at the end of this entry): on unless it is `0`; acts for every task
  with `-db`, in the search of all queries and in the blocks of `LOSAT_X_BPQBLOCK`, not with ABI v1's titles) the reports of a `-db` search read the letters of
  the listed subjects (every subject of the preliminary hit lists, `ListedSubjects`) from their stored residues and hold
  the titles of the reported subjects only (`blastp/blast_engine/database/x_listed_lean.rs`). Without it
  `ListedSubjects::read_titles` makes one record per listed subject, with its letters in ASCII and its title (empty for a
  subject in no final hit list), and outfmt 0 and the `stitle` field each copy every title into an `Arc`. With it the
  same titles are read in the same order with the same errors into one `Arc` per reported subject (the others share one
  empty title), outfmt 0 and `stitle` borrow them, and the two readers of letters (the recount of the HSPs with
  pyrrolysine, `recount_pyrrolysine_identities`, and the O of rendered rows, `restore_displayed_o`) read a stored residue
  as `ncbistdaa_to_ascii` gives it, which is how the records' letters are made (`blast_engine.rs` `ShownSubjects`), so
  reports and errors do not change. Evidence 2026-10-20 (perf/bp-db-tail, workstation,
  probe build): `blastp -db` UniRef50 q100 at 32 threads lists 70,701 subjects (50.8 MB of letters) for the first batch
  of 91 queries; the records added 41 MiB and the two title copies 10 MiB of anonymous memory in the main heap.
  Measured (head de4b6401, workstation): that run on / off in one job (node 1) peak RssAnon 286 / 337 MiB (max RSS
  282 / 327 MiB; NCBI 195 MiB), wall 538 / 543 s, the same output (= NCBI). Instructions (`perf stat` `:u`, 8 threads,
  on / off) x0.997 OrthoFinder blastp-fast Species4 x BlastDBSpecies6, x0.9995 Swiss-Prot q100 outfmt 6, x0.9993 the
  same with outfmt 0, x1.000 `uniref50.04` q100. Timed (same binary on / off, 8 cores of node 0 with their siblings
  idle, warmed, order alternated): median ratio of wall time x0.992 (0.977-1.016, 8 pairs) OrthoFinder, x0.999
  (0.980-1.008, 10) Swiss-Prot q100, x0.935 (0.878-1.034, 6) the same with outfmt 0, x1.000 (0.946-1.074, 8)
  `uniref50.04` q100; max RSS 123 / 166 MiB Swiss-Prot q100, 139 / 182 MiB with outfmt 0. Gates: `cargo test` (87
  binaries, 1,993 passed); ci_fast all programs with the switch, with `=0` and with `LOSAT_X_DEFAULT=0` (223 cases, 0
  failures each); blastp `-db` fixtures at 1 and 4 threads under the shadow (52 rows, 0 differing, the exit line
  `identical` in every run that ends without an error).
- BPREDODBSUBJ (graduated 2026-10-27, in the set since 2026-10-22 (perf/bp-redo-db-subjects, evidence at the end of
  this entry); no variable; always on; acts for every task with `-db`) (BLASTP `-db`, memory) the composition redo
  reads the subject of each match from the database when
  it redoes the match, as NCBI does (blast_kappa.c:1381-1385), instead of a copy of every listed subject made before the
  redo (`blastp/x_redo_subjects.rs`, `blastp/blast_engine/database/x_redo_db_subjects.rs`). The previous path,
  `ListedSubjects::read`, copied the residues of every subject of the preliminary hit lists, which the redo and the
  reports read; it is gone (a `cfg(test)` form of `read` remains for the unit tests of the titles).
  `ListedSubjects::x_list` keeps their ordinal ids and lengths (from the index; the estimate of
  `LOSAT_X_BPREDOLPT` reads the lengths), every redo path and the first Newton problems of the Newton lanes read a
  match's stored residues with `SeqDb::read_sequence` and drop them after the match (the Newton problems read ahead of
  the redo and the redo reads again: keeping the reads ahead for the redo raised the maximum resident set; the
  one-thread loop reads after the early termination test, as NCBI, and after the Newton problems; the parallel paths
  read the matches they compute ahead of the test's replay), a read that finds a pyrrolysine (O) sets the subject's
  flag, and after the redo `ListedSubjects::x_read_shown` reads the subjects of the final hit lists again for
  the reports. Where every report is tabular rows of the default fields or Seq-aligns, which read the letters only to
  recount the HSPs whose query or subject has a pyrrolysine (O), it keeps the residues of the subjects with an O (read
  again only where the flag is set: the redo read every final subject) and of the subjects of the HSPs whose query has
  one; otherwise those of every final subject; those of every listed subject when the reports make records
  (`LOSAT_X_BPLISTEDLEAN=0`, ABI v1's titles); the others keep their lengths. The redo reads the same
  stored residues and the reports the same letters and lengths, so reports do not change. A reader error of the redo's
  reads is the error that the copy's read gave, reached in the redo and only for a subject that the redo reads; only a database
  changed after its preliminary read can fail there. The reports are held by the frozen NCBI `-db` fixtures
  (`db.o_stored.*` for stored O) and by `tests/blastp_redo_db_subjects.rs`, which compares the reports of every redo path and
  residue set with the one-thread search. Evidence 2026-10-21 (perf/bp-redo-db-subjects, head 90bb60f8, workstation; the same binary on / off in one job):
  `blastp -db` UniRef50 q100 at 32 threads (node 1, untimed) peak RssAnon 235 / 289 MiB (max RSS 229 / 280 MiB; NCBI
  195 MiB), wall 541 / 543 s, the same output (= NCBI); Swiss-Prot q100 peak RssAnon 91.5 / 121.0 MiB at 8 threads and
  101.4 / 130.1 MiB at 16 (medians of 2 pairs; NCBI 125 and 136 MiB). Instructions (`perf stat` `:u`, 8 threads, on /
  off) x1.0007 OrthoFinder blastp-fast Species4 x BlastDBSpecies6, x1.0002 Swiss-Prot q100 outfmt 6, x1.0003 the same
  with outfmt 0, x1.000 `uniref50.04` q100. Timed (8 cores of node 0 with their siblings idle, warmed, order
  alternated): median ratio of wall time x0.997 (0.991-1.001, 8 pairs) OrthoFinder, x0.995 (0.986-1.083, 10)
  Swiss-Prot q100, x1.002 (0.973-1.066, 8) `uniref50.04` q100; max RSS 94 / 123 MiB Swiss-Prot q100, 121 / 121 MiB
  OrthoFinder, 54 / 54 MiB `uniref50.04` q100 (untimed, 3 rounds: 54.6 / 52.8 MiB; keeping the reads ahead of the
  Newton problems for the redo read 61.9 MiB). Gates: `cargo test` (94 binaries, 2,065 passed); ci_fast all programs
  with the switch, with `=0` and with `LOSAT_X_DEFAULT=0` (223 cases, 0 failures each); blastp `-db` fixtures at 1 and
  4 threads under the shadow (59 rows, 0 differing, the exit line `identical` in every run that ends without an
  error).
- `LOSAT_X_BPREDOPAR=1` (BLASTP, speed; in the strict set since 2026-10-17 (perf/perfmb-default-bundle, evidence at the end of this entry): on unless it is `0`) the composition-based redo of a search with more than one thread and
  at least two matches splits its matches between the threads, as NCBI's `omp for` over its array of
  matches, in jobs of up to four consecutive matches of one query, instead of one job per query (or the
  matches of the only query with matches); the early termination test, the heap insertions and the error
  that ends the redo are then replayed on the calling thread in the order of the one-thread loop, on the
  same heaps (`blastp/x_redo_par.rs`, exactness argument in its header). The jobs run in windows of
  consecutive jobs, one `pool.install` each, that may end inside a query: at most 64 jobs per thread
  (`LOSAT_X_BPREDOPARWINDOW` below), so at most 4 x 64 x threads redone results are held before their
  replay. A job builds its query's workspace once per window (`OnceLock`; a query that spans windows has it
  built again from the same inputs, a small cost that grows with the query's length) and its own
  parameters, prefetches the first Newton problems of its own matches with the Newton lanes (one vector
  of four, against the heaps as the window's start left them), and redoes its matches with per-worker
  scratch and the search's shared SEG store (`SEGSHARE`). A
  query's preliminary list is released after the replay of the window that holds its last match, as the
  one-thread loop releases it (BPPFREE, BPPCOMPACT). The early termination test is
  monotone over the replay (once true for a match's E-value, true afterwards), so a match whose test is
  already true at its window's start is not computed. A match can be ended at all only for a query whose
  batch's heaps can all be filled (with fresh heaps, the last query of a batch whose queries all have at
  least `max_target_seqs` matches) and with a best E-value above 5 x the cutoff; at a window's start, not
  when a query of its batch whose matches are all replayed has a heap that is not filled, nor among its
  query's next `max_target_seqs - n` matches (`n`: the lists its heap holds). A job with a match left open
  is at risk. A window holds at most cap(q) jobs at risk of a query q (it may hold jobs at risk of several
  queries): `threads` at first and after a window that dropped a computed match of q; otherwise cap(q)
  grows by q's jobs at risk that the window computed (doubling when it closed the window), up to the job
  budget. A computed match is dropped only in the window during which its batch's heaps fill, so at most
  one window drops matches of q, and it drops at most 4 x cap(q) of them, cap(q) <= threads + q's jobs at
  risk computed, and all used, before (the one-thread loop computes none); where q's batch's heaps never
  fill, q's windows are back to their budget after a number of windows logarithmic in q's work at risk.
  With
  `LOSAT_TIMING=1` each search the path runs prints
  `[TIMING] x_bpredopar: matches=<n> jobs=<j> windows=<w> split=<s> ended=<e> skipped=<k> wasted=<x>
  par_s=<s> replay_s=<s>` (on one line; `split`: windows that end inside a query; `ended`: the matches the
  early termination test ended, as many as in the one-thread loop; `skipped`: those not computed;
  `wasted`: those computed and dropped). Diagnostics only, the output bytes do not change: the counters
  and traces of the redo count the `wasted` matches as well (`[TIMING] x_bpredofit: lists= hsps=`, the
  kappa range counters, `LOSAT_DEBUG_BLASTP_TRACE_HSP` lines), and `[X_STATS] NEWTONLANES` counts the
  prefetches of the jobs (four matches each, the one-thread loop prefetches 16 ahead). With
  the eight lanes (aarch64) a job's prefetch of at most four problems never fills eight lanes.
  Reason: the query-level redo waits for the largest query of a search (in the blocks of the human x yeast
  blastp-fast job with `LOSAT_X_BPQBLOCK=20000`, 33 % of a block's redo work at the median, 54 % at p90).
  Default ON since 2026-10-17 (perf/perfmb-default-bundle; workstation (2x Xeon Gold 6338, 16 physical cores of one
  node with their siblings idle, alternating rounds, every output the expected one; perf-mb
  `logs/defaultset_ws_report.md` unless named); a reduction is 1 - with / without the switch, the other switches of
  the candidate set on): blastp-fast of 23,879 human proteins against the yeast proteome at 16 threads -32.3 % wall
  with `-subject` (four rounds, -28.4..-39.3 %) and -19.5 % with `-db` (-13.1..-32.7 %), at 4 threads -7.0 / -4.1 %;
  exact blastp `-db` inside the rounds' spread (-7.2 %; rounds at a matched clock +4.3 %); one thread unchanged; max
  RSS -5 MiB.
- `LOSAT_X_BPREDOPARWINDOW=<jobs per thread>` (tuning variable of `LOSAT_X_BPREDOPAR`; not in the strict set,
  read once when the path first runs) the job budget of a window per thread, a positive integer. Unset or
  any other value (0, a negative number, not a number) gives 64, without a message, as the value switches of
  `utils/xswitch.rs` give their built-in value. The tests set 1 to make windows that end inside queries. The
  output does not depend on it.
- `LOSAT_X_BPREDOLPT=1` (BLASTP, speed; in the strict set since 2026-10-17 (perf/perfmb-default-bundle, evidence at the end of this entry): on unless it is `0`; acts only where `LOSAT_X_BPREDOPAR` runs,
  with more than one thread, for a window of two or more jobs) each window of the redo is computed in one
  `rayon::scope` of `min(threads, jobs)` tasks instead of a `par_iter` over its jobs (`blastp/x_redo_par.rs`,
  `x_compute_window_lpt`, exactness argument on it). Each task makes one set of scratch for its first job and
  takes the window's jobs one at a time from one counter, in the order of their estimated cost: the sum over
  the job's matches not ended before the window of query letters times subject letters, largest first, jobs of
  equal cost in index order. Each job's results go to its own slot, and the window's results are read by job
  index as without the switch, so the replay, the matches the early termination test ends and the error that
  ends the redo are those of the `par_iter` path. Without it, rayon's split of the jobs into halves and the
  spread of the last jobs' ends leave threads idle at the end of each window. Measured on the workstation
  (16 cores of one node, HT siblings idle; the human x yeast blastp-fast job at 16 threads with
  `LOSAT_X_BPSUBJRUN=1` and the redo and memory switches, without `LOSAT_X_BPLKPFLAT`; 4 alternating pairs,
  the same outputs): with the switch off the redo windows (`par_s` of the `[TIMING]` line) take 1.56 s of the
  7.84 s with `LOSAT_X_BPQBLOCK=50000` and 2.05 s of 10.60 s with 20,000; with it on, 7.84 -> 7.56 s
  (x0.964) and 10.60 -> 10.17 s (x0.960), the windows' time 1.56 -> 1.28 s and 2.05 -> 1.63 s; max RSS
  101.8-103.0 -> 99.9-102.6 MiB and 82.6-85.1 -> 80.5-81.9 MiB. (The idle time at the windows' ends that
  motivated it was measured before, on a build with `LOSAT_X_BPLKPFLAT=1`, #221, not on main, where the
  job takes 6.9 s / 8.9 s: 0.39 / 0.77 s, of which the largest single job explains 0.01 / 0.10 s;
  perf-mb-20261014 `logs/split2_report.md`.) With `LOSAT_TIMING=1` the `[TIMING] x_bpredopar:` line ends
  with ` lpt_windows=<w>`, the number of windows that ran with the switch (windows of one job and runs at
  one thread have none). Memory: one scratch per task instead of one per rayon work unit.
  Default ON since 2026-10-17 (perf/perfmb-default-bundle; workstation (2x Xeon Gold 6338, 16 physical cores of one
  node with their siblings idle, alternating rounds, every output the expected one; perf-mb
  `logs/defaultset_ws_report.md` unless named); a reduction is 1 - with / without the switch, the other switches of
  the candidate set on): blastp-fast of 23,879 human proteins against the yeast proteome (`-subject`) at 16 threads
  -8.6 % wall (four rounds at a matched clock, -3.8..-10.5 %), at 4 threads -1.3 %; with `-db` inside the rounds'
  spread (rounds at a matched clock -10.2 % blastp-fast, +1.3 % exact); max RSS within 3 MiB. It acts only with
  `LOSAT_X_BPREDOPAR`.
- BPREDOWINEXPAND (graduated 2026-10-27, in the set since 2026-10-22 (perf/bp-redo-window-expand, evidence at the
  end of this entry); no variable; always on; acts only where `LOSAT_X_BPREDOPAR` runs, with the compact lists of
  BPPCOMPACT) (BLASTP, memory) the redo in windows expands the compact preliminary lists as its windows reach them instead of
  the lists of all queries before the redo (`blastp/x_redo_par.rs`, `XExpandAhead`). Before each window the lists from
  the first one not expanded are expanded in query order, by the same `x_expand_hit_list`, until the queries expanded
  from the window's start hold more jobs than the window's budget (or every list is expanded); a query's jobs are its
  matches in jobs of up to four, counted from its compact list as from its expansion. A window holds at most its
  budget and then reads the next job (one query's), so the window and the first match after it lie within expanded
  lists (asserted after each window), and the windows, the replay, the matches ended and the lists released are those of the full
  expansion. The match counts that the early termination bound (`x_end_above`) and the choice of the path read are
  counted from the compact lists. Where the path does not run after all (fewer than two matches), every list is expanded before
  the other paths (the only expansion before the previous switch existed). The expanded lists held at a time are those of the queries of about one
  window, and each query's compact list is freed when it is expanded.
  Evidence 2026-10-20 (v4-D12 s3 probe, workstation): in `blastp -db` UniRef50 q100 at 32 threads the expansion of all
  91 lists of the first query batch (194,347 HSPs of 168 bytes) added 31 MiB to the main heap before the redo, which
  ran in 10 windows. Measured (head 6ef84094, workstation, on / off in one job): that run (node 1, the branch merged
  with `LOSAT_X_BPLISTEDLEAN`'s pull request) peak RssAnon 268.0 / 286.2 MiB (max RSS 258.1 / 276.2 MiB; NCBI 195
  MiB), wall 539.7 / 539.9 s, the same output (= NCBI), the same 10 windows; the end-stage rise 67.9 / 86.6 MB. Less
  than the 31 MiB of the full expansion: with it, the allocations after the redo fitted into the memory that the
  released lists left in the main heap. Swiss-Prot q100 at 8 threads, max RSS (median of 3 pairs): 119.6 / 123.6 MiB
  merged as above, 139.6 / 139.6 MiB with outfmt 0, and 166.4 / 165.4 MiB on the branch alone (without the lean titles
  the ASCII records made after the redo reused the released lists' memory). Instructions (`perf stat` `:u`, 8 threads,
  on / off) x1.000011 OrthoFinder blastp-fast Species4 x BlastDBSpecies6, x0.999987 Swiss-Prot q100 outfmt 6,
  x0.999980 the same with outfmt 0, x0.999998 `uniref50.04` q100. Timed (same binary on / off, 8 cores of node 0 with
  their siblings idle, warmed, order alternated): median ratio of wall time x1.001 (0.958-1.048, 8 pairs)
  `uniref50.04` q100 (max RSS 49 / 53 MiB), x1.000 (0.986-1.010, 10) Swiss-Prot q100, x1.000 (0.995-1.003, 8)
  OrthoFinder. Gates (head 1273374d): `cargo test` (88 binaries, 2,018 passed); ci_fast all programs with the switch,
  with `=0` and with `LOSAT_X_DEFAULT=0` (223 cases, 0 failures each); blastp `-db` fixtures at 1 and 4 threads under
  the shadow (56 rows, 0 differing, the exit line `identical` in every run that ends without an error).
- `LOSAT_X_BPSEGWIN=<entries>` (BLASTP, memory; in the strict set from the pull request that added it (perf/bp-seg-win,
  evidence at the end of this entry), a value switch: an integer of 2 or more, else 4,096; `=0` turns it off; acts only
  where the redo runs in windows, `LOSAT_X_BPREDOPAR` (more than one thread and at least two matches)) after a window of
  the redo is computed, the search's store of SEG-masked subject ranges (`blastp/kappa.rs` `XSharedSubjectRanges`,
  SEGSHARE) is emptied when it holds that many entries or more (`blastp/x_redo_par.rs` `x_bp_seg_win`,
  `XSharedSubjectRanges::x_len` / `x_clear`) instead of holding every range the search masked until the redo ends; a
  smaller store is kept for the next windows, as without the switch. A later window that reaches a range again masks it
  again and computes its composition again (the SEG cache entry keeps it for the redo): an entry is a function of its key (the subject, the
  range and the SEG parameters), so it gets the same bytes, bias flag and composition, and reports do not change (lib
  test `windows_with_the_seg_store_emptied_after_each_equal_windows_with_one_store`: windows computed with the store
  emptied after each and with one store give the same results, every field; `tests/blastp_redo_memory.rs` with the
  default, `=2` and `=0`: the CI fixtures' stores stay below 4,096 entries, so they reach the emptying only with `=2`).
  In `blastp -db` q100 x `uniref50.03` at 8 threads the store held 68,795 ranges (12.6 MiB of letters) at the end of the
  first query block's redo, and emptying it masks 16 % more ranges there (v4-MF REPORT §6.4); with the rule the store
  there is emptied every few windows (the first block's 40 windows hold about 2,000 entries each: 78,501 SEG calls
  against 79,963 with `=2`), and its peak RssAnon does not differ from `=2`'s beyond the runs' spread (4,096 entries
  61.6 MiB, 1,024 entries 62.8, `=2` 63.6; median of 3 runs each, 60.9-64.5). The 177 WSSV proteins x the 86 PajaWSV
  proteins at 8 threads (2 windows) hold 86 entries, which the rule keeps: emptied after every window (`=2`) they free
  no memory (26.6 MiB either way, s12) and the SEG calls go from 89 to 174. Why the rule: emptying that store read x1.024
  wall (x1.012 user) against both switches `=0` on the workstation (v4-D12 s12, node 1, warmed 3 s, 20 rounds, `wait4`;
  in s12b's run x0.997 wall, x1.006 user against the default, which keeps it), so it is skipped where nothing comes
  back. Evidence 2026-10-25 (v4-D12 s12b, workstation, head 0bd4f0a4 on dev 51234c2a, every report identical and the
  same as NCBI's), with `LOSAT_X_BPLKDROP` below; peak RssAnon, untimed, node 0, median of 2 runs (q100) or 3 (WSSV),
  MiB, dev / both switches / both `=0`: q100 x `uniref50.03` at 8 threads 91.9 / **62.2** (5 runs, 60.9-64.0) / 92.5
  (NCBI 150), at 16 threads 105.3 / **95.0** / 106.9 (NCBI 155); WSSV x PajaWSV at 8 threads 26.3 / **24.6** / 26.7
  (NCBI 15.8), at 16 threads 34.1 / **31.0** / 34.1 (NCBI 19.6-20.0). Max RSS (max of `time` %M and VmHWM) q100 x
  `uniref50.03` 100.0 -> 71.4 MiB at 8 threads, 113.4 -> 103.2 at 16. Each switch alone (v4-D12 s12, head 9380a251 on
  dev 31f36ab0, the store emptied after every window, median of 2): q100 at 8 threads 71.8 (`LOSAT_X_BPSEGWIN` only) and
  85.8 (`LOSAT_X_BPLKDROP` only) against dev 92.7, at 16 threads 104.2 and 105.1 against 105.4. One thread does not run
  the redo in windows, and its peak is not in the redo: WSSV 12.5 MiB in every arm (NCBI 10.3), 3,000 human proteins x
  yeast (`-evalue 0.001`) 25.1 -> 24.6 (NCBI 12.2); `perf stat` `:u` (s12, node 0, median of 3, against dev):
  instructions x0.99996 WSSV, x1.00000 human x yeast. Timed (node 1, warmed 3 s, order alternated, against the dev
  binary): WSSV x PajaWSV at 8 threads x0.997 (0.953-1.096, 10/20 below 1; user x1.005; s12b, `wait4`); s12 (`time` %e,
  both switches emptying after every window): human x yeast x0.998 (0.992-1.010, 8 pairs) at 8 threads and x0.999
  (0.986-1.023, 8) at 16, q100 x `uniref50.03` x0.995 (0.993-1.004, 3) at 16, WSSV x PajaWSV x1.000 (0.938-1.000, 16) at
  16. Gates (s12b, 0bd4f0a4): `cargo test --all-features` (11 binaries, 2,534 passed); ci_fast all programs with the
  switches and with both `=0` (223 cases, 0 failures each) and blastp with `LOSAT_X_BPSEGWIN=2` (28 cases, 0 failures);
  blastp `-db` fixtures at 1 and 4 threads (162 rows, 0 differing); `switch_reach_check` 0 failures.
- `LOSAT_X_BPLKDROP=1` (BLASTP, memory; in the strict set from the pull request that added it (perf/bp-seg-win, evidence
  in the `LOSAT_X_BPSEGWIN` entry above): on unless it is `0`; acts for every task, with `-subject` and `-db`) the
  search's lookup table is dropped after the preliminary stage, before the composition-based redo
  (`blastp/blast_engine.rs` `x_lookup_after_prelim`), instead of at the search's end: nothing reads it after the
  preliminary stage. It is kept to the end, as without the switch, only where the search gives its vectors to a spare of
  the block search that a later block's build takes (`LOSAT_X_BPLKPFLAT`, `LOSAT_X_BPLKPAHEAD`,
  `blast_engine/x_lookup_ahead.rs` `XAheadBlock::x_table_taken_later`): in the block search of `LOSAT_X_BPQBLOCK` with
  more than one thread, every block but the last two where the build ahead acts (block i's vectors go to the build of
  block i + 2, which starts when block i + 1's table is ready) and every block but the last where only
  `LOSAT_X_BPLKPFLAT` acts (the next block's own build takes them). So the builds reuse the vectors as before, and only
  the tables that nothing reuses leave early. Reports do not change (`tests/blastp_redo_memory.rs`, which also shows by
  `LOSAT_X_REACH` that a block search of several blocks keeps the tables of its first blocks).
- `LOSAT_X_BPLKREUSE=1` (BLASTP, memory; in the strict set from the pull request that added it (perf/bp-n1-lkreuse):
  on unless it is `0`; exact blastp (`-task blastp`, the table built in place), `-subject` and `-db`)
  in the block search of `LOSAT_X_BPQBLOCK` where no table is built ahead of its block (`LOSAT_X_BPLKPAHEAD` acts with
  more than one thread; so: one thread, Wasm without threads, or `LOSAT_X_BPLKPAHEAD=0`), each block's search gives its
  table's backbone, overflow and pv vectors back to the block search (`x_query_blocks::x_bp_aa_block_spare`) at its
  end, and the next block's build takes them (`lookup/x_lut_direct.rs` `XAaSpares::x_build_on_block_spare`): cleared
  and zero-filled to the new table's sizes, and a vector short of its size grows by a reallocation of its own
  allocation (`reserve_exact`), never by a new allocation and a free. Without the switch each block allocated its table
  and `LOSAT_X_BPLKDROP` freed it before the redo; glibc raises its dynamic mmap threshold to the size of a freed mapped
  allocation, so the next block's table (up to that size) went to the main heap, and the heap's free space stayed
  resident below what was allocated after it (`GLIBC_TUNABLES=glibc.malloc.mmap_threshold=131072`, which keeps the
  threshold fixed, gave about the same peak as the switch: 17.5 MiB below). A reallocation does not move the threshold
  (glibc grows a mapped allocation with `mremap`). With the switch the table of every block but the last is kept until
  its search's end, as where the build ahead takes it (`LOSAT_X_BPLKDROP`); the last block's is dropped before the redo.
  Reports do not change: the vectors taken are cleared and zero-filled to the new table's sizes before the build writes
  them, as a new table is (lib test `aa_block_spare_builds_each_block_on_the_vectors_given_back` compares each table
  with a fresh build; `tests/blastp_redo_memory.rs` compares reports with `=0` at 1 and 2 threads). Evidence 2026-10-26
  (v4-D12 s14, workstation, node 0, head 591e4fc4 on dev 70d75af5, every report identical and the same as NCBI's):
  peak RssAnon at one thread, MiB, switch / `=0`, 2 runs each (Swiss-Prot 1): 3,000 human proteins x yeast
  (`-evalue 0.001`, 23 blocks) **17.6, 18.6** / 25.1, 25.1 (NCBI 12.2; `=0` with the threshold held fixed 17.5); WSSV x
  PajaWSV (2 blocks) 12.5, 12.5 / 12.5, 12.5 (NCBI 10.3); Swiss-Prot q100 (2 blocks) 73.6 / 73.4. `perf stat` `:u`
  against `=0`, rounds rotated with an A/A arm (`=0` again), load 8-39: human x yeast 12 rounds instructions x1.000006
  (the zero fill), cycles x1.000 (0.961-1.060) where the A/A arm read x1.001 (0.939-1.034), x0.999-1.002 in the 4
  quietest rounds (load 8-16); WSSV 8 rounds instructions x0.999998, cycles x1.002 (0.986-1.008), A/A x0.999
  (0.985-1.010); minor faults 13.2k -> 5.2k (human x yeast).
- `LOSAT_X_BPIHPOOL=1` (BLASTP, memory, every task, `-subject` and `-db`; in the strict set from the pull request that
  added it (perf/bp-ih-pool, v4-D12 s21-s22, `DECISIONS` 62-63): on unless it is `0`) a subject's array of initial hits
  (`BlastpSubjectScratch::init_hsps`, the ungapped hits `blastp_save_init_hsp` saves) past NCBI's first array size lives
  in one pool shared by the threads (`blastp/blast_engine/x_ih_pool.rs`). NCBI starts each thread's list at
  `MIN_INIT_HITLIST_SIZE` (100) and doubles it with `realloc` (`blast_extend.c:221-224`, `336-341`); the list is reset,
  never shrunk, for each subject, so each thread keeps the largest array any of its subjects needed. The scratch did
  the same, and each worker's malloc arena also kept the arrays its vector outgrew, so the preliminary stage held the
  sum over the threads of each thread's largest (UniRef50 q100 at 32 threads: 53 MiB of capacity, of which a subject
  used at most 0.2 MiB at any time). With the switch on the scratch holds an array of 100; a subject that fills it
  continues in an array of twice the length from the pool (the smallest that holds it, or a new one) and the scratch's
  array goes to the pool; a pool array that fills is exchanged for a larger one from the pool, or doubled in place when
  the pool has none; at the subject's end the array goes back to the pool and the scratch takes an array of 100 from
  it. The hits, their order and every value are those of the scratch's own vector (lib tests
  `hits_through_the_pool_equal_the_scratch_vector`, `threads_sharing_the_pool_keep_their_own_hits`;
  `tests/blastp_ih_pool.rs` compares reports with `=0`: `-subject` and `-db`, 1 and 4 threads, a subject past 100,
  200 and 400 hits). At one thread the pool holds the largest array, as the scratch did. `LOSAT_TIMING=1` prints
  `[TIMING] x_bpihpool: pooled=<n>` per search (the arrays that continued through the pool). After a search's
  preliminary stage the pool's arrays are freed (`x_ih_pool::release`, next to the lookup table's drop), as the
  subject loop's scratches freed theirs; without that, Swiss-Prot q100 at one thread, whose peak comes after the
  preliminary stage, kept the largest array. The comparison before the push changes the code of the fused two-hit
  loop: `=0` runs the previous allocation path, not the dev base's instructions (x0.9954-0.9959 of them, as the default).
  Evidence 2026-10-28 (v4-D12 s22, workstation, heads 4bab83fa and fece4d68 (with `release`) on dev 4403cc12; arms default,
  `=0`, and the dev base built the same way; every report identical): peak RssAnon of UniRef50 q100 `-evalue 0.001`
  at 32 threads (node 1, 50 ms samples, runs alternated) 158.5, 160.0, 161.4, 162.1 MiB against `=0` 183.1, 185.1,
  189.3 MiB (wall 402.9-403.4 against 403.3-408.3 s); q100 against a third of UniRef50 at 16 threads (4 pairs side
  by side) 63.9-72.7 against 66.0-74.5 MiB; one thread (node 0, 2 ms samples) WSSV x PajaWSV `-db` 9,192-9,240
  against the dev base 9,188-9,228 KiB, `-subject` 9,004-9,008 against 9,008-9,012, 3,000 human x
  `BlastDBSpecies6` 13,240-13,264 against 13,264-13,300, Swiss-Prot q100 76,600-78,012 against 75,340-76,720 and
  `=0` 75,628-77,748 (its peak is the last 0.4 s, the output stage, and spreads as much in every arm). Time: one
  thread (`perf stat` `:u`, node 0, arms rotated) instructions default / `=0` x0.99998-1.00000, default / dev base
  x0.9953 (WSSV), x0.9959 (human), x0.9956 (Swiss-Prot); cycles default / dev base x0.996 and x1.003 (WSSV, 4 and 8
  rounds), x0.997 and x0.995 (human, 3 and 6), x0.995 and x0.995 (Swiss-Prot, 2 and 3); 16 threads (node 1 CPUs
  32-47 warmed, 8 rounds rotated with an A/A arm, head 4bab83fa) Swiss-Prot q100 wall x1.0000 (A/A x0.9992), cycles
  x1.0008 (A/A x0.9986), `bp_hsap` wall x1.0000 (A/A x1.0000), cycles x1.0004 (A/A x1.0017).
- `LOSAT_X_BPHITSFLAT=1` (BLASTP, memory, every task, `-subject` and `-db`; in the strict set from the pull request that
  added it (perf/bp-hits-flat, v4-D12 s26, `DECISIONS` 71-72): on unless it is `0`) a subject's preliminary HSPs go to
  their queries' compact lists through one array of their positions in query order (`blastp/blast_engine/x_hits_flat.rs`)
  instead of one vector per query kept in the subject scratch (`BlastpSubjectScratch::preliminary_hits_by_query`, cleared,
  never shrunk). NCBI makes the per-query lists for each subject (`calloc` of `num_queries` pointers and
  `Blast_HSPListNew` per query with HSPs, `hspfilter_collector.c:112-153`) and hands them to the hit lists, so a thread keeps
  nothing of them; the scratch kept, per query, the most HSPs that query had in any one subject of the thread, so the
  preliminary stage held, per thread, the sum over the queries of each query's largest list (UniRef50 q100 at 32 threads:
  14.0 MiB of capacity). With the switch a stable counting sort
  (`order_by_query`: counts, ends, filled from the last HSP back) puts the positions in query order, and each query's HSPs are
  read in place from the subject's list in the order of that list, the order the per-query vectors held them in
  (`x_build_compact_preliminary_hsp_list` takes an iterator); a thread keeps 4 bytes per HSP of its largest subject and 4 per
  query. The lists, their order and every value are unchanged (lib test `x_hits_flat_order_matches_the_per_query_vectors`:
  200 random subjects, 1-5 queries, a query index past `num_queries`).
  Evidence 2026-10-30 (v4-D12 s26, workstation, head e03234b2 on dev a46c6897c; arms default, `=0`, and the dev base
  built the same way, 8 rounds rotated, peak RssAnon from 2 ms samples; every report identical across arms): UniRef50
  q100 `-evalue 0.001` at 32 threads (node 1, membind 1) 129.8 [126.9-134.7] MiB against the dev base 143.1
  [141.7-144.2] and `=0` 142.5 [139.3-144.7] (9.4-15.4 MiB lower in 8 of 8 rounds; wall x0.9982 [0.9905-1.0066]); q100
  against a third of UniRef50 at 16 threads (two runs side by side, 12 slots): the highest sample of the preliminary plateau (30-110 s) 39.9
  [39.7-40.1] against 42.0 [41.8-42.1] and `=0` 42.0 [41.8-42.1], while the peak, set by the jump after that plateau
  (about 146-151 s), spreads over 64-81 MiB in every arm (69.9 against 69.1 and `=0` 66.9);
  3,000 human x `BlastDBSpecies6` at 16 threads 37.4 [37.1-37.5] against 37.4 [37.2-37.6]; one thread (node 0) human
  12.9 against 13.0, OrthoFinder Blast0_6 19.2 against 19.2 MiB. Census (probe exp/d12s26-treeprobe 5abfcfc3, UniRef50
  q100 at 32 threads): the per-query vectors held 14.0 MiB of capacity.
  Time: one thread (`perf stat` `:u`, node 0) instructions x1.0003 of the dev base in both cases (`=0` x1.0000), cycles human x0.9985 (`=0` x0.9970), Blast0_6 x0.9971
  (`=0` x0.9984); 16 threads (node 1 CPUs 32-47 warmed 3 s, 8 rounds) Swiss-Prot q100 instructions x1.0005, cycles
  x0.9980 (8 of 8 below 1; `=0` x0.9963), wall x0.9966; human instructions x1.0003, cycles x0.9972 (`=0` x0.9966), wall
  x0.9956.
- `LOSAT_X_BPTREECELLS=1` (BLASTP, memory, every task, `-subject` and `-db`; in the strict set from the pull request that
  added it (perf/bp-tree-cells, v4-D12 s28): on unless it is `0`) the grid index of blastp's preliminary interval tree
  takes at most 1 << 14 cells per subject (`X_BP_TREE_MAX_CELLS` in `blastp/blast_engine.rs`, through
  `BlastIntervalTree::reset_with_bounds_cells`) instead of the tree's own limit of 1 << 16 (`X_MAX_CELLS` in
  `blastn/interval_tree.rs`; blastn keeps it through the unchanged `reset_with_bounds`). The grid is LOSAT's side index
  of the tree (ITREEFAST); NCBI walks the tree (`BlastIntervalTreeContainsHSP`, `blast_itree.c`) and has no grid. The
  subject scratch keeps the tree, and the grid's per-cell heads (4 bytes each) stay as large as the largest grid the
  thread has needed, so the limit bounds them at 64 KiB per thread instead of 256 KiB. A grid that would exceed the
  limit gets cells twice as wide in both offsets, as often as needed; the cell of a start point still lists every leaf
  whose box holds that point, and the unchanged predicate tests them, so the index answers as before (lib test
  `x_side_indexes_agree_with_tree_walk_on_random_hsps` runs the grid at 1 << 16, 1 << 14 and 64 cells against the walk).
  Subjects whose grid had fewer than 1 << 14 cells keep the same grid (probe exp/d12s27-tree-hist, 8 threads, subjects
  with HSPs: human x `BlastDBSpecies6` 89.6 %, Swiss-Prot q100 98.0 %, UniRef50 q100 `-evalue 0.001` 87.95 %, q100
  against a third of UniRef50 100 %). Every report is unchanged.
  Evidence 2026-11-01 (v4-D12 s28, workstation, head b20c3448 on dev b94e9844f; arms default, `=0`, and the dev base
  built the same way, rotated, peak RssAnon from 2 ms samples; every report identical across arms): UniRef50 q100
  `-evalue 0.001` at 32 threads (node 1, membind 1, 8 rounds) 118.5 [116.2-119.9] MiB against the dev base 130.8
  [128.9-133.7] and `=0` 130.2 [128.6-133.2] (11.4-14.4 MiB lower in 8 of 8 rounds; wall x0.9973 [0.9926-1.0051]); q100
  against a third of UniRef50 at 16 threads (two runs side by side, 12 slots, no subject with HSPs above the limit in the
  probe): the highest sample of the preliminary plateau (30-110 s) 39.9 [39.8-40.1] against 39.9 [39.9-40.0] and `=0` 40.0
  [39.9-40.2], and the peak 69.6 [63.1-73.9] against 71.2 [63.9-86.0] and `=0` 69.5 [60.4-78.5], so no effect shows; 3,000 human x
  `BlastDBSpecies6` at 16 threads (8 rounds) 35.5 [35.3-35.6] against 37.7 [37.2-38.3] and `=0` 37.5 [37.2-37.9]
  (1.7-3.0 MiB lower in 8 of 8); one thread (node 0, 8 rounds) human 12.9 against 13.0, OrthoFinder Blast0_6 19.1
  against 19.3 MiB (8 of 8 lower in both).
  Time: one thread (`perf stat` `:u`, node 0, 8 rounds) instructions human x0.99997 of the dev base (`=0` x1.00000),
  Blast0_6 x0.99999 (`=0` x1.00001), cycles human x0.9943 (`=0` x0.9941), Blast0_6 x0.9985 (`=0` x0.9982); 16 threads
  (node 1 CPUs 32-47 warmed 3 s, 8 rounds) Swiss-Prot q100 instructions x1.0000, cycles x0.9961 (`=0` x0.9964), wall
  x0.9949; human instructions x0.99997, cycles x0.9966 (`=0` x0.9962), wall x0.9955. blastn's machine code does not
  change: the 729 functions of the two binaries named for blastn or the tree disassemble alike with addresses
  normalized (only a blastp closure and 3 `OnceLock` helpers differ); megablast Sakai x MG1655 at one thread: instructions
  within 1,200 of 4.52 G (the run-to-run spread).
- `LOSAT_X_BPLK16=1` (BLASTP, memory; in the strict set from the pull request that added it (perf/bp-lk16, v4-D12 s19-s20,
  `DECISIONS` 59-61): on unless it is `0`; recorded for `LOSAT_X_REACH` once per table, taken where the table has 16-bit
  hits): exact blastp's block tables (`-task blastp`, the table built
  in place) hold 16-bit hits in 8-byte cells where the block's lookup query is at most 65,535 letters
  (`BlastAaLookupTable<BackboneCell16>`, `lookup/backbone.rs`: `num_used: u16` and `entries: [u16; 3]`; a cell with more
  than three hits keeps its cursor in the overflow array as a `u32` over `entries[0]` (low) and `entries[1]` (high)).
  A block's lookup query is at most about 50,000 letters (`LOSAT_X_BPQBLOCK`), so every query offset and every cell's
  hit count (at most one hit per offset) fits 16 bits; the overflow array and the cells halve (one thread, ws:
  `bp_wssv` overflow 6.20 MiB, `bp_hsap` 8.56 MiB with 32-bit hits). The table holds exactly the hits of the table of
  32-bit hits in the same order, and the fused two-hit loop (`blastp/x_seed_fused.rs`, generic over the hit type
  `XHit`: `i32`, `u16`) reads them with one `vpmovzxwd ymm, m128` for a batch of eight and the same widening load for
  the tail of fewer than eight hits (the lanes past the hits are blended away as before). Every table of every other
  program, the `-window_size 0` search (its offset-pair scan reads the 32-bit table), the query chunks of a split
  batch, `-task blastp-fast`, the reference builder (`LOSAT_DIAGNOSTICS`) and a block above 65,535 letters keep the
  32-bit table. The 16-byte tail load of a `u16` hit list may read up to 14 bytes past the list; the table keeps two zero
  cells after its backbone and eight zero hits after its overflow array so that the load stays inside the allocation
  (the padding is in no count, `longest_chain` or presence bit). The spares of the block search keep the vectors of one
  width (`XAaSpares::give16`, `x_make_room` with the width); a spare of the other width is dropped, not reused.
  `LOSAT_X_BPLK16SHADOW=1` (diagnostic, never in the set): for each table of 16-bit hits the reference builder also
  builds the table of 32-bit hits and the two are compared cell by cell (same `num_used`, same hit list, same pv,
  `longest_chain` and fields; the padding is zero), panicking with the cell's index on a difference. Lib tests
  `table_of_u16_hits_equals_the_table_of_i32_hits` (short queries, a low-complexity query, the block at 65,535 and
  65,536 letters), `spares_serve_the_builds_of_both_widths` and
  `scans_of_the_table_of_u16_hits_equal_the_scans_of_the_table_of_i32_hits` (every cell test and scan of the fused
  loop). The tables and loops of 32-bit hits are the generic instances for `BackboneCell`/`i32`, whose reads are the
  expressions the code had before; with the switch off, every program runs the instructions of the code before
  (`perf stat` instructions:u of the dev base and of `=0` x1.0000 on three blastp cases and a blastx one).
  Evidence (ws, binary 4bcbaf70 on dev 03528945 against `=0` and the dev base, reports identical; v4-D12 s20): peak
  RssAnon at one thread (node 0, samples every 2 ms) `bp_wssv` (WSSV x PajaWSV) `-subject` 8,976-8,980 against
  12,376-12,380 KiB, `-db` 9,208-9,212 against 12,608-12,612, `bp_hsap` (3,000 human x `BlastDBSpecies6`) 13,240-13,244
  against 17,852-17,860, Swiss-Prot q100 75,316 against 75,288 (its peak is the composition stage); UniRef50 q100
  `-evalue 0.001` at 32 threads (node 1) 181.8 against 191.9 MiB (wall 406.8 against 415.3 s); q100 against a third of
  UniRef50 at 16 threads (six pairs side by side): the search stage 42.8-43.1 against 47.4-47.8 MiB, the peak (the
  traceback stage, after the table is freed) 63.2-82.1 (median 70.9) against 67.1-72.6 (median 69.4), no difference
  in rank (sum 41, 39 expected), wider spread. Time: one thread (`perf stat` :u, node 0, A/A arm) `bp_hsap` cycles
  x0.973 (8 of 8 below 1; A/A x1.008), wall x0.980; Swiss-Prot q100 cycles x0.981 (4 of 4; A/A x0.996), wall x0.979;
  `bp_wssv` cycles x0.994 (8 of 8; A/A x1.003); instructions x0.9992, x0.9983, x1.0003; 16 threads (node 1 CPUs
  32-47 warmed, 16 rounds rotated with an A/A arm) Swiss-Prot q100 wall x0.991 (15 of 16; A/A x0.999), cycles x0.990
  (16 of 16), `bp_hsap` wall x0.980 (16 of 16; A/A x1.000), cycles x0.985.
- `LOSAT_X_BPDBWIDE=1` (BLASTP, every task, speed; in the strict set from the pull request that added it
  (perf/bp-db-wide, evidence at the end of this entry): on unless it is `0`; acts only where `LOSAT_X_BPDBKEEP`'s window
  runs: `-db`, more than one thread, and a database of more blocks than 8 per thread) the window of `LOSAT_X_BPDBKEEP` is 64 blocks per thread (`X_BPDBWIDE_BLOCKS_PER_THREAD`), and past the old 8
  per thread a block is taken only while at most 1 MiB per thread of results wait for the merge (`XWindow::full`,
  `blastp/blast_engine/database/x_db_keep.rs`). The old window was shorter than the slow blocks at the start of a pass:
  UniRef50 is sorted longest first (up to 47,538 residues), and with blastp-fast a block of the first tenth took 39 ms at p99 and up to 1.2 s,
  against about 38 ms of the wall clock for 64 blocks at 8 threads and a median of 4.8 ms per block. Three rounds on 16
  idle cores of node 0 (`dropin-survey-20261029/BIGDB/s3/bpab/ab.tsv`, the same output in every arm), on / `=0` / dev:
  10 queries blastp-fast x UniRef50 at 16 threads 5.46 / 7.49 / 7.40 s, max RSS 45-47 / 42-43 / 41-44 MiB; 100 queries 48.4
  / 57.8 / 56.3 s, 123-126 / 113-116 / 113-117 MiB; at 8 threads 9.45 / 11.1 / 11.1 s and 88.5 / 95.5 / 95.5 s. Default
  blastp, where the waits are few, within the spread: 10 queries 73.7 / 74.9 / 71.1 s. A database of 7 proteomes (69 MB)
  blastp-fast 15.6 / 16.8 / 16.8 s, 164-166 / 155-169 / 165-172 MiB. NCBI: 14.4, 177.9 and 57.9 s for the three blastp-fast
  cases at 16 threads. One thread does not run the window.
- `LOSAT_X_BPDBKEEP=1` (BLASTP, memory; in the strict set from the pull request that added it (perf/bp-db-thresh):
  on unless it is `0`; the preliminary stage of every BLASTP task with `-db`, at any thread count; also the chains of
  every BLASTP lookup build, `-subject` too, and the reads of the index check and the title checks wherever a program
  opens a database) keeps glibc's dynamic mmap threshold at its first value (128 KiB) through the preliminary stage
  of a query block. glibc maps an allocation of 128 KiB or more and raises the threshold to the size of each mapped
  allocation that is freed (up to 32 MiB); after that, allocations below the threshold come from the thread arenas,
  and the space they leave there stays resident. Without the switch every task of a round made a search scratch
  (offset pairs 256 KiB, a diagonal table of up to 512 KiB) and read its block into a new buffer of up to 1 MiB, and
  freed both; before the stage the index check read 1 MiB at a time, the title checks read 1 MiB pieces and made a
  list of all ranges (397 KiB with UniRef50), and the lookup build freed its chains (448 KiB): each free raised the
  threshold. With the switch:
  (1) each worker of the pool keeps one scratch and one block buffer for the whole stage
  (`blastp/blast_engine/database/x_db_keep.rs`), and the blocks are cut to at most 96 KiB (`X_KEEP_BLOCK_BYTES`, at
  one thread too), so that a kept buffer stays in its worker's arena instead of a mapping of its own beside it;
  (2) the workers take the blocks one at a time in ordinal order and put each block's results into a window of 8
  blocks per worker; the worker that ends the oldest block not yet merged merges it and the later blocks already
  ended, in ordinal order, and a worker takes block j only once the blocks before j - W are merged (W, the window).
  No round ends, where every worker waited for the round's last block and its merge: with 96 KiB blocks the rounds
  were 11 times as many (Swiss-Prot at 16 threads 2 -> 17, UniRef50 at 32 threads about 74 -> 794);
  (3) the index check (`seqdb/index.rs`) and the header pieces of the title checks read 64 KiB at a time, and
  BLASTP's title checks read each range from the index (`SeqDb::x_header_range`) instead of a list of all ranges;
  (4) the cells and chains of BLASTP's lookup build (`x_bp_build_direct`) come from allocations above 32 MiB
  (`utils/x_mapped.rs` `x_unthresholded_capacity`), which glibc maps whatever the threshold and whose free does not
  move it (the pages not written are not resident); the chains' slabs are cut from one such region while it has room
  (`tblastx/lookup/backbone.rs`), so a build reserves one region of address space per 32 MiB of chains, not one per
  256 KiB slab (64-bit Linux with glibc only; elsewhere the sizes are the plain ones).
  Reports do not change: each block is searched as a task of the rounds searches it, from the one-thread order's
  diagonal offset of its first sequence, and the blocks' results go to the hit lists in ordinal order; a failed
  search stops at the first error in ordinal order, and `x_db_cut.rs` returns the one-thread order's error. The
  waits end: the oldest block's search waits for nothing a waiting worker holds and runs no rayon work, and a panic
  stops the window (`x_db_keep.rs` header). Lib tests `kept_scratches_and_buffers_give_the_search_without_them`,
  `the_window_waits_for_its_oldest_block_and_ends_on_a_panic` and
  `a_damaged_database_gives_the_one_thread_error_with_the_cut` compare the sequences, offsets, hit lists and errors
  with the rounds without the switch at 1-7 threads. Evidence 2026-10-27 (v4-D12 s18, workstation, node 1 with
  membind 1 at 16 and 32 threads, node 0 at one thread, binary f465af2d on dev f5326a59, every report identical and
  the same as NCBI's), switch / `=0`: peak RssAnon UniRef50 q100 `-evalue 0.001` `-db` at 32 threads **188.1 /
  212.1 MiB**, wall 411 / 425 s; `-db uniref50.03` at 16 threads **59.0, 59.2 / 95.7, 93.0 MiB**, wall 161.3, 161.6 /
  169.1, 169.4 s; WSSV x PajaWSV `-db` at one thread **12,452, 12,528 / 12,784, 12,780 KiB**; 3,000 human proteins x
  `BlastDBSpecies6` at one thread 17,792, 17,824 / 19,076, 18,068 KiB; Swiss-Prot q100 at one thread 75,040 /
  75,232 KiB. 16 threads, 16 rounds rotated with an A/A arm (`=0` again), CPUs warmed 3 s, load 16-47: Swiss-Prot
  q100 wall x0.943 (16 of 16 below 1; A/A x1.000), cycles:u x1.003, instructions x1.0001, minor faults x0.83;
  human x `BlastDBSpecies6` wall x0.979 (15 of 16; A/A x1.003), cycles:u x0.999. One thread `perf stat` `:u` (8 and
  4 rounds): instructions x0.9999 / x1.0000, cycles x0.999 / x1.003 where the A/A arm read x1.001 / x1.005. Forms
  measured on the way (v4-D12 s16-s18, same cases): kits with 1 MiB blocks, UniRef50 n32 234.9 against 210.6 MiB (a
  kept 1 MiB buffer is a mapping beside the arena); the fixes before the stage without kits, 210.2 against 209.0 (the
  first task's block free raises the threshold to 1 MiB); kits with 96 KiB blocks in rounds, 184.0-186.3 against
  206.6-209.2 MiB but Swiss-Prot n16 cycles:u x1.002-1.004 and wall x1.000-1.004; rounds of the cut's bytes,
  199.7 MiB and wall x0.949. The window holds about 3 MiB more than those rounds with UniRef50 at 32 threads (a
  worker, not the main thread, merges). A threshold held by `GLIBC_TUNABLES` (mmap and trim at 128 KiB, a diagnostic)
  gave 136 MiB with 144 times the minor faults.
- BPSEED (graduated 2026-10-26, in the set since 2026-10-14; no variable, no `[X_SHADOW] BPSEED` line; always on;
  `LOSAT_X_BPSEED_PREFETCH` became the constant 16) (BLASTP, speed) the two-hit seed
  stage with an amino-acid lookup table
  (`blastp/x_seed_fused.rs`, exactness argument in its header): the subject scan and the diagonal test in one
  loop that reads each lookup cell's list in place, in the reference's order, and hands only the hits that need
  the context check and the ungapped extension to the reference's code; the context of a candidate comes from a
  two-level table (a coarse table of at most 16,384 entries and a walk over the frame bases) instead of the
  binary search; the diagonal entry of a later hit of the cell is prefetched 16 hits
  ahead. Every two-hit blastp search runs it (`-window_size 0` keeps the offset-pair buffer of
  `blastp_for_each_offset_pair`, NCBI's one-hit word finder). Measured 2026-10-13 on the workstation
  (OrthoFinder-style job,
  fly x human, every 10th protein, one thread, medians of 6 rotating rounds under load): 51.4 s -> 34.2 s,
  output unchanged. The loop also covers the compressed lookup table of blastp-fast
  (`blastp/x_seed_fused_compressed.rs`, exactness argument in its header): the scan (NCBI's
  `s_BlastCompressedAaScanSubject` statements) lists the subject offsets whose word passes the PV test, then
  reads those cells in scan order, the backbone cell 8 cells ahead prefetched, and hands each cell's hits in
  NCBI's order (the backbone's offsets, then the overflow list from the cell filled last) to the same test and
  closure; the offset-pair array's rounds are not repeated (a resumed round reads the index the rolling scan has,
  as NCBI's reciprocal division is exact below the table size). Measured 2026-10-14 on the WSL machine (blastp-fast,
  2M query letters x 1M subject letters, `LOSAT_X_BPQBLOCK=20000` on both sides, one thread, 8 alternating pairs,
  load 7.5-10.9): `blastp_scan_ungapped` summed over the blocks 2.05 -> 1.40 s (median ratio 0.67), wall
  8.94 -> 8.22 s (0.89), output unchanged. The unit tests of `x_seed_fused.rs` compare the cell tests with the
  per-hit test, and the frozen NCBI fixtures hold the reports.
- `LOSAT_X_BPCLUT=1` (BLASTP `-task blastp-fast`, speed and memory; in the strict set since 2026-10-17, on by
  default together with `LOSAT_X_BPLKPFLAT`) the compressed lookup table (`CompressedAaLookupTable`) is
  built flat in three passes (`tblastx/lookup/x_clut_direct.rs`, exactness argument in its header): a count pass
  makes the reference builder's appends in its order (the words of `add_word_hits`, the neighbors' last letter
  taken as a prefix of the sorted row) and records the cell of each append; a layout pass places each cell's hits
  in one slice of one array (an array of starts, and the reference's PV array and longest chain from the counts);
  a fill pass writes the hits in the order of their appends and then moves the hits of each cell of more than 5
  to NCBI's scan order (`x_emission_position`: the first three, the overflow cell filled last, then the older
  cells newest first). The test for NCBI's overflow banks (the per-batch refusal) is the reference's, with the
  same message. Each block of `LOSAT_X_BPQBLOCK` builds its own flat table. The scans read a cell as its slice:
  the fused loop (`blastp/x_seed_fused_clut.rs`) (`-window_size 0`: the offset-pair scan,
  `scan_subject`'s statements). The table holds 4 bytes per hit plus 4 per backbone cell (the
  reference: 32 bytes per backbone cell plus 8 per overflowed hit), and its build holds another 4 bytes per hit.
  Measured 2026-10-14 on the WSL machine (blastp-fast, 2M query letters x 1M subject letters,
  `LOSAT_X_BPQBLOCK=20000` on both sides, one thread, 8 alternating pairs, load 5.9-8.2):
  `blastp_build_lookup` summed over the blocks 4.55 -> 2.27 s (median ratio 0.54), `blastp_scan_ungapped`
  1.41 -> 1.08 s (0.84), wall 8.21 -> 5.43 s (0.73), max RSS 77 -> 40 MB, output unchanged.
  With more than one thread (`-num_threads`) the search's thread pool builds the same table
  (`tblastx/lookup/x_clut_par.rs`, exactness argument in its header): the query offsets of each count range
  are cut into chunks in order, each chunk makes its appends (the serial count pass's statements) and sorts
  them stably by stripe (2^shift backbone cells, at least the largest PV shift); each stripe then reads its
  part of every chunk in chunk order, which gives each cell its offsets in the serial order, counts them,
  lays out its cells' starts from its first hit (the hits of the stripes before) and sets its PV words, and
  runs the serial fill and the move to NCBI's scan order on its cells. The overflow cells that each of
  NCBI's tables takes are counted per stripe and summed, and each table is tested at its end by the serial
  builder's rule (same refusal, same message). At one thread the serial builder runs, unchanged. The build
  holds 4 bytes per hit for the sorted appends, as the serial build's record of the appends.
  Measured 2026-10-14 on the WSL machine (i9-14900HX, 8 performance and 16 efficiency cores, shared with
  other sessions; same input, `LOSAT_X_BPQBLOCK=20000 LOSAT_X_BPCLUT=1` on both sides,
  with this change and the per-block set-up on the pool of `LOSAT_X_BPQBLOCK` below; alternating pairs,
  medians; load 15-32): wall at 1 / 4 / 8 / 16 threads 7.36 / 5.00 / 4.55 / 7.92 -> 7.19 / 2.74 / 2.33 /
  3.97 s (median pair ratios 0.98 / 0.55 / 0.51 / 0.49), `blastp_build_lookup` 3.04 / 3.07 / 3.00 / 4.91
  -> 3.03 / 1.03 / 0.86 / 1.41 s, max RSS 39 / 45 / 46 / 53 -> 42 / 42 / 45 / 63 MB, output unchanged.
  Measured 2026-10-14 on the workstation (Xeon Gold 6338, node 1, load about 40 from other jobs; human x
  yeast proteomes, `-subject`, `LOSAT_X_BPQBLOCK=20000`; binaries of that day, not the
  current ones): the table built serially (#181) -> on the pool (4d80432f) took wall 96.0 -> 92.5 s at 1
  thread, 62.7 -> 33.7 at 4, 57.5 -> 31.5 at 8 and 55.4 -> 26.2 at 16, max RSS 68 -> 71, 66 -> 72, 67 -> 78
  and 73 -> 83 MiB, every output the bytes of NCBI BLAST+ 2.17.0 `-subject`; the build alone 28.3 -> 27.0 s
  at 1 thread and 31.2 -> 5.9 s at 16. At one thread without and with the flat table (without -> with
  `BPCLUT`; four jobs run at once on four cores): human x yeast 110.9 -> 94.3 s and 103 -> 68 MB,
  Arabidopsis x zebrafish 366.9 -> 352.7 s and 96 -> 88 MB.
- `LOSAT_X_BPCLUTSCANBF=1` (BLASTP `-task blastp-fast`, speed; in the strict set from the pull request that
  added it, 2026-10-23, line v4-G #11; `=0` turns it off; acts in the fused loop on the flat
  compressed table of `LOSAT_X_BPCLUT`, on tables with at least `LOSAT_X_BPCLUTSCANBFMIN` used cells per 10,000;
  the reference compressed table keeps its scan) the flat table's scan (`tblastx/lookup/x_clut_direct.rs`
  `for_each_cell::<true, _>`) writes every subject offset that reaches `PV_TEST` at the end of its list and
  advances the list's length by the presence bit, the form of the standard table's `x_scan_cells_bf` (BPSCANBF). The
  offsets, their indices and their order are those of the loop with the branch, so the cells visited and the
  report are the same (the lib test `flat_scans_give_the_reference_scans_hits` compares both lists with the
  reference scan, `the_flat_loop_is_the_fused_loop_on_the_reference_table` both instances of the fused loop with
  the reference table's loop, and `tests/blastp_fast_scan_bf.rs` the reports with `=0`'s). The table records its
  cells with hits (`used_cells`, `finalize`'s count, set by the serial and the parallel builder), and the
  branch-free list is taken for a table with at least `LOSAT_X_BPCLUTSCANBFMIN` used cells per 10,000 (built-in
  500 natively and 0 on wasm32; `0` every table): natively, with fewer, few offsets pass `PV_TEST`, the branch is
  predicted, and the store at every offset costs more than the mispredicts it saves (one core on the workstation,
  branch-free / branch cycles: 1.071 at 2.3 % of the cells used, 1.001 at 4.4 %, 0.956 at 7.3 %, 0.845 at 14.9 %,
  0.800 at 22.6 %); in WebAssembly the loop with the branch is the slower one at every density measured (the
  serial WASI command in Node 26, one instance per search in one process, 3 rounds, workstation node 0: branch-free /
  branch time 0.81 at 1.2 %, 0.80 at 2.3 %, 0.70 at 7.3 %, 0.81 at 31.6 %; the 357-residue protein 13.33 -> 10.62 s,
  10 proteins 51.87 -> 42.10 s). The decision is
  recorded per query block (`LOSAT_X_REACH`). The branch-free instance (`blastp/x_seed_fused_clut.rs`
  `x_scan_flat_bf`) is out of line with a closure of its own type, so the instance with the branch is compiled as
  without it; the forms that put the scan in a function of its own, or both lists inline before one visit loop,
  ran 0.5-10 % more cycles with the branch. Timed 2026-10-23 on the workstation (node 0 CPUs 8-11; dev cd8e0f89
  and this change's code at 4cd3b735 in turn, its arms `LOSAT_X_BPSCANBF=0` (the branch, at 4cd3b735 still that
  switch) and `LOSAT_X_BPSCANBFMIN=0` (the branch-free list) standing for the list the built-in 500 chooses; 3
  rounds, one output hash per case; UniRef50 2026_03 volume 00 `-db`, `-evalue 1e-5`): one thread, a 42-residue
  protein (1.2 % of the cells used) 7.37 -> 7.26 s and the 357-residue protein (2.3 %) 8.35 -> 8.18 s (the branch;
  that instance runs 6 G fewer instructions than dev's), 3 proteins of 1,005 residues (7.3 %) 10.33 -> 9.68 s, 10
  proteins of 4,639 residues (31.6 %) 34.27 -> 27.95 s (branch misses 1.18 G -> 0.55 G); 100 proteins of 58,139
  residues at 4 threads (two blocks, 91 % and 23 %) 70.92 -> 66.84 s; max RSS unchanged. In the strict set from
  its first pull request (no measured case slower; AGENTS.md "Speed-up switches").
- `LOSAT_X_BPFFILL=1` (BLASTP `-task blastp-fast`, speed; in the strict set from the pull request that added it,
  2026-10-25, line v4-BPF K2-2; `=0` turns it off; acts where the serial builder makes the flat compressed table of
  `LOSAT_X_BPCLUT`: one thread, natively and in WebAssembly, and any thread count in builds without threads, such as
  the serial WebAssembly module; with two or more threads the pool builds the table, `tblastx/lookup/x_clut_par.rs`,
  already by stripes of cells) fills the table by buckets of 2^10 cells
  (`tblastx/lookup/x_clut_direct.rs` `x_try_build_flat_by_buckets`, exactness argument in that file's header):
  the appends are recorded as before but not counted while they are made; they are sorted by bucket, stably,
  into the slots of `hits` that the bucket's cells will take (each packed as the query offset above the cell's
  index in its bucket), and each bucket then counts its cells, lays out their starts, PV bits and statistics,
  places its appends at their cells' cursors and puts the cells of more than 5 offsets in the scan order, within
  4 KiB of counts and the bucket's slots. The count pass's random updates of 3 MB of counts and the fill pass's
  random writes into the starts and hits (48 cycles per append, 63 % of them stalls, line v4-BPF K1) become
  passes within the cache (about 30 cycles per append, of which the partition about 15; floor H about 6). The
  record of the appends serves as the buckets' scratch and is dropped where the count-and-fill builder drops it,
  the vectors are allocated in its order and the appends per query offset are dropped after the partition, so max
  RSS does not change. An append takes at most one overflow cell, so a table whose appends stay within NCBI's
  overflow capacity cannot be refused and its count is not needed while the appends are made; the bucket fill
  declines (the count-and-fill builder runs, `LOSAT_X_REACH`) a lookup table whose appends pass the capacity, a
  query of 2^22 letters or more given to one table (the packed offset; with `LOSAT_X_BPQBLOCK` a table holds one
  block of about 50,000 letters) and the word sizes other than 5, which `-task blastp-fast` does not take. The
  table is the count-and-fill builder's slot for slot (the lib test `x_bucket_fill_is_the_count_and_fill_build`
  compares them on random queries, tables and capacities; `tests/blastp_fast_fill_buckets.rs` the reports with
  `=0`'s). The function is out of line, so the count-and-fill builder keeps its code. Not taken (line v4-BPF K2-2
  probes, one core on the workstation): the histogram counted by the appends' sink (17 more instructions per
  append in the neighbor enumeration, which is no longer vectorized), the scan position computed in the fill with
  selects instead of the move over the cells of more than 5 offsets (x1.01-1.02 the cycles), and buckets of 2^8,
  2^12 or 2^14 cells (the same or slower); in WebAssembly (the serial WASI command in Node 26, one process, 3
  rounds) the probe took x0.88 the time of the 1/10 of the human proteins x the yeast proteome and x0.98-0.99 of the
  1/10 of the Arabidopsis proteins x the zebrafish proteome, with 2 MiB less linear memory. Timed 2026-10-25 on the
  workstation (node 0, one physical core, its HT sibling idle, warmed 3 s before every run, 8 rounds of the pairs in
  alternating order; other lines' gates ran on node 0, load 14-93, so the cycle and instruction counts are the
  firmer measure; dev 17a08882; every report identical between the arms; medians, `cycles:u`): OrthoFinder's 1/10
  of the human proteins x the yeast proteome (`-db`, `-evalue 0.001`) x0.865 (instructions x1.037), of the
  Arabidopsis proteins x the yeast proteome x0.908, of the yeast proteins x the yeast proteome x0.932, of the
  Arabidopsis proteins x the zebrafish proteome x0.969; once each, the whole human proteome x the yeast proteome
  68.64 -> 64.19 s (cycles x0.933) and the whole Arabidopsis proteome x the zebrafish proteome 233.20 -> 230.09 s
  (x0.989), max RSS 83.9 -> 83.9 MiB and 76.9 -> 76.8 MiB. Where a table is small the build is a few milliseconds
  and the instructions are dev's: UniRef50 volume 00 (`-db`, `-evalue 1e-5`) with 42 residues, 3 proteins and 10
  proteins x1.0000 (cycles x1.002-1.004), exact blastp (1/10 of the Arabidopsis proteins x the yeast proteome)
  x1.0000; `=0` runs dev's instructions (x1.0000, cycles x0.992-1.031 of dev's in the same rounds). With two or
  more threads the pool builds the table and the switch does not act. In the strict set from its first pull
  request (no measured case slower; AGENTS.md "Speed-up switches").
- `LOSAT_X_BPFSEED=1` (BLASTP `-task blastp-fast`, speed; in the strict set from the pull request that added it,
  2026-10-24, line v4-BPF K2-1; `=0` turns it off; acts in the fused loop on the flat compressed
  table of `LOSAT_X_BPCLUT` where the scan lists the cells without a branch, `LOSAT_X_BPCLUTSCANBF`, natively; a
  sparser table, the reference compressed table and wasm32 keep their loops) the loop on such a table
  (`blastp/x_seed_fused_clut.rs` `x_scan_flat_fseed`) takes two changes: each hit goes first through blastp's
  branch-free two-hit test (`XHitBf::test`, the test of `LOSAT_X_BPHITBF`, which does not act here: too far and
  overlap without a branch, a flag-1 entry or a candidate on to `x_test_hit`; exact where `XHitBf::fits` holds for
  the subject, which is decided and recorded per subject, `LOSAT_X_REACH`), and no hit reads ahead a later hit's diagonal entry (with about 3.5 hits
  per cell the read-ahead rarely found one, and its bound test stayed). The cells, their order and the hits are
  `for_each_cell::<true, _>`'s, the diagonal table holds after each hit what `x_test_hit` leaves there, and the
  extension is called for the same hits with the same values, so the report is the same (the lib test
  `the_flat_loop_is_the_fused_loop_on_the_reference_table` compares the loop with the reference table's loop, and
  `tests/blastp_fast_seed_bf.rs` the reports with `=0`'s, with both lists and with the reference table). The loop
  is out of line, so the loops without it keep their code; only the choice in `x_scan_two_hit_compressed_flat`
  is new (`=0` runs the instructions of the binary before the switch on the dense tables, x1.0000, and `=0` and
  the default x0.9998 of them on a sparse one). Not taken (line v4-BPF K2 probes, one
  core on the workstation): the scan reading each listed cell's slice bounds (a list of (offset, start, end)) took
  x1.50 the cycles of 3 proteins x UniRef50 volume 00, as the branch-free list then loads the bounds at every
  offset; on the tables with the branch the loop gained nothing (x0.98-1.01 on a 42-residue protein x UniRef50);
  AVX2 batches of the cells of eight or more hits (`LOSAT_X_BPDIAG` with `LOSAT_X_BPDEFER`) x1.06 the cycles of
  the loop without them; reading the hits 4 or 16 cells ahead instead of 8, the slice bounds 8 or 32 ahead instead
  of 16, and slices without bound tests changed nothing measurable. In WebAssembly (the serial WASI command in
  Node 26, one process, 3 rounds) the loop took x1.055-1.086 the time of 3 proteins x UniRef50 (x0.945-0.979 for 10
  proteins, x0.99 for 42 residues), so wasm32 keeps `x_scan_flat_bf`. Timed 2026-10-24 on the workstation (node 0, one physical core, its
  HT sibling idle, warmed 3 s before every run, 8 rounds of the pairs in alternating order; other lines' gates ran on
  node 0, load 22-88, so the cycle and instruction counts are the firmer measure; dev a806eaa5; every report
  identical between the arms; medians, `cycles:u`): OrthoFinder's 1/10 of the Arabidopsis proteins x the zebrafish
  proteome (`-db`, `-evalue 0.001`) x0.901 (instructions x0.989, mispredicted branches x0.701), the 1/10 of the human
  proteins x the yeast proteome x0.939, the 1/10 of the yeast proteins x the yeast proteome x0.935; UniRef50 volume 00
  (`-db`, `-evalue 1e-5`) 3 proteins x0.957, 10 proteins x0.912; once each, the whole human proteome x the yeast
  proteome 74.16 -> 69.63 s (cycles x0.937) and the whole Arabidopsis proteome x the zebrafish proteome 256.49 ->
  235.51 s (x0.915), max RSS 93.1 -> 94.9 MiB and 75.6 -> 75.7 MiB (one run each). Where the switch does not act the
  instructions are those of dev: proteins of 42 and 357 residues x UniRef50 (the list with the branch) x0.9998 with
  cycles x1.003 (the `=0` arm read x1.009 in a probe round: code layout), exact blastp (1/10 of the Arabidopsis
  proteins x the yeast proteome) x1.0000. In the strict set from its first pull request (no measured case slower
  natively; AGENTS.md "Speed-up switches").
- `LOSAT_X_BPFSTAGE=1` (BLASTP `-task blastp-fast`, speed; in the strict set from the pull request that added it,
  2026-10-26, line v4-BPF s4; `=0` turns it off; acts where `LOSAT_X_BPFSEED`'s loop runs, natively, on the subjects
  that fit its scan: compressed alphabet of 15 letters, word length 2 to 6, one presence bit per cell, as on every
  table dense enough for the branch-free list, and a byte past the last word, which a search's subjects of a word or
  more keep, `LOSAT_X_REACH`; wasm32, where `LOSAT_X_BPFSEED`'s loop does not run,
  keeps its loops) changes that loop (`blastp/x_seed_fused_clut.rs` `x_scan_flat_stage`,
  exactness argument in that file's header) in three ways. The scan is out of line and takes two offsets per step
  (`tblastx/lookup/x_clut_direct.rs` `x_list_pairs`): with both new letters in the alphabet the index of `s + 1` is
  NCBI's step from the index of `s`, and the preshift after `s + 1` is the index of `s` over 225 plus the second
  letter's scaled value over 15, so the dependency chain from preshift to preshift, which set the scan's time (about
  7 cycles per offset at 26-36 instructions, inlined in the larger loop with its state on the stack), takes one
  multiply per two offsets; the divisions are exact below 15^6 (lib test), the letter's over 15 is a multiply by the
  inverse of 15 modulo 2^32, and the presence bits' shift is the constant one. Each cell's hits are staged, in scan
  order, in runs of 8 copied whole as `q - s` and `s + diagonal offset`, so that no loop ends at a cell's last hit
  (about one mispredicted branch per cell, half of the search's mispredicts on a dense table); when fewer than 8 of
  the 1,024 places are left, and at the end, one loop, out of line, gives the staged hits in order to the loop's
  test: the diagonal is `(q - s) & mask` and the query offset `(q - s) + s`. The read-ahead of the slice bounds and of the hits is that of `for_each_cell`, without its
  bound tests (the list holds 16 entries of cell 0 past its end). The hits, their order, the tests, the entries and
  the calls are those of the loop without the switch, so the report is the same (the lib test
  `the_flat_loop_is_the_fused_loop_on_the_reference_table` compares the loop with the reference table's loop,
  `the_fused_loop_takes_either_loop_with_the_same_outcomes` the choice with `x_scan_flat_bf` on subjects of up to
  3,000 letters, several flushes, and `the_staged_loop_copies_the_last_cell_of_the_hits` the copy at the end of the
  hits array with the loop without the switch; `flat_scans_give_the_reference_scans_hits` and
  `the_pairs_scan_lists_word_length_six_and_declines_outside_its_range` the list with `for_each_cell`'s, and
  `the_pairs_scan_divisions_are_exact_below_fifteen_to_the_sixth` the divisions; `tests/blastp_fast_stage.rs` the
  reports with `=0`'s). The loops are out of line,
  so the loops without the switch keep their code. Not taken (line v4-BPF s4 probes, one core on the workstation):
  the staged loop inlined with the scan (its copies were not vectorized and its state went to the stack: +18
  instructions per hit, no gain although the mispredicts halved), runs of 16 hits, the test per 256 staged hits, no
  read-ahead of the hits, and the reviewed first form of this pull request, with the letter over 15 by NCBI's
  reciprocal multiply and the presence bits' shift read from the table (x1.014 the cycles of this form on the 1/10
  of the Arabidopsis proteins x the zebrafish proteome, x1.12 on 3 proteins x UniRef50: the scan's state went to the
  stack). In WebAssembly (probe, the serial WASI command in Node 26, one process, 3 rounds, in place
  of `x_scan_flat_bf`, which every wasm32 table takes) the loop took x0.867 and x0.892 the CPU time of the 1/10 of
  the Arabidopsis and human proteins x the zebrafish and yeast proteomes, but x1.131 of 3 proteins and x1.182 of 42
  residues x UniRef50 volume 00 (x0.982 for 10 proteins), so wasm32 keeps its loops. Timed 2026-10-26 on the
  workstation (node 0, one physical core, its HT sibling idle, warmed 3 s before every run, 8 rounds of the pairs in
  alternating order; load 0.4-17; dev 70d75af5; every report identical between the arms; medians, `cycles:u`):
  OrthoFinder's 1/10 of the Arabidopsis proteins x the zebrafish proteome (`-db`, `-evalue 0.001`) x0.923
  (instructions x0.992, mispredicted branches x0.519), of the human proteins x the yeast proteome x0.943, of the
  yeast proteins x the yeast proteome x0.946; UniRef50 volume 00 (`-db`, `-evalue 1e-5`) 3 proteins x0.795
  (instructions x0.740: the scan is most of that search), 10 proteins x0.895; once each, the whole human proteome x
  the yeast proteome 62.07 -> 58.80 s (cycles x0.946) and the whole Arabidopsis proteome x the zebrafish proteome
  222.57 -> 204.97 s (x0.920), max RSS 82.7 -> 83.9 MiB and 75.2 -> 75.7 MiB (NCBI's 95 and 115 MiB). Where the
  switch does not act the instructions are dev's: a 42-residue protein x UniRef50 (a sparse table, the list with the
  branch) x1.0000 (cycles x0.999), exact blastp (1/10 of the Arabidopsis proteins x the yeast proteome) x1.0000;
  `=0` runs dev's instructions (x1.0001-1.0006, cycles x0.993-1.001 in the same rounds). With 16 threads
  (node 1, CPUs 32-47, once, 3 rounds) the whole human proteome x the yeast proteome took 5.47 -> 5.23 s (cycles
  x0.951) and the whole Arabidopsis proteome x the zebrafish proteome 18.19 -> 16.98 s (x0.921). In the strict set
  from its first pull request (no measured case slower natively; AGENTS.md "Speed-up switches").
- `LOSAT_X_BPFLEAN=1` (BLASTP `-task blastp-fast`, speed; in the strict set from the pull request that added it,
  2026-10-26, line v4-BPF s6; `=0` turns it off; acts where `LOSAT_X_BPFSTAGE`'s loop runs, on the same subjects,
  and records its sides where that switch does, `LOSAT_X_REACH`; with `LOSAT_X_BPFSTAGE=0` it does not act) runs
  that loop as `x_scan_flat_lean` (`blastp/x_seed_fused_clut.rs`, exactness argument in that file's header), which
  reads differently but stages and tests the same. The scan's two-offset steps are a loop of their own, out of line
  (`tblastx/lookup/x_clut_direct.rs` `x_list_pairs_lean`, `x_pairs_loop`), so that the offset, the list's length,
  the preshift and the next letter stay in registers (in one loop with the other steps the list's length was kept on
  the stack, loaded and stored at every pair), with NCBI's reciprocal for 15 as the constant `x_pairs_fit` asks for,
  `u32` arithmetic and the subject read through a pointer a word length in; the first word, a letter outside the
  alphabet and the last offset take `x_list_pairs`'s steps, so the list is the same. The copy is out of line and
  makes no call (`x_lean_copy`): it stages runs until fewer than 8 of the 1,024 places are left or the list ends and
  returns where it stopped (within a cell too); the caller gives the staged hits to the same test loop
  (`x_stage_test`) and calls it again. A run's second half of 4 hits is read from `k + 4` only where the cell has
  more than 4 hits left, else from `k` again (places past the run's length, left unread as before): in the probe
  72 % of the copy's L2 misses were that half's line, unused by cells of 3.5 hits on average. The read-ahead also
  takes the slice bound after the one it reads ahead (a line crossed 1 time in 16). The runs, their places, the
  flushes, and so the hits, their order, the tests, the entries and the calls are those of `x_scan_flat_stage` (the
  lib tests `the_lean_loop_is_the_staged_loop`, which replays the staging to show a flush within a cell and runs on
  either side of the 4-hit rule, `the_flat_loop_is_the_fused_loop_on_the_reference_table` and
  `the_staged_loop_copies_the_last_cell_of_the_hits` compare the loop with the staged loop, the reference table's loop
  and the loop without `LOSAT_X_BPFSTAGE`; `flat_scans_give_the_reference_scans_hits` and
  `the_pairs_scan_lists_word_length_six_and_declines_outside_its_range` the list with `for_each_cell`'s;
  `tests/blastp_fast_stage.rs` the reports with `=0`'s). The new loops are out of line and the lean loop makes its
  stage itself, so the loops without the switch keep their machine code: `x_scan_flat_stage`, `x_stage_test`,
  `x_list_pairs`, `x_scan_flat_fseed_loop` and `x_scan_flat_bf` have dev's instructions (objdump of both release
  builds; a first form that also made its stage with `XStage::new` left that function out of line in
  `x_scan_flat_stage`, 29 instructions fewer there). Not taken (line v4-BPF s5 probes, one core on the workstation):
  the steps' loop inlined in the scan (the list's length on the stack: x0.93 of the cycles on 3 proteins x UniRef50,
  against x0.82 out of line, and within 2 % on the OrthoFinder tables), a read-ahead of the line of a later cell's last hit (x1.04 on 3 proteins x UniRef50), and
  read-ahead distances of 16 and 32 cells (no difference). Timed 2026-10-26 on the workstation (node 0, one physical
  core, its HT sibling idle, warmed 3 s before every run, 8 rounds of the pairs in alternating order; load 15-69, other
  lines' builds and gates on node 0; dev d40021ba; every report identical between the arms; medians, `cycles:u`):
  OrthoFinder's 1/10 of the Arabidopsis proteins x the zebrafish proteome (`-db`, `-evalue 0.001`) x0.981
  (instructions x0.993), of the human proteins x the yeast proteome x0.988, of the yeast proteins x the yeast proteome
  x0.992; UniRef50 volume 00 (`-db`, `-evalue 1e-5`) 3 proteins x0.817 (instructions x0.807), 10 proteins x0.927;
  once each, the whole human proteome x the yeast proteome 62.18 -> 62.09 s (cycles x0.984) and the whole Arabidopsis
  proteome x the zebrafish proteome 218.30 -> 215.26 s (x0.986), max RSS 83.9 -> 84.0 MiB and 75.4 -> 77.5 MiB
  (NCBI's 95 and 115 MiB). On the 1/10 of the Arabidopsis proteins (perf record, the same core) the copy took
  8.62 -> 7.20 G cycles with 97.7 -> 53.9 M L2 misses, the scan 2.03 -> 1.46 G. Where the switch does not act the
  instructions are dev's: a 42-residue protein x UniRef50 (a sparse table) x1.0000 (cycles x1.000), exact blastp
  (1/10 of the Arabidopsis proteins x the yeast proteome) x1.0000; `=0` runs dev's instructions (x1.0000-1.0004, cycles x0.994-0.999 in a confirming run of 4 rounds of the final
  form, load 31-37, where the switch took x0.987, x0.992 and x0.815 on the first three cases above). With 16
  threads (node 1, CPUs 32-47, once, 3 rounds) the whole human proteome x the yeast proteome took 5.25 -> 5.19 s
  (cycles x0.991) and the whole Arabidopsis proteome x the zebrafish proteome 17.01 -> 16.75 s (x0.989), max RSS
  123.5 -> 124.2 MiB and 106.0 -> 108.2 MiB (medians; per run 122.8-125.5 and 104.1-108.8 MiB, dev's 122.7-123.7
  and 105.0-106.9). In the strict set from its first pull request (no measured case slower natively; AGENTS.md
  "Speed-up switches").
- `LOSAT_X_BPLKPFLAT=1` (BLASTP `-task blastp-fast`, speed; in the strict set since 2026-10-17, on by default
  together with `LOSAT_X_BPCLUT`; acts only where the search's pool builds the flat table of
  `LOSAT_X_BPCLUT`, with more than one thread) changes how that build is scheduled
  and where its vectors come from, not the table (`tblastx/lookup/x_clut_par.rs`, exactness argument in its
  header): (1) the appends pass and the stripes pass each run in one `rayon::scope` with one task per thread,
  each task taking the next chunk or stripe from one counter and storing its result at its index, instead of
  `par_iter`'s tree of joins (on the workstation at 16 threads its last wake-ups took 0.6-0.8 ms per pass); (2) in the block
  search of `LOSAT_X_BPQBLOCK` a block's table takes the `starts` and `hits` vectors of the previous block's
  table (kept by the block search, `blastp/x_query_blocks.rs` `XBlockSubjects`), cut or grown to their
  lengths and not zeroed (the build writes every entry it keeps before any is read), instead of a zeroed
  allocation per block; the query chunks of a split batch, which run inside their block's search and drop
  their tables there, take none, and a block whose appends pass the bound (`TooLarge`, the reference
  builder's table) drops them; (3) with `LOSAT_X_BPQBLOCK` the contexts' Karlin blocks (the values of
  `prepare_blosum62_lookup_query_for_word_size`, one `x_context_karlin` per context) are extra tasks of the
  appends pass, after its chunks, instead of a pass of their own on the pool; the contexts whose tasks did not
  run (where the build does not take these passes) are computed after the build by the thread that called
  it. The table's vectors are allocated by the thread that runs the build or recycled, never by the passes'
  tasks. Memory: one table's `starts` and `hits` are kept between blocks (at most the largest block's). At one
  thread nothing changes.
  Measured 2026-10-16 on the workstation (node 0 CPUs 0-15, medians of alternating rounds, every output hash
  as expected; binary 31d2bddf, a merge of this change and `LOSAT_X_BPSUBJRUN`, `LOSAT_X_BPCLUT` on in every
  arm; `logs/combo_ws_tables.md` of the perf-mb task): human x yeast `-subject`, `blastp-fast`,
  `LOSAT_X_BPSUBJRUN` off in both arms, off -> on: at 16 threads with blocks of 50,000 letters wall 8.713 ->
  7.787 s (0.893, 6 rounds), the lookup stage 2.812 -> 1.885 s, max RSS 103 -> 102 MiB; with blocks of 20,000
  letters 11.922 -> 10.236 s (0.858), RSS 85 -> 75 MiB; at 8 threads 13.445 -> 12.602 s (0.938, 4 rounds) and
  16.755 -> 15.498 s (0.928), RSS 101 -> 93 and 77 -> 67 MiB. With `LOSAT_X_BPSUBJRUN` on as well (both
  together vs both off): 0.784 (blocks of 50,000) and 0.741 (20,000) at 16 threads, 0.849 and 0.831 at 8.
  `-db` at 16 threads, both on vs both off: 0.950 (2 rounds). Exact `blastp` 1.012 (2 rounds; no compressed
  table). At one thread the serial build runs and the switch changes nothing.
  With `LOSAT_X_STATS` one line at exit when the pool built a table:
  `[X_STATS] BPLKPFLAT tables= recycled= recycled_hits=` (`tables`: flat tables the switch's pool build made
  (each search's, each query chunk's and, with `LOSAT_X_BPLKPAHEAD`, each build ahead); `recycled` /
  `recycled_hits`: builds whose `starts` / whose `hits` are the buffer of an earlier table, `hits` only where
  it had the room; with `LOSAT_X_BPLKPAHEAD` on, also builds ahead on the second spare's vectors, which
  `x_make_room` allocated or moved there, since the reuse test sees only capacity). The unit tests of
  `x_clut_par.rs` build on recycled vectors and on vectors of garbage (`u32::MAX`, `i32::MIN`) and compare
  every table with the serial build's and the reference's; those of `x_query_blocks.rs` and
  `x_lookup_ahead.rs` compare the contexts made with the Karlin tasks.
- `LOSAT_X_BPLKPAHEAD=1` (BLASTP `-task blastp-fast`, and exact blastp since 2026-10-20 (perf/bp-lkp-exact), speed; in the strict set since 2026-10-17 (perf/perfmb-default-bundle, evidence at the end of this entry): on unless it is `0`; acts only in the block search of `LOSAT_X_BPQBLOCK`, with `-subject`
  or `-db`, with more than one thread: for blastp-fast where `LOSAT_X_BPLKPFLAT` acts on the pool's flat table
  of `LOSAT_X_BPCLUT`; for exact blastp where its table is built in place
  (not with `LOSAT_DIAGNOSTICS`, whose reference
  builder prints its phases)) builds the next query block's query frames, lookup
  table and contexts during the block before it, as one task on the search's pool, instead of at the start
  of the next block's own search (`blastp/blast_engine/x_lookup_ahead.rs`, exactness argument in its
  header): the block loop runs each block's search in a `rayon::in_place_scope`, and when the block's own
  table is ready (at once when it was itself built ahead) the next block's build is spawned there; the next
  block's search then takes the values instead of computing them. They are the values of the same calls
  (`blastp_query_frames`, `blastp_compressed_lookup`; exact blastp: `blastp_aa_lookup`, a leaf task without
  rayon) on the same records, batches and pool, so the reports do not change. A build ahead that fails, refuses, would need the reference builder (`TooLarge`) or lacks
  the room in its vectors gives nothing, and the block builds inline, so an error is raised where it is
  without the switch; a panic of the build is kept and resumed where the block's search would start (with
  `panic = "abort"`, the release profile, it ends the process when it happens). Memory: one more table
  (`starts`, `hits`, PV array) and one block's frames and contexts are alive while the next block's table
  is built; every vector of that table is allocated by the block search's thread, never by a pool task:
  that thread keeps a second spare (`x_query_blocks.rs` `XBlockSubjects`) with the room for the largest
  table so far plus a half (new vectors of twice its hits where they are short; untouched capacity is not
  resident), which the build ahead takes whole (`XFlatArgs::x_ahead`). Exact blastp: one more table
  (backbone 0.44 MiB, overflow 6-8.5 MiB at 50,000-letter blocks, PV 3.5 KiB), on vectors the block search's
  thread allocates and recycles (`lookup/x_lut_direct.rs` `XAaSpares`: the same room rule for the overflow;
  each block's table is given back at the end of its search), with that thread's neighbour lists moved to the
  builds ahead (one set; given back to it before a block that builds its own table and at the last block). At
  one thread nothing changes.
  Measured 2026-10-16 on the workstation (node 0, 16 physical cores, siblings idle; Species4 hsap x Species6 scer
  blastp-fast `-subject`, #199's blastp switches + `LOSAT_X_BPCLUT` `LOSAT_X_BPREDOPAR` `LOSAT_X_BPLKPFLAT`
  `LOSAT_X_BPSUBJRUN` in every run; alternating pairs, medians; binary built from a local merge of this branch
  with main 61b2a124, for `LOSAT_X_BPSUBJRUN`): `LOSAT_X_BPQBLOCK=50000` 16 threads 6.858 -> 6.476 s
  (x0.944), max RSS 101.8-103.9 -> 125.9-128.1 MiB, `blastp_build_lookup` 1.89 -> 0.01 s; `=20000` 8.904 -> 8.139 s
  (x0.914), 75.4-77.1 -> 93.3-95.4 MiB; 8 threads x0.980 / x0.965. The second table is 21.4 MiB at 50,000-letter
  blocks (largest table: `hits` 18.4 MiB, `starts` 2.9 MiB, PV 0.1 MiB) and 13.2 MiB at 20,000 (largest
  `hits` 10.3 MiB), counted on the hsap x scer queries with 8 threads. That accounts for the 50,000 arm's +24 MiB
  (within about 3 MiB), but not for all of the 20,000 arm's +18 MiB (3.5-5 MiB more than one table) nor for
  `-db`'s +34 MiB without `LOSAT_X_BPDBCUT` (13 MiB more). The rest is not attributed by measurement; the most
  likely cause is the build's own transient buffers (up to 19 MiB of packed appends at 50,000 letters, frames,
  contexts), which now coincide with the current block's search peak (with `-db`, including the scan's
  per-thread blocks). `-db` (`BlastDBSpecies6`, 1 MiB database blocks) without `LOSAT_X_BPDBCUT`: 2 pairs
  inconclusive (pairs x1.143 and x0.987; the job has two speed modes; the scan runs 3 of 16 threads), max RSS
  93-95 -> 128-129 MiB. With `LOSAT_X_BPDBCUT=1` (binary from a local merge of this branch, main 61b2a124 and
  #231; `-db`, 16 threads, 4 alternating pairs): 7.072 -> 6.602 s (x0.932, pairs 0.928-0.937), `blastp_build_lookup`
  summed over the blocks 1.85 -> 0.015 s, max RSS 92.9-95.8 -> 116.4-120.6 MiB (+23 MiB), all outputs
  equal (sha f4ecd9f36f7f2ea1).
  Default ON since 2026-10-17 (perf/perfmb-default-bundle; workstation (2x Xeon Gold 6338, 16 physical cores of one
  node with their siblings idle, alternating rounds, every output the expected one; perf-mb
  `logs/defaultset_ws_report.md` unless named); a reduction is 1 - with / without the switch, the other switches of
  the candidate set on): blastp-fast of 23,879 human proteins against the yeast proteome (`-subject`) at 16 threads
  -6.1 % wall (all four rounds at a matched clock, -4.1..-14.6 %), at 4 threads -1.3 %; with `-db` and
  `LOSAT_X_BPDBCUT` x0.932 (four pairs, 0.928-0.937, perf-mb `logs/lkpdb_report.md`). Max RSS +23-28 MiB, below
  NCBI's (123 vs 328 MiB `-subject`, 116 vs 582 MiB `-db`). Exact blastp since 2026-10-20 (perf/bp-lkp-exact;
  workstation node 1, 16 physical cores, siblings idle, 3 s warm-up of the booked cores before each run,
  alternating pairs, medians, every output the same; `v4-B/results/lkpab/`): before, each block's table was
  built on the block's thread alone while the pool waited (about 14 ms per 50,000-letter block, 60 ms for the
  first with a cold neighbour memo: 0.37 s of 4.10 s for 3,000 human proteins x the yeast proteome at 16
  threads); with the switch only the first block builds its own (0.06 s). 3,000 human proteins x yeast
  (`-subject`, outfmt 6) 16 threads 4.105 -> 3.750 s (x0.914, 8 pairs 0.912-0.919), 8 threads 7.600 -> 7.285 s
  (x0.959, 8 pairs), 2 threads 28.81 -> 28.59 s (x0.992, 3 pairs); outfmt 0 at 16 threads 4.170 -> 3.820 s
  (x0.916, 8 pairs); the yeast proteome against its own `-db` 16 threads 21.10 -> 20.25 s (x0.959, 4 pairs),
  8 threads 36.46 -> 35.72 s (x0.980, 4 pairs); every switch-on run faster than every switch-off run of its
  case. Max RSS +3-5 MiB (16 threads 48.7-50.4 -> 53.2-55.0 MiB; NCBI 327 MiB). The builds ran on one pool
  thread, took as long as inline (median 13.9 ms) and ended at least 81 ms (median 145 ms) before the end of
  their block's search, so they did not overlap its report (the binary had #328's `LOSAT_X_PAIRWISEPAR`,
  which writes outfmt 0 on the pool) and the block never waited for one. Not measured: blastp-fast `-db` on a database that `LOSAT_X_BPDBCUT` leaves in 1 MiB blocks (at least 4 MiB per thread); the only pairs near it (the 2.9 MB database with the cut off) read x1.143 and x0.987 (perf-mb `logs/lkpahead_report.md`); a Swiss-Prot blastp-fast pair is in the pending workstation plan (`ws/defbundle_plan2.txt`).
  With `LOSAT_X_STATS` one line at exit when a counter is not zero: `[X_STATS] BPLKPAHEAD blocks= prebuilt=
  dropped=` (`blocks`: blocks searched where the switch acts; `prebuilt`: blocks that took values built ahead;
  `dropped`: builds ahead that gave none, whose blocks built inline; every block but the last starts one build
  once its own table is ready, if that table is a flat one or (exact blastp) one built in place, as it is
  unless the reference builder made it, so `prebuilt + dropped = blocks - 1`). A unit test of
  `x_lookup_ahead.rs` builds a block's values ahead and inline, for exact blastp and blastp-fast, and compares
  the query frames, the table and the contexts (the Karlin blocks bit for bit).
- `LOSAT_X_BPDIAG=1` (BLASTP, speed; in the strict set since 2026-10-17; `=0` turns it off; acts only in
  the fused loop) the two-hit
  diagonal test of the hits of one lookup cell in batches of eight instead of one after the other
  (`blastp/x_seed_fused.rs`, exactness argument in its header): for each hit the entry's new value and whether
  the hit is a candidate are computed without branches on the entry in AVX2 registers (the eight entries read
  with scalar loads, compared in i32 lanes, blended; a cell's last batch of fewer than eight hits takes a
  scalar batch), all eight entries are stored, then the candidates go to the context check and the extension
  in the cell's order. Exact because a cell lists a query offset at most once (`s_AddWordHits`; checked for
  both BLASTP builders by `x_lut_direct` test `every_cell_lists_a_query_offset_once`, and per batch by a
  `debug_assert`), so a batch's hits fall on distinct entries, and a candidate's entry is stored back unchanged
  before its call. The 31-bit truncation and the i32 wrapping are those of the per-hit test. Without AVX2
  (another architecture, or an x86-64 CPU without it) the switch keeps the per-hit loop, because the scalar
  batch was slower than the per-hit loop on both machines (workstation 12.56 s against 11.59 s, local about
  7.3 s against 6.8 s); `LOSAT_X_BPDIAG_SIMD=0` forces the scalar batch (diagnostic; the tests use it). The
  AVX2 batch does not use the AVX2 gather: on the workstation (Xeon 6338, Gather Data Sampling microcode
  mitigation on) the gather was slower than eight scalar loads; locally both were equal. The batches prefetch
  nothing (the prefetch of 16 hits ahead is the per-hit loop's). With `LOSAT_TIMING=1`
  each search of a query set prints `[TIMING] x_bpdiag: mode=<avx2|avx2-defer|avx2-defer-masked|scalar|per-hit> cells=<n>
  batched_hits=<m> candidates=<k>` (the path that ran, cells tested in batches, their hits, the hits passed on; zeros with
  `per-hit`). Measured 2026-10-14 on the workstation (one pinned core, alternating, the other 8 exact switches
  and `LOSAT_X_BPQBLOCK=50000`, x human proteome by `-subject`, rows identical in every arm; `blastp_scan_ungapped`
  (wall) in s): human block 1 per-hit 11.59 (20.61), gather 10.17 (18.90), scalar loads 9.13 (17.86), scalar
  batch 12.56 (21.53); C. elegans block 5 per-hit 32.68 (70.68), gather 28.55 (66.16), scalar loads 26.45
  (64.62).
- `LOSAT_X_BPDEFER=1` (BLASTP, speed; in the strict set since 2026-10-17; `=0` turns it off; acts only
  where `LOSAT_X_BPDIAG`'s AVX2 batch runs, not with `LOSAT_X_BPDIAG_SIMD=0` or
  without AVX2) the candidates of a lookup cell's batches of eight are called after the cell's last batch of
  eight instead of at the end of each batch (`blastp/x_seed_fused.rs`, exactness argument in its header): each
  batch moves its candidates' query offsets to the front of its register (a permutation from a table of the 256
  lane masks), stores the register into a buffer and advances the buffer by the number of candidates, so the
  batch loop has no branch on the candidates; then the closure is called for the buffer in order, and the
  cell's last batch of fewer than eight hits is the scalar batch as before. Exact because the hits of one cell
  lie on distinct entries and a call reads and writes only its candidate's entry and the HSP list, so the later
  batches of the cell touch nothing a call reads or writes, and the calls keep their order. With
  `LOSAT_TIMING=1` the `x_bpdiag` line names `mode=avx2-defer` (with `LOSAT_X_BPMASKED=0`; by default
  `avx2-defer-masked`, see that entry), with the counts of the AVX2 batch. Measured 2026-10-16 on the
  workstation (perf stat, user-space cycles; the other nine blastp switches and `LOSAT_X_BPQBLOCK=50000` in
  both arms, x the subject proteome by `-subject`, one thread, outputs identical in every arm): on stage7's 49
  OrthoFinder jobs (every 10th query protein of each species x the whole subject proteome, 28 jobs side by
  side, both arms in one run) the summed cycles were 0.938 of the run without it (per job 0.925-0.965, median
  0.936) and the instructions 0.985; on one core in alternated pairs Blast0_6 0.920 (`blastp_scan_ungapped`
  45.0 -> 40.2 s), human block 1 0.950 and C. elegans block 5 0.963 (its gapped stages are over half of the
  job), with 13-18 % fewer branch mispredicts.
- `LOSAT_X_BPMASKED=1` (BLASTP, speed; in the strict set from the pull request that added it, 2026-10-24;
  `=0` turns it off; acts only where `LOSAT_X_BPDEFER`'s deferred calls run: not with `LOSAT_X_BPDEFER=0` or
  `LOSAT_X_BPDIAG_SIMD=0`, not without AVX2 (aarch64, WebAssembly), not with
  `LOSAT_X_BPDIAG=0`, never in blastp-fast's loops on the compressed table (they have no AVX2 batches); and
  only on lookup tables with at least `LOSAT_X_BPMASKEDMIN` hits per used backbone cell, in hundredths,
  built-in 250, `0` every table) blastp's
  fused loop takes the masked batch of blastx's fused loop (`LOSAT_X_BXDEFER` until its fold on
  2026-11-02) in place of `LOSAT_X_BPDEFER`'s deferred calls
  (`blastp/x_seed_fused.rs`, `x_cell_avx2_defer_masked`; exactness argument in its module header): the same buffer and calls, with the batches of eight taking their lanes' entry coordinates from
  scalar loads of the hits, and a cell's last batch of fewer than eight hits as one masked AVX2 batch whose
  candidates join the buffer instead of the scalar batch after the calls. The calls (each with the entry its
  test read) and the table are those of the deferred calls. A table with fewer hits per used cell than
  `LOSAT_X_BPMASKEDMIN` (a short query) keeps the deferred calls: there the masked batch and the buffer cost
  more than the scalar batch they replace. The choice is made once per thread and table, with
  `LOSAT_X_BXDIAGMIN`'s measure. With `LOSAT_TIMING=1` the `x_bpdiag` line names `mode=avx2-defer-masked`
  (the counts of every table), and a second line `[TIMING] x_bpmasked: min_hits_per_cell=<n>
  tables_masked=<a> tables_deferred=<b> hits_per_cell=<least>-<most>` gives the rule's (thread, table)
  choices. Measured 2026-10-24 on the workstation (one pinned core, user-space cycles, medians of paired
  per-round ratios, every output identical in every arm), against `dev` before the switch: stage7 Blast6_6 and
  Blast0_6 (`-subject`, about 100 hits per used cell) x0.954 and x0.941 (the cell test 7.28 -> 6.32 cycles per
  hit), 300 yeast proteins x the fly proteome (`-db`) x0.955, 567 human proteins x Species4 x0.957, WSSV x
  PajaWSV x0.984 (its redo dominates), 3 proteins (2.60 hits per used cell) x UniRef50 volume 00 (`-db`) x0.933
  and x the fly proteome x0.989; the tables below the cutoff take the deferred calls and run as before: one
  42-residue (1.11) and one 357-residue protein (1.48) x UniRef50 volume 00 x1.006 and x1.000, 2 proteins
  (2.08) x the fly proteome x1.006 (0.982-1.024, a 0.55 s run). With the masked batch on every table
  (`LOSAT_X_BPMASKEDMIN=0`) against `LOSAT_X_BPMASKED=0` the first k of 10 UniRef50 query proteins against
  the fly proteome took x1.148 at 1.11 hits per used cell (x1.192 against UniRef50 volume 00), x1.123 at 1.48
  (x1.071), x1.020 at 2.08, x0.987 at 2.60, x0.924 at 4.52, x0.928 at 7.07 and x0.935 at 10.38 (x0.919 against
  UniRef50 volume 00), so the built-in cutoff lies between 2.08 and 2.60.
- BPSCANBF (graduated 2026-10-26, in the set since 2026-10-17; no variable; always on) (BLASTP, speed; acts in
  the fused loop on the standard lookup table; the compressed table keeps
  its scan) blastp's fused loop finds the lookup cells with hits by the
  branch-free cell mode of BXSCANBF (graduated 2026-10-26; blastx always takes it) (`blastp/x_seed_fused.rs` `x_scan_cells_bf`: the positions of
  256 words listed with tblastn's `x_block_candidates`, the list advanced by the presence bit, and the listed
  cells passed to the cell test in position order), in every cell-test mode (`LOSAT_X_BPDIAG`'s AVX2 batch
  with or without `LOSAT_X_BPDEFER`, its scalar batch, the per-hit test). The cells, their order and their
  hit lists are the presence-test loop's (the tests of BXSCANBF), so the report does not change.
  The instance (`x_scan_two_hit_bf`) is out of line and calls the cell test through a closure of
  its own type, as it did while the loop with the branch was a second instance (with one closure type the
  cell test had two callers, the compiler no longer inlined it, and that instance ran 1.9-6 % more
  instructions); the loop with the branch is gone for the standard table, and `x_scan_two_hit_in_with::<false, _>`
  serves blastx's loop and the tests.
  Timed 2026-10-17 on the workstation (perf-si-20261012 session SCANBF-1, node 1, one binary without and
  with the switch and the build without it, in turn; one output hash per case, NCBI BLAST+ 2.17.0's
  included where run): one 357-residue human protein against UniRef50 2026_03 volume 00 (`-db`, 2,178,834
  sequences, 3.0 G residues, `-evalue 1e-5`) 0.75 at one thread (69.70 -> 52.42 s; NCBI 71.67 s) and 0.76
  at 8 (9.19 -> 6.94 s; NCBI 9.76 s), branch misses 2.64 G -> 1.44 G, cycles 220 G -> 165-171 G; a
  42-residue protein against the same volume 0.64 at one thread (21.03 -> 13.40 s, cycles 65.2 G -> 41.0
  G; NCBI 11.17 s). With many query proteins the presence bit is set at most words and the branch is
  predicted: 567 human proteins x the human proteome (an OrthoFinder job's first block, `-subject`,
  `-evalue 0.001`) 1.00 at one thread, every 10th Arabidopsis protein (2,763) x the yeast proteome 1.00,
  10 proteins x volume 00 at 8 threads 0.995. In the strict set from its first pull request (no measured
  case slower; Owner decision of 2026-10-17 that the branch-free scans join the set once each shows no
  loss, and plan v4's rule that a switch that wins everywhere joins in the pull request that adds it).
- `LOSAT_X_BPHITBF=1` (BLASTP, speed; in the strict set from the pull request that added it, 2026-10-20;
  `=0` turns it off) the per-hit loop of the fused seed loop (`blastp/x_seed_fused.rs` `x_cell_per_hit`, the
  loop without `LOSAT_X_BPDIAG`'s batches: the default where the CPU has no AVX2, aarch64 and WebAssembly;
  exactness argument in the module header) takes the two outcomes most hits of a flag-0 entry end in, too far
  (`diff >= window`: `last_hit = s + diag_offset`) and overlap (`diff < wordsize`: nothing), without a branch,
  in doubled arithmetic: with the packed entry `raw`, `d2 = 2 (s + diag_offset) - (raw <<
  1)` is twice the reference's `diff` (a cleared entry holds `-window`), a candidate is `(d2 - 2 wordsize) as
  u32 < 2 (window - wordsize)`, and the entry's next value is selected. A flag-1 entry or a candidate goes on
  to `x_test_hit`. The loop also reads ahead (16 hits, x86-64) in a loop of its own over the
  hits that have a later hit, without a bound check per hit. Acts per subject where `1 <= wordsize < window`,
  `0 <= diag_offset` and `diag_offset + length + 2 window < 2^30` (otherwise the loop without the test runs: a
  window equal to the word size, or a subject of about 2^29 letters or a `-window_size` near 2^29). The
  AVX2 batches of `LOSAT_X_BPDIAG` (the x86-64 default) and blastp-fast's loop on the reference compressed
  table keep their own tests; blastp-fast's loop on the flat table takes the test under `LOSAT_X_BPFSEED`
  (`LOSAT_X_BPHITBF=0` does not turn it off there). With `LOSAT_X_BPDIAG=0` 72-74 % of the hits of WSSV x
  PajaWSV and of 300 yeast
  proteins x the fly proteome are too far, 18.6-18.8 % overlap, 2.5-3.7 % meet a flag-1 entry and 4.7-5.3 %
  are candidates (v4-A10 probe bp1). With `LOSAT_X_STATS` (and the `xstats` feature) `main` prints
  `[X_STATS] BPHITBF subjects=N` at exit when N is not zero (N: the subjects this loop scanned).
  Measured 2026-10-20 on the workstation (node 0, one physical core, its HT sibling idle, warmed 3 s before
  every run, arm order reversed every other round, one warm-up round + 4 rounds, medians of the per-round
  ratios against dev 23b4ef09, the binary before the switch; other lines' gates ran on node 0 meanwhile, load
  12-27, so the cycle and instruction counts are the firmer measure; every report identical between the
  arms). With `LOSAT_X_BPDIAG=0` (the per-hit loop on every subject, as without AVX2): 30 yeast proteins x the
  fly proteome (`-db`, `-evalue 1e-5`) x0.824 (5.40 -> 4.45 s; cycles x0.823, instructions x0.988,
  mispredicted branches x0.546); WSSV x PajaWSV (`-subject`) x0.966 (cycles x0.965, instructions x0.993,
  mispredicted branches x0.768). At 8 threads (node 1, CPUs 32-39, their HT siblings idle, the same form, 4
  rounds): x0.835 (0.730 -> 0.611 s; cycles x0.825, instructions x0.988, mispredicted branches x0.548) and
  x0.974 (0.195 -> 0.191 s; cycles x0.959, instructions x0.993); `=0` x0.995-0.997, the default x0.995, both
  with instructions x1.0000. By default with AVX2 the batches take every subject: x0.996-1.002,
  instructions x1.0000. With `=0` the instruction counts equal those of the binary before the switch
  (x1.0000-1.0001, with and without `LOSAT_X_BPDIAG=0`). On WebAssembly (LOSAT Web's serial reactor, ABI v2,
  under Node 26 in one process: an instance per arm, the default and `=0` alternating, one warm-up search +
  8 rounds, on one booked workstation core): AP027131 x AP027133 proteins x0.977 (19.51 -> 19.11 s, 7 of 8
  rounds below 1), WSSV x PajaWSV x0.991 (6 of 8), the same output bytes.
- BPUNGAP (graduated 2026-10-26, in the set since 2026-10-14; no variable, no `[X_SHADOW] BPUNGAP` line; always
  on) (BLASTP and BLASTX one-hit and two-hit BLOSUM62 ungapped extension, speed) the word loop and the
  X-drop walks of `s_BlastAaExtendLeft` / `s_BlastAaExtendRight` read the scores from one flat 32 x 32 byte
  table (`q << 5 | s`; the sentinel's `BLAST_SCORE_MIN` stored as -128) with raw-pointer walks over the
  positions the reference indexes, NCBI's statements and exit tests unchanged (`blastp/x_ungapped_fast.rs`,
  exactness argument in its header). Acts in `extend_one_hit_blosum62` / `extend_two_hit_blosum62` (TBLASTN
  calls the matrix-generic functions; its word finder takes the same walks under `LOSAT_X_TBNUNGAP`); a
  dropoff outside 1..=127, a word size outside 1..=11 or a residue code >= 32 run the reference
  (`extend_one_hit_impl` / `extend_two_hit_impl`). The unit tests of `x_ungapped_fast.rs` compare the walks
  with the reference functions on random inputs, boundaries, ties, very long walks and sentinels in words.
- `LOSAT_X_BPUNGAPSENT=1` (BLASTP, speed; in the strict set from the pull request that added it, 2026-10-23;
  `=0` turns it off) the two-hit ungapped extension of blastp's fused seed loop for word size 3
  (`blastp/x_ungapped_sent.rs`, exactness argument in its header): a closure of its own, with the extension
  inlined, takes the candidates of the table of single letters. The checks `x_ungapped_fast.rs` makes on every
  call are made once: residue codes below 32 and every context's X-drop in 1..=127 for a query block, codes
  below 32 and a sentinel at each end for a subject. The walks have no bound test: the query block is copied
  with a sentinel added at each end (2 bytes more than the block), the subject buffer has its own, and a step
  on a sentinel (-128 in the byte table) stops the walk at the step where the reference's bound stops it, with
  the same best, length and stop. The word loop is unrolled for three letters and the walks take two steps a
  loop. A block with another word size (blastp-fast, `-word_size`), an X-drop outside 1..=127 or a code
  >= 32, and a subject with a code >= 32, run the plain closure (`extend_two_hit_blosum62`); so does `-window_size 0` (no fused loop). BLASTX has its
  own switch, `LOSAT_X_BXUNGAPSENT`; its reference extension is `extend_two_hit_blosum62`.
  Measured 2026-10-23 on the workstation (node 1, CPU 40 for one thread and CPUs 40-47 for eight, their HT siblings
  idle, warmed 3 s before every run, arm order reversed every other round, one warm-up round + 8 rounds at one
  thread and 2 rounds at eight, medians of the per-round ratios against dev 7cfec941, the binary before the switch;
  every report identical between the arms and equal to NCBI BLAST+ 2.17.0's). OrthoFinder jobs of 7 species at one
  thread (`-subject`, `-evalue 0.001`): Blast0_6 x0.880 (54.55 -> 48.05 s, 8 of 8 rounds below 1; cycles x0.880,
  instructions x0.748), Blast6_6 x0.903 (18.33 -> 16.57 s; cycles x0.903); WSSV x PajaWSV x0.978 (cycles x0.977);
  5,000 x 5,000 UniRef50 proteins (`-db`, `-evalue 0.001`, 8 threads) x0.897 (239.7 -> 214.9 s; cycles x0.891). In
  Blast0_6 the extension went from 249 to 177 cycles per candidate (288 million candidates; its share of the dev
  binary's cycles 41.5 % in `perf record`). With `=0` the instruction counts equal the binary's before the switch
  (x1.0000-1.0002) at x0.996-0.998 of its time. Max RSS within the spread of repeated runs of one binary. On
  WebAssembly (LOSAT Web's serial reactor, ABI v2, under Node 26 in one process: an instance per arm, the default
  and `=0` alternating, one warm-up search + 8 rounds, on one booked workstation core): WSSV x PajaWSV x0.984 (7 of
  8 rounds below 1), AP027131 x AP027133 proteins x0.991 (18.96 -> 18.86 s, 6 of 8), the same output bytes.
- LINKSORT (graduated 2026-10-22, in the set since 2026-10-14; no variable, no `LOSAT_X_LINKSORTSHADOW`;
  always on) (TBLASTX, speed) the index sorts of the in-place HSP list path (`tblastx/x_link_sort.rs`,
  exactness argument in its header).
  Two changes, both of how an index order is computed, not of the order: (1) the linking call sorts its
  output order once, stably, with `s_FwdCompareHSPsTransl` and ties by `s_RevCompareHSPsTransl`, where it
  sorted twice (a stable sort by B followed by a stable sort by A is one stable sort by "A, then B"); this
  also saves a sort at one thread. (2) An index sort of 32,768 indices or more runs as rayon's stable
  parallel sort (`par_sort_by`) when its call may use the search pool: the two sorts of a linking call
  under the call's `allow_parallel`, as its per-group linking; the other HSP sorts (the score sorts after
  each frame, which with `LOSAT_X_TBXKMERGE` only the reference frame loop makes, after the first
  linking call and after the re-evaluation) when the calling thread is a worker
  of a rayon pool, which in a search started by the command line or the Web adapter it is exactly when
  `-num_threads` is 2 or more (the header names the one case this repository's code does not show: a
  program that calls a one-thread search from a thread of its own rayon pool). A stable sort with a
  comparator that is a total preorder has one result, so the index order is the same at any number of
  threads; the comparators are NCBI's four, each a lexicographic comparison of integer fields. With one thread, in a
  build without the `parallel` feature and on wasm32 without `wasm-threads` every sort is the serial
  statement (`sort_today`). The old statements in `x_hsp_inplace.rs` (an index `sort_by` and two stable
  sorts) and the shadow went; the reference linking body of `LOSAT_DEBUG_CUTOFFS` / `LOSAT_DIAGNOSTICS`
  keeps its own two sorts of references.
  `LOSAT_X_STATS=1` prints
  `[X_STATS] LINKSORT sorts= indices= parallel_sorts= parallel_indices= composite_sorts=`
  (index sorts through the path and their lengths added; those that ran on the pool; final sorts made as
  one). With several subjects on a pool, a thread that waits for the parts of its sort can run another
  subject meanwhile, as it can while it waits for its linking groups, so the `sort`, `final_sorts`,
  `merge` and `reap` times of `LOSAT_TIMING` can then hold time of another subject. Measured once on
  2026-10-13 (medians of three alternating runs before and after, one binary, a loaded
  machine; `[TIMING] x_link1` and `x_subject`, seconds, before -> after): AP027131 x AP027133 (genetic
  code 4, 338,859 HSPs in the first linking call) at 8 threads `sort` 0.085 -> 0.039, `final_sorts`
  0.067 -> 0.009, `reap` 0.073 -> 0.045, `merge` 0.167 -> 0.156, and at one thread `final_sorts`
  0.060 -> 0.035 (the other sorts are the serial statement there); NZ_CP006932 against itself (516,395
  HSPs) at 8 threads `sort` 0.139 -> 0.075, `final_sorts` 0.100 -> 0.014, `reap` 0.098 -> 0.068, and at
  one thread `final_sorts` 0.105 -> 0.047. Max RSS did not change (164 -> 161 and 242 -> 243 MiB at 8
  threads, one run each). What is left of `sort` and `reap` includes the move of the records into the
  sorted order, which stays on one thread. In a search of 200 short subjects no list reached 32,768 HSPs
  and no sort went to the pool.
- `LOSAT_X_TBXKMERGE=1` (TBLASTX, speed; in the strict set since 2026-10-18, on by default; acts where
  the PAIRPAR unit path runs the units) the merge of a subject's unit lists (`XPairPar::merge_subject`,
  `tblastx/x_kmerge.rs`, exactness argument in its header) without its full sorts: each unit sorts its own
  list in its task on the pool (a compact list between the copy it returns and the task's scratch
  vector, so the sort allocates nothing); a chunk's init list is the k-way merge of its units' lists, ties
  to the lower unit, where it was the stable sort of their concatenation; the conversion's final sort
  (`Blast_HSPListSortByScore`) is a stable sort of each run of HSPs equal in score, subject start and
  subject end (the list leaves the init sort in that order; a list that is not falls back to the full
  sort); each frame list is merged into the subject's list by a backward two-way merge, where it was
  appended and the whole list sorted. A stable sort with a comparator that is a lexicographic comparison of
  integer fields has one result, so the list is the same element for element. `LOSAT_X_STATS=1` prints,
  whenever the switch is on (zeros where no unit path runs),
  `[X_STATS] TBXKMERGE subjects= chunks= tie_runs= k1_fallback= frame_unsorted=` (the two fallbacks have
  stayed 0). `LOSAT/tests/tblastx_kmerge.rs` runs it at 1-16 threads, with the units' full lists
  and on tied HSPs of two queries, at three chunks per subject frame too. Measured on the
  workstation (node 1, 16 cores with idle siblings, dev 79fab14e against the switch, alternating pairs,
  reports identical): `[TIMING] x_pairpar` merge at 16 threads (the units' own sorts are counted in
  `units` now) 2.30 -> 0.92 s for AP027131+AP027132+AP027133 x AP027133+AP027078+AP027131 (code 4),
  0.13 -> 0.05 s for AP027131 x AP027133, 0.22 -> 0.09 s for NZ_CP006932 x itself, 5.82 -> 1.66 s for
  Sakai x MG1655 (code 11); wall time, medians: the first 11.81 -> 11.68 s at 16 threads (5 pairs, ranges
  apart) and 17.29 -> 16.86 s at 8 (3 pairs), AP027131 x AP027133 2.59 -> 2.48 s at 8 (5 pairs), Sakai x
  MG1655 137 -> 124 s at 16 (3 pairs, a loaded node: 115-141 s); with the units sorting by
  `slice::sort_by` (otherwise the same) NZ_CP006932 x itself 6.35 -> 6.21 s and 6.84 -> 6.70 s at 16 and
  8 threads, and at one thread (3 pairs) x0.991, x0.982 and x0.995 for the first three inputs. Max RSS
  in the same runs: Sakai x MG1655 1,967 -> 1,945 MiB at 16 threads (2,008 with `slice::sort_by`, whose
  buffers the 16 tasks held at once), the first input 592 -> 540 and 439 -> 434 MiB, AP027131 x AP027133
  117 -> 120 MiB at 8 threads (NCBI 233 MiB at one thread). The switch-off build runs 0.001 % more
  instructions than dev (`perf stat`, one thread).
  The unit test `compact_check` in `x_pair_par.rs` compares the k-merge with the plain merge (both are the unit
  path); unit path against frame loop is checked only through the frozen NCBI bytes of the fast fixtures (at
  one thread the frame loop for batches under 5,000 letters, `LOSAT_X_PAIRPARN1`, the units above; the units at
  4 threads).
- `LOSAT_X_TBXPIPE=1` (TBLASTX, speed; in the strict set since 2026-10-19, on by default; acts only with
  two or more threads in a native build and two or more query batches, and
  not with debugging output (`LOSAT_DEBUG*`, `LOSAT_TRACE_HSP` (diagnostic build), `LOSAT_DIAGNOSTICS`,
  `LOSAT_WASI_THREADS_DEBUG`, the stage dumps (diagnostic build)); in a default build any `LOSAT_DEBUG*`
  variable still keeps the pipe off (prefix scan in `x_pipe_wanted`), while the other diagnostic-build
  variables are ignored)
  the next query batch starts when the batch before has ended its seed stage, while that batch still
  links its subjects (`tblastx/x_batch_pipe.rs`, exactness argument in its header). With `LOSAT_X_TBXNEXT`
  every subject of a batch starts at once and its units fill the pool, so after a batch's last unit only
  its subjects' chains are left (merge, both linking calls, the re-evaluation): on AP027131+AP027132+AP027133
  x AP027133+AP027078+AP027131 (code 4, three batches) at 16 threads that was 3.6 s of 11.2 s with 13 of
  16 threads waiting. The batch loop of `run_in_pool` is unchanged (checks, errors, writes, their order);
  it takes each batch's result from the pipe, which searched it when the loop reached it or started it at
  the seed-stage end of the batch before (the last subject's dispatch point 5, where
  the seed release has freed that batch's tables), on the same pool. Admission: the next batch starts
  early only while the batches past their seed-stage end hold at most 4,194,304 HSPs
  (`LOSAT_X_TBXPIPEHSPS=<n>` changes it, for tests); otherwise when the loop reaches it, as before. While a
  later batch may run (the batch is not the last one and its held HSPs are within the limit), a batch links
  each subject's HSPs on the subject's own thread instead of on the pool: a thread that waits in a pool
  call runs any job of the pool (rayon), so a linking call on the pool could wait inside another batch's
  whole subject chain and hold its own batch's end back (16 threads below: 9.54 s instead of 8.59 s). The
  last batch, and a batch over the limit, links on the pool as before. A
  batch's `LOSAT_TIMING` block is printed when the loop takes the batch. `LOSAT_X_STATS=1` prints
  `[X_STATS] TBXPIPE planned= taken= early= loop= waited= key_mismatch= max_held_hsps=`.
  On the workstation (16 CPUs of one socket, medians of alternating runs, reports equal), the three
  batches above took 8.59 s with it and 10.51 s without at 16 threads (peak memory 558 / 559 MiB), 14.00 /
  15.60 s at 8 (456 / 434 MiB) and 24.69 / 26.02 s at 4 (365 / 327 MiB); two batches (AP027131+AP027132 x
  AP027133+AP027078) 3.84 / 4.17 s at 16 (294 / 330 MiB). One-batch inputs (MelaMJNV x 9 MJNV genomes at
  16 threads, a one-genome pair at 8) are equal in time and memory. The memory cost is the next batch's
  tables and HSPs while the batch before links (up to +38 MiB at 4 threads here).
- `LOSAT_X_TBXBOARD=1` (TBLASTX, speed; in the strict set from the pull request that added it,
  2026-10-27; `=0` turns it off; acts in a native build with two or more threads on every query batch: one
  subject, or several in `LOSAT_X_TBXNEXT`'s loop (with `LOSAT_X_TBXNEXT=0` a batch of several subjects
  runs as before); no size or count rule) one work board per search
  (`tblastx/x_work_board.rs`, exactness argument in its header) hands out the batch's pure work items:
  every linking call's groups (`x_link_slices.rs` `x_link_board`), the chunk sorts and merges of a chain's
  stable sorts of at least the parallel size (`x_link_sort.rs` `x_board_sort`), a subject's frame merges
  (`x_kmerge.rs`), and the subjects' units (`x_pair_par.rs` `run_units`). A thread that takes from the
  whole board (one of the subject loop without a subject, or one that waits for a batch) takes groups,
  sorts and merges of any batch first, then units by batch and subject index (a batch without the unit
  path has only groups, merges and sorts). A subject's thread puts its units on the board (a batch's only
  subject also spawns up to threads - 1 helper jobs that take its units), takes its own units only, then
  blocks until the units others took have run (taking groups and earlier subjects' units first left it
  blocked on its own units that others had taken: slower than without the board on 18,838 subjects of
  100 bp at 16 threads); a thread of the subject loop without a subject takes items until the batch's
  units are all taken, and sleeps (a condition variable, 1 ms at most) while nothing is left to take and a
  subject has not yet put its units on the board; a chain's thread does the same with its groups and sorts,
  with `rayon::spawn` helper jobs that take one group of any set each and are spawned one after the other
  (at most one queued per board). An
  item makes no pool call and waits for nothing, so a thread that waits for the items others took blocks
  outside the pool, and no subject's chain runs nested inside an item (the caller of a batch's subject
  loop still waits for its broadcast in the pool, as before): before, the chain of a batch that was not
  the last one (`LOSAT_X_TBXPIPE`'s pipe chains, which link without the pool) linked its groups one after
  the other on one thread, and on AP027131+AP027132+AP027133 x AP027133+AP027078+AP027131 at 16 threads
  the search ended with batch 0's self pair (AP027131 x AP027131: groups of 165,637, 161,243, 131,690 and
  122,349 HSPs) linking on one thread after every unit had run. Its units run the same function on the
  same unit and are kept by unit index, groups are disjoint slices whose kernels read their group only, a
  stable sort has one order (equal chunks sorted with `sort_by`, then two-way merges that take the left
  element on ties), and frames are put together in frame order as before, so every list is the same. With
  `LOSAT_X_TBXPIPE`, the next batch starts when every unit of the batch before has been taken (not at its
  seed-stage end) while the HSPs held by the batches past their seed stage are below TBXPIPE's limit,
  builds its tables and waits at a gate before its subjects until the batch before has ended its seed
  stage and the held HSPs leave room, or every unit of the batch before has run and the records those
  units returned (a bound on the HSPs that batch keeps) fit beside the held HSPs, or the batch before
  returned, taking board items meanwhile; the pipe loop waits for a batch by taking board items. A wait
  outside the pool also runs this thread's own pool jobs and its part of a `rayon::broadcast`
  (`rayon::yield_local`): a broadcast (a later batch's subject loop) ends only after every thread ran its
  part. wasm32 (LOSAT Web) keeps the previous loop (plan v4 line W).
  Memory: the groups of one linking call are linked at once (as the pool path already did), so the
  kernels' per-HSP arrays (about 203 bytes per HSP, `link_group`) of a self pair's four large groups are
  live together instead of one at a time; and the next batch's tables are built while the batch
  before's last units run, so two batches' tables are held for that time. A frame merge's chunk init
  list and converted list get glibc's mapped capacity from 1 MiB on (`x_kmerge.rs` `x_frame_capacity`):
  any of the search's threads now merges a frame, and a 1-32 MiB list freed in a thread's arena stayed
  resident there (without it the three genomes below peaked at +161 MiB over dev at 16 threads). The
  capacity is on every path, also at one thread and with `=0`: the three genomes at one thread take the
  same time with 14 MiB less memory. On the workstation (node 1, 16 CPUs of one socket, 3 s warm-up per
  run, medians of 8 rounds in rotated order with dev and `=0` in the same run, 2 ms VmHWM, every report
  identical; dev = dev + `TBXFEWGROUPS`, 2026-10-31; the times and peaks below are with / without) the three
  genomes took 4.81 / 5.79 s at 16 threads (x0.831; one round at 32: 3.47 / 4.75 s), 8.69 / 9.33 s at 8 and,
  on one node-0 core, 70.5 / 70.0 s at 1 with 1.00008 x the instructions; the same queries as a database
  (`-db`) 4.82 / 5.82 s at 16; the three queries x AP027133 4.54 / 8.23 s at 16; the 18,838 pieces 4.82 /
  5.03 s at 16 and 9.14 / 9.84 s at 8; MelaMJNV x 9 MJNV genomes 1.402 / 1.411 s at 8 and 2.500 / 2.538 s at
  4; AP027131 x AP027133 1.175 / 1.182 s at 16; the Autographivirales set of ViPTreeGen (BLOSUM45, 40
  genomes, `-db`) 6.70 / 6.65 s at 8 (x1.006, the one row above dev: x1.013 in a second run with
  1.003 x the instructions, not found in a profile of the stages); Begomovirus.40 0.29 / 0.28 s at 8
  (runs spread +-15 %). Peak memory (MiB): the three genomes 578 / 536 at 16 (one round at 32: 777 / 657),
  420 / 473 at 8 and 158 / 172 at 1; as a database 574 / 562 at 16; the three queries x AP027133 248 / 237 at
  16; the 18,838 pieces 136 / 133 at 16 and 117 / 116 at 8; MelaMJNV x 9 genomes 128 / 132 and 98 / 105 at 8
  and 4; AP027131 x AP027133 124 / 124 at 16; Autographivirales 50 / 57 at 8; Begomovirus 27 / 26 at 8.
- `LOSAT_X_TBXRPTPAR=1` (TBLASTX, speed; in the strict set from the pull request that added it,
  2026-10-19, on by default; acts only with two or more threads on a pool (native, or Wasm with threads))
  the report stage after the last query batch (`write_tblastx_outputs`) runs its per-HSP work on the
  search's pool (`tblastx/x_rpt_par.rs`, exactness argument in its header): the displayed identities
  (each HSP's codons translated again), the outfmt 0 rows (translated again, masks lowercased; also the
  Web adapter's hit records), the Seq-aligns' engine counts (outfmt 8 and 12, in chunks of the HSP
  list joined in order) and the outfmt 6 / 7 rows of a query with more than 1,024 rows, formatted
  into one buffer per 1,024 rows by the pool's threads, in windows of 4 chunks per thread that the loop's
  thread writes in order before the next window (warnings, outfmt 7 headers and the last line stay on the
  loop's thread). Each value is the serial function of one HSP or row; the list and row order are the
  serial order. With the Web adapter's format observer the rows stay serial. Without it the whole stage
  ran on one thread while the pool waited: on AP027131+AP027132+AP027133 x AP027133+AP027078+AP027131
  (code 4) at 16 threads (ws, dev d7e13f99) the stage took 1.31 s of 8.6 s with outfmt 6 (identities 0.99 s,
  rows 0.17 s, final order 0.15 s) and 3.37 s of 10.6 s with outfmt 0 (outfmt 0 rows 1.80 s, formatting
  0.41 s), the same at 8 threads. The final order, the outfmt 0 text and the writes stay on the loop's
  thread. On the workstation (16 CPUs of one socket, 3 s warm-up per run, medians of alternating runs,
  every report identical; off = `LOSAT_X_TBXRPTPAR=0`), the three genomes above took 7.38 s with it and
  8.39 s without at 16 threads with outfmt 6 (stage 0.25 / 1.31 s), 7.88 / 10.48 s with outfmt 0 (stage
  0.78 / 3.39 s), 12.59 / 13.65 s at 8 threads and 23.42 / 24.05 s at 4; MelaMJNV x 9 MJNV genomes at 16
  threads 1.16 / 1.44 s (outfmt 6) and 1.29 / 1.84 s (outfmt 0); AP027131 x AP027133 at 8 threads 1.96 /
  2.07 s and 2.05 / 2.35 s, at 1 thread 10.43 / 10.41 s (the serial stage). Peak memory is equal within the
  runs' spread (three genomes at 16 threads 529 / 527 MiB with outfmt 6, 485 / 568 with outfmt 0).
- TBXMASKSEEK (graduated 2026-10-27, in the set from the pull request that added it, 2026-10-24; no variable; always
  on; the scan of every range stays as the reference of the unit test) (TBLASTX, the outfmt 0 rows and the hit
  records of LOSAT Web; speed) `PairwiseRows` keeps the shown SEG masks of each
  query frame sorted by start, with the running maximum of their ends (`FrameMasks`, `tblastx/report.rs`), and
  `lowercase_query_row` starts the scan of a row at the first range whose end, or the end of a range before it,
  reaches the row (`partition_point`) and stops at the first range that starts after the row. A letter is lowercase
  when a range of the row's frame covers it, in any order of the ranges, so the rows are the same. Before, every row
  scanned all ranges of its frame: in LOSAT Web `PairwiseRows::hit` took 334 M of the 21.6 G cycles of TBLASTX
  AP027131 x AP027133 (genetic code 4), 33 M with the switch. With outfmt 6 no row is made and the switch does not
  act. Workstation, the same binary with `=0` and without, outfmt 0, alternating order (2026-10-24, outputs
  identical in every run): AP027131 x AP027133 at one thread x0.983 of the time (cycles x0.982, instructions
  x0.976; 8 rounds), MelaMJNV x PemoMJNVA x0.993 (cycles x0.993); eight threads (eight cores, siblings idle, 16
  rounds) x0.990 and x0.986 (cycles x0.981 and x0.995). LOSAT Web (Node, serial reactor, 8 rounds in one process,
  the same reactor with `=0`): the whole search x0.987 and x0.992 of the cycles (instructions x0.969 and x0.995).
- LINKSLICES (graduated 2026-10-22, in the set since 2026-10-15; no variable, no `LOSAT_X_LINKSLICESSHADOW`;
  always on) (TBLASTX, memory) the linking call of the in-place HSP list path links its groups as disjoint
  sub-slices of the one sorted list (`tblastx/x_link_slices.rs`, exactness argument in its header). The
  earlier form cut the sorted list into one Vec per group, linked each and copied the linked groups into a
  new list of the whole size, so the list existed twice while the new one was filled (AP027131 x AP027133,
  genetic code 4: 49 MiB twice in the first linking call, where the heap of that search peaks). The groups
  are the same runs (the same start test, `group_starts`), each kernel writes fields of its own records in
  place and changes no order, and the list after the kernels is the reference's concatenated list, record
  for record; no second list is made. The kernels, their arguments, the per-thread scratch and the pool
  condition (`allow_parallel` and more than one group) are the reference's; `LOSAT_LINKING_LEGACY` and the
  linking traces link each group's copy with the literal port and write it back, since that port returns its
  group in the same order and length. The cut, `link_one_group`, the assembly, the shadow and its
  `[X_SHADOW] LINKSLICES` line went; the reference linking body of `LOSAT_DEBUG_CUTOFFS` /
  `LOSAT_DIAGNOSTICS` keeps its push grouping. `LOSAT_TIMING=1`: `assemble` of the `x_link1` / `x_link2` lines
  is about zero. `LOSAT_X_STATS=1` prints
  `[X_STATS] LINKSLICES linkings= groups= parallel_linkings= legacy_fallbacks=` (linking calls
  linked as slices, their groups, those on the pool, and calls handed back to a cut, which never
  happens: always 0). Measured once on 2026-10-13 (one binary, a loaded machine): massif of AP027131 x
  AP027133 at one thread, heap requested, 107.3 MiB at the peak in the first linking call with the cut;
  with the slices that call reaches 77.3 MiB and the search peaks at 102.4 MiB at the end of the seed
  stage, as with the cut. Max RSS (medians of three alternating runs, cut -> slices): that pair 110 ->
  111 MiB at one thread and 162 -> 163 at eight (its RSS peak is in the seed stage); NZ_CP006932
  against itself 142 -> 132 at one thread and 242 -> 222 at eight. `groups` 0.011 -> 0.002 s and
  `assemble` 0.010 -> 0.000 s in the first linking call of the pair at one thread.
- `LOSAT_X_BPQBLOCK=<letters>` (BLASTP, memory and speed; in the strict set since 2026-10-14, a value switch:
  a positive number of letters, built-in 50,000; `=0` turns it off) the queries are searched in blocks of whole query batches
  (`blastp/x_query_blocks.rs`, exactness argument in its header): a block closes at the first batch end where
  its letters reach the value, and each block is a search of its own (its lookup table, its diagonal table,
  its hit lists) whose rows are written before the next block is searched, as NCBI searches and reports one
  batch at a time. Both lookup tables take it: the default one, and the compressed one of `-task blastp-fast`
  (word size 5), whose table each block builds for its own batches and tests for NCBI's overflow banks with
  them. Taken only where a block's report is its rows: one output of outfmt 0, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15 or
  16, or a `run_local` call with several outputs, hit records or an observer (no `-db`; no outfmt 11 or 13 to 16 in
  such a call), no `-query_loc`, every batch read, more than one block, no subject whose tabular ID the reference path rejects, and none of the
  reference path's other refusals that come before its rows (a `-window_size` whose diagonal table with all
  queries passes 2^30, a query split that LOSAT rejects, or with the compressed table a query batch whose
  offsets add more than four times NCBI's overflow cells, about 859 million words, counted before the search);
  otherwise the search of all queries. With `-db` the same, and the blocks run only when every definition line
  of the database reads and none of LOSAT's refusals of a reported title applies (`x_titles_allow_blocks`, one
  pass over the header files in the search pool before the first block); a query whose effective search space
  could pass NCBI's `Int8` (with `-db` or `-dbsize`) or an empty database also takes the search of all queries.
  Refusals by format, all checked before the first block with the writers' own functions, so that a block never
  fails after the rows of the blocks before; where one applies, the search of all queries runs and refuses (or
  completes) exactly as it does:
  outfmt 0 needs every shown subject title writable (not a title of punctuation that NCBI's `x_CleanAndCompress`
  reads past the end of; not ABI v1's titles with `-db`); outfmt 8, 9, 11 and 12 with `-db` need every definition
  line to read, to be one line and to hold no byte outside ' ' ..= '~' (`database::x_titles_read_allow_blocks`, one
  pass in the search pool; with `-subject` no title is read).
  Outfmt 5 needs the database name and every query's definition line and message id writable and every `-subject`
  record's names (conservative: the writer refuses only a subject with hits), and with `-db` the definition lines
  as for outfmt 8; the fallbacks include a refusable `-subject` title on a subject without hits, a query id that
  would fail in a message never written, and a `-db` sequence with an unreadable header or more than one
  definition line. Outfmt 13 to 16 add the XML's checks; outfmt 15 and 16 with `-db` also need the Blast-def-line
  taxid of every sequence to read and a taxonomy database that cannot refuse a lookup
  (`database::x_titles_taxids_allow_blocks`, `common::database::taxdb_may_refuse`; the pre-check costs about +1.5 s
  at 1 thread and +0.4 s at 4 per 2 M sequences, user time, all of it the pass over the headers); outfmt 13 and 14
  write a side file per query and the master with the last block (an unopenable side file n stops the run at
  report n as before, rc 1; a write error on it is rc 6, DECISIONS 35); outfmt 11 (one BLAST archive document per
  NCBI batch) and 13 to 16 are taken only in a call of one output without hit records and observer (Owner-delegated
  2026-10-22 and 2026-10-23: their run-level frame and carry were verified for that shape only).
  A `run_local` call needs every output to be a writer of the caller, no `-db`, and an observer only with the hit
  records; the title and name checks of every kind of report among the outputs run before the first block, the
  records are handed over once after the last block (each `q_idx` moved by the block's first query), the observer's
  HSP indices are moved past the blocks before, and a run that fails in a later block hands over no records. A
  `-subject` search takes the block path only when no subject's tabular ID fails; the IDs and titles are then made
  once for all blocks (`blast_engine/x_subject_meta.rs`, about 80 B per subject, freed when the block search ends).
  Evidence of the block search against the search of all queries (every report identical; ws, 2026-10-19 to
  2026-10-23, medians; the outfmt 15 / 16 and 11 rows local, WSL): sp4 3000 x Species6 (`-subject`) outfmt 0 at 16 threads 4.43 against 6.55 s, peak 36-58
  against 217-350 MiB; the yeast proteome against its own `-db` outfmt 0 24.63 against 34.40 s, peak 36-58 against
  662-1,266 MiB; bp_hsap `-db` outfmt 8 / 9 / 12 at 1 / 16 threads 32 / 48.6 MiB against 212.8 / 351.9 (NCBI 74-76 /
  131-151), wall x0.81-0.83 at 1 thread and x0.64-0.66 at 16; 3000 x 6,600-protein Species6 `-db` outfmt 15 / 16 at
  1 / 4 threads 25.7 / 37.6 MiB against 219.7 / 259.9, user time x0.87; outfmt 11 (2,200 queries of 120 letters, one
  thread) 23.2 against 53.9 MiB (NCBI 71.2); in `run_local` (the Web build's search, one Node process, 8 alternating
  rounds) NZ_CP006932's proteome (199,959 letters, four blocks) against itself 24.42 against 24.60 s at 1 thread
  and 6.51 against 6.56 s at 4; against 89 SicyWSV proteins x0.931-0.939 at 4 threads. A subject set made once for
  all blocks: human x yeast blastp-fast `-subject` (254 blocks of 50,000 letters) at 16 threads x0.952 of the
  per-block form. The per-format switches `LOSAT_X_BPQBLOCKFMT`, `SA`, `XML`, `XML2`, `WEB` and `LOSAT_X_BPSUBJMETA`
  were folded into this one on 2026-10-31 (structural line S3): they only cut the block search off (`LOSAT_X_BPSUBJMETA` made the IDs and titles once per block instead).
  The output file is created when the first block's rows are written; a run that fails in a later block leaves
  the rows of the blocks before.
  With `LOSAT_TIMING=1` each block prints its own `[TIMING]` lines.
  The block search's set-up uses its thread pool (no value changes): the count of the query batches' words
  for the overflow banks runs on the pool before the first block (groups of consecutive batches, one counter
  per group; when it fails nothing has been written and the search of all queries runs); the subjects are
  encoded once for all blocks (attached to the pool, `x_query_blocks::XBlockSubjects`); and with the
  compressed table and more than one thread, the contexts' Karlin blocks of each block are computed on the
  pool's threads (`tblastx/lookup/backbone.rs` `x_prepare_blosum62_lookup_query_in_pool`: each context's
  block reads only its own letters, `x_context_karlin`).
- `LOSAT_X_BPQSLOTS` (two query blocks of `-task blastp-fast` searched at a time on the pool, `blastp/x_query_slots.rs`) was
  removed on 2026-10-28 (Owner): at 16 threads on the workstation it saved 7.0 % of the time and raised the max RSS by
  9.1 % (136 -> 150 MiB), which breaks the rule that the rise in max RSS in percent must not exceed the time saved in
  percent. A run ignores the variable and searches the blocks of `LOSAT_X_BPQBLOCK` one at a time, as it did with
  `LOSAT_X_BPQSLOTS=0`; the reports are unchanged. The entry and its evidence are in the history of this file
  (dev b1af130b).
- `LOSAT_X_BPSTRIPEFILL=1` (BLASTP `-task blastp-fast`, speed; in the strict set from the pull request that added
  it, 2026-10-28, on unless it is `0`; acts wherever `LOSAT_X_BPLKPFLAT`'s one-latch build of the flat table of
  `LOSAT_X_BPCLUT` runs, so with more than one thread) builds each stripe of that table from one read of the query
  chunks' parts (`tblastx/lookup/x_clut_par.rs` `x_stripe_fill`, exactness argument in its header): the stripe
  copies its part of each chunk, in chunk order, while it counts its cells, then lays out and fills from the copy;
  before, `x_stripe_build` read the 128 parts (at 16 threads) twice, once to count and once to fill, and the second
  read of short parts written by other cores stalled (line BPF `REPORT.md` §13-§14: the stripes pass was about 14 %
  of the cycles of the human x yeast search at 16 threads, with a cycles per instruction of 0.565 against 0.395 for
  the serial bucket fill). The overflow cells of NCBI's table are counted only when some table's appends pass
  its capacity (an append takes at most one, so otherwise the capacity test cannot refuse). Memory: each task of
  the pass keeps a copy as large as the largest stripe it took (4 bytes an append) beside the stripe's counts.
  At one thread nothing changes (the serial builder runs). A form that wrote the chunks'
  appends straight into the table's hits instead (line BPF s10, `exp/bpf-s8` 3c405fdb) was as fast and took
  20-35 MiB more (its scattered stores stalled in turn).
  Measured 2026-10-28 on the workstation (node 1 CPUs 32-47 / 32-39, siblings idle, `--membind=1`, CPUs warmed 3 s
  before every run, 8 rounds in alternating order, dev 4403cc12 against 1cbb5109, one output hash per case and
  thread count; line BPF `results/s10ab`), `-task blastp-fast`, wall at 16 / 8 threads: human x yeast (Species4 x
  Species6 `-subject`) 4.550 -> 4.447 s (x0.977) / 8.270 -> 8.044 s (x0.973); its first 3,000 proteins 0.458 ->
  0.448 (x0.978) / 0.778 -> 0.758 (x0.974); Species4 x Species6 `-db` 4.723 -> 4.607 (x0.976) / 8.389 -> 8.157
  (x0.972); Species0 x Species3 `-db` 14.204 -> 14.161 (x0.997) / 27.049 -> 26.852 (x0.993). Cycles x0.964-0.970
  (x0.993 on Species0 x Species3), instructions x0.989-0.991; the median of the per-round max RSS differences
  -5.2 to +1.5 MiB (no rise); one-thread instructions the same (the serial builder).
- `LOSAT_X_BNTABPAR=1` (BLASTN, all tasks, outfmt 6 / 7, two or more threads, speed; in the strict set since
  2026-10-24 (perf/bn-tab-par, evidence at the end of this entry): on unless it is `0`) the rows of a report are
  formatted on the search's pool: they are cut into pieces of at most 512 rows of one query (the first piece of a
  query with its outfmt 7 header; a query without rows is one piece), the pieces are taken in windows of at least
  2,048 rows per thread (a piece without rows counts as one), each piece of a window is formatted into a buffer of
  its own by `write_query_header` and `write_hsp_rows` (the two parts of the serial writer's function of one query,
  `write_query_rows`, `blastn/hsp.rs`), and the thread of the batch loop writes the window's buffers in order, each
  query's warnings before its first piece, before the next window is formatted; outfmt 7's last line stays on that
  thread (`blastn/x_tab_par.rs`, called from `post_process_hits_and_write` in `blastn/blast_engine/run.rs`). Same
  bytes: a query's pieces in order are its header and its HSP lists' rows in order, which is the serial call for
  that query, and the serial writer is those calls in query order. Where: the search's pool, which `write_reports` gets in the batch loop
  and after it; serial with one thread, with the Web adapter's format observer and with the BLASTN trace.
  Memory: the rows of one window (about 6 MiB at 32 threads with ~100-byte rows), also for one query of many rows
  or many queries without rows. `LOSAT_X_REACH` records the decision once per report call of at least one query
  (`switch_reach_expected.tsv`: conditional). Tests: `x_tab_par` unit tests (random hit lists, outfmt 6 and 7, with
  and without the epilog, warnings, 2 / 3 / 8 threads, a query larger than a window: the serial writer's report,
  flushes and warnings in the same order; every row in exactly one piece, in order; a format observer keeps the
  serial writer), integration test `blastn_tab_par` (default vs `=0` reports equal, the new side reached).
  Premise (2026-10-24, v4-D s12): 31,168 cDNAs x 3 Mb `-task blastn` `-outfmt 6` (1.02 M rows) with
  `LOSAT_X_QBATCHPAR=0` (the batch loop, now the only form) on ws spent 0.88 s in the report calls at 8, 16 and 32 threads alike (`[TIMING] common`
  output 0.879 / 0.886 / 0.884 s), on one thread while the pool waited; local perf: 87 % of the call in
  `write_hit_fields_bytes`. Evidence (2026-10-24, ws node 1, ratios per round, median): window v4ds13-w2, this
  code (216ede53) vs dev 499cd95d, 106 runs, one sha per case: AB with `LOSAT_X_QBATCHPAR=0` x0.929 / 0.882 /
  0.866 at 8 / 16 / 32 threads (10.83 -> 10.12, 7.18 -> 6.33, 5.70 -> 4.93 s; `=0` 0.998 / 0.985 / 0.993 of dev),
  C5 (5,000 cDNAs) loop n16 x0.884, `-db` (C5's queries x the 3 Mb database) n16 x0.887, A1 (1,775 cDNAs, blastn,
  by default) n16 x0.909, AB by default (slots, removed on 2026-10-28: serial) 1.000; cycles and instructions per run equal (0.996-1.011);
  max RSS AB loop n32 124 -> 127 MiB; `[TIMING] common` output AB loop 0.904 / 0.936 / 0.947 -> 0.241 / 0.230 /
  0.219 s. Window v4ds13-w1 (the first version, one buffer per query, 10f8aa9e vs dev 7a731a15, 370 runs): the
  same AB / C5 loop ratios (x0.932 / 0.892 / 0.869, x0.917 / 0.892 / 0.877), A1 by default x0.933 / 0.904 / 0.862,
  M1 / D1 by default 0.96-1.01, `-db` loop n16 x0.904, and the defaults with slots (G2, MB, C5, AD, AB at 8 / 16 /
  32 threads) 0.96-1.01.
- `LOSAT_X_PAIRWISEPAR=1` (BLASTN, TBLASTX, BLASTP, speed; in the strict set from the pull request that added it,
  2026-10-19, on by default; acts only with two or more threads on a pool (native, or Wasm with threads)) the pairwise
  report (outfmt 0) of BLASTN (every task, `-subject` and `-db`), TBLASTX and BLASTP is formatted on the search's pool
  (`report/x_pairwise_par.rs`, exactness argument in its header). The three writers of `report/pairwise.rs` have one
  shape: after the prolog, per query its warnings, its head (header, then "No hits found" and the footer, or the
  description table), the alignments of the first `-num_alignments` subjects and its footer; the epilog after the last.
  With the switch the queries are taken in windows of about 1 MiB of estimated text per thread; the pool groups each
  window's hits by subject and formats each query's head, runs of about 256 KiB of its alignments (cut between HSPs,
  also within one subject's) and its footer into one buffer each, in chunks of the same estimated size; the loop's
  thread writes a chunk's buffers in the serial order, each query's warnings before its head, before the next chunk is
  formatted. BLASTN's and BLASTP's outfmt 0 rows (`pairwise_hits`, `build_pairwise_hits`) are made on the pool too, in
  the list order (the first error in that order). Each piece is the serial writer's calls on the same arguments, into a
  buffer instead of the writer; a piece that fails has its bytes up to the failure written, then the error. The serial
  writer stays: `LOSAT_X_PAIRWISEPAR=0`, one thread, a format observer (the Web adapter's; text only). Without it the writers ran on the loop's thread while the pool waited: megablast 100 queries x 16S rRNA
  (`-subject`, 12 NCBI query batches, `LOSAT_X_RPTSTREAM` writes each batch's reports after its search) spent 1.4-2.0 s
  of 7-10 s there at 8 and 16 threads, TBLASTX's outfmt 0 text 0.41 s of the three genomes below at 16 threads. On the
  workstation (16 CPUs of one socket, 3 s warm-up per run, medians of alternating runs, every report identical; off =
  `LOSAT_X_PAIRWISEPAR=0`): megablast 100 queries x 16S 4.96 / 6.41 s at 16 threads, 9.36 / 10.70 s at 8, 67.46 /
  67.45 s at 1; TBLASTX AP027131+AP027132+AP027133 x AP027133+AP027078+AP027131 (code 4) 7.67 / 7.99 s at 16 threads,
  12.95 / 13.23 s at 8; MelaMJNV x 9 MJNV genomes 1.23 / 1.30 s at 16; AP027131 x AP027133 2.01 / 2.07 s at 8;
  BLASTP 3,000 proteins x Species6 (`-subject`) 4.59 / 4.85 s at 16 threads and 7.71 / 7.89 s at 8, the
  yeast proteome against its own `-db` 23.13 / 25.72 s at 16 and 37.32 / 39.59 s at 8 (in `LOSAT_X_BPQBLOCK`'s
  blocks), equal at 1; BLASTN cDNA 31,168 x 3 Mb, where `LOSAT_X_QBATCHPAR`'s slots ran (removed on 2026-10-28; the switch did not act), 6.26 / 6.24 s.
  Peak memory: a window's text is at most 1 MiB of estimate per thread plus one piece; megablast 85-87 / 89-94 MiB at
  16 threads and 91-93 / 78-81 MiB at 8 (the pool threads' glibc arenas keep the freed buffers), TBLASTX and BLASTP within the
  runs' spread (BLASTP 40-51 / 42-53 MiB).
- BPPFREE (graduated 2026-10-26, in the set since 2026-10-14; no variable; always on) (BLASTP, memory) the
  preliminary hit lists hold no more than their HSPs need
  (`blastp/x_prelim_free.rs`, exactness argument in its header): a subject's HSP list for a query gives back
  the room of the HSPs that the preliminary E-value reaped, and the composition-based redo frees each
  query's preliminary list after its redo (the serial redo query by query, the windows of `x_redo_par.rs`
  after the replay of the window that holds the list's last match), as NCBI's redo frees each list it has
  redone. The match-parallel and query-parallel redo paths keep the lists until they end. No value changes, so
  there was no shadow. With
  `LOSAT_TIMING=1` each search of a query set prints `[TIMING] x_bppfree: freed_lists=<n>`. Measured once
  (C. elegans block of 150 proteins, 50,869 letters, x human proteome, one thread, max RSS): 153,744 ->
  117,968 kB (NCBI -db 109,080 in the same run).
- BPPCOMPACT (graduated 2026-10-26, in the set since 2026-10-14; no variable, no `LOSAT_X_BPPCOMPACT=shadow`,
  no `[X_SHADOW] BPPCOMPACT` line; always on) (BLASTP, memory) the subject
  scan keeps each preliminary HSP as a 64-byte element (the
  `BlastpPreliminaryHsp` and its E-value) instead of a 168-byte `BlastpHsp`, in the same hit-list helpers
  made generic over the element (`blastp/x_prelim_compact.rs`, exactness argument in its header); the redo
  reads lists of full HSPs rebuilt from them: the one-thread redo rebuilds each query's list when it
  reaches it and drops it after its redo, a search with more than one thread rebuilds all lists before the
  redo (or, where `LOSAT_X_BPREDOPAR` runs, window by window), query chunks hand on rebuilt lists; with `-db`
  the compact lists are listed and renamed to ranks (`ListedSubjects::read`). The unit test
  `hit_list_helpers_keep_the_same_lists_on_both_types` compares the hit-list helpers on both element types.
  With `LOSAT_TIMING=1` each search of a query
  set prints `[TIMING] x_bppcompact: expanded_lists=<n>`.
- `LOSAT_X_BNSTREAM=1` (BLASTN, all tasks, `-subject` with two or more threads, memory; in the strict set
  since 2026-10-16, on by default) the threads claim the subjects
  of a query batch one at a time in OID order and one of them (the owner) hands the results to the hit lists in
  OID order while the others go on with the next subjects; a thread claims a subject only while it is fewer
  than `threads x 64` subjects past the last one handed over, so at most that many subjects' results wait
  (`x_bnstream::search_in_order`, called from `search_query_batch` in `blastn/blast_engine/run.rs`; each
  thread runs its loop as one `rayon::broadcast` job, so that a thread waiting inside a nested parallel stage
  of a search never starts a second loop, and since 2026-10-23 no subject is searched before every thread
  has started its loop: rayon runs a thread's own broadcast job inside a wait for a nested job too, so a
  thread that had taken a nested job of a search before its loop started ran the loop above that job and the
  pass hung, rarely, on a loaded host; the ignored probe `probe_b17_many_passes_with_a_nested_first_item` of
  `x_bnstream.rs`, 20,000 passes on more threads than CPUs, hung in 12 of 12 runs before and in none of 12
  after), in both passes: the preliminary HSPs to the preliminary hit lists
  (`PrelimHitCollector`, which keeps those of at most `prelim_hitlist_size` subjects per query) and the
  traceback's HSPs to the final hit lists (`update_hitlists_with_subject_hits`). Without it two or more threads
  search every subject of the batch and then hand all results over, so the HSPs of every subject are held at
  once. The hit lists take the subjects one at a time in OID order either way, as NCBI's
  `Blast_HitListUpdate` takes one HSP list at a time, so the report is the same. One thread always hands each subject over as soon as it is searched (no
  switch). No value changes, so there is no shadow; `LOSAT/tests/blastn_subject_stream.rs` compares the
  reports of one thread and of two and four threads with and without the switch.
  Measured 2026-10-14 on the WSL machine (megablast, 100 16S queries against 28,036 16S subjects, 8 threads,
  output unchanged): max RSS with the blastna and packing changes of #187 and #188 applied, one run each in one
  series, 207.6 MiB without the switch and 145.3 MiB with it (NCBI BLAST+ 2.17.0 `-db` 115.6 MiB); time on this
  branch, 8 alternating pairs (load 5-19), a median of 8.09 s without and 8.91 s with it (pair ratio median 1.13,
  1.04-1.25), and 0.39 / 0.47 s with one query. Blocks of 256 or 1,024 subjects per thread were not faster within
  the noise (6 rounds each: ratios 1.23 and 1.14 against 1.19 for 64), and 1,024 kept 160 MiB. That form
  searched blocks of `threads x 64` subjects and handed a block over before the next started; on the
  workstation (2026-10-15, 8 pairs) it took 1.11 times the wall time, the threads idle at ~110 block ends of
  ~10 ms each, which the window above removes. The window, on the same workstation (2026-10-16, node 0
  pinned, rounds of four arms in rotating order, output unchanged): 8 threads, 8 rounds, wall time 0.944 of the
  default (pair ratios 0.900-0.957; the block form 1.103), max RSS 202 -> 141 MiB; 16 threads, 4 rounds, 0.882
  (0.834-0.890; block form 1.145), 211 -> 145 MiB; one query at 8 threads 0.998, 79 -> 67 MiB. It is faster than
  the default, probably because the hand-over runs while the threads search instead of after them and fewer
  results are held (not measured separately).
- `LOSAT_X_BNSTORE=1` (BLASTN, all tasks, `-subject` without `-subject_loc`, memory; in the strict set since
  2026-10-18 (perf/bn-memory-rules, evidence at the end of this entry): on unless it is `0`)
  each subject record's letters are packed as
  soon as the reader has read it (`sequence/packed_subject.rs`, the module of `LOSAT_X_TBNSTORE`: ncbi2na 4 bases
  per byte with the remainder count in the last byte, the runs of letters other than `ACGT` as ncbi4na codes, the
  lower-case intervals), so the letters of the whole set are never held; the records keep their deflines. The
  search reads each form from `SubjectLetters` (`blastn/x_subject_store.rs`): the length, the lower-case intervals
  (`-lcase_masking`), the packed bases of the preliminary stage (those of the packed subject, not a second copy;
  the bytes of `encode_subject_ncbi2na_packed`), the blastna letters of a subject in the traceback
  (`blastna_into`) and the upper-case letters of each outfmt 0 row (`letters_into`), each the bytes the record's
  letters give (the comparison tests of `packed_subject.rs`). A subject range keeps the letters (the records
  are cut before the search), and so does a set with a byte that has no ncbi4na code (the CLI reader admits
  none); the library entry `run_local` (and LOSAT Web) reads the caller's records. With `LOSAT_X_STATS=1` a
  `[X_STATS] BNSTORE subjects= stored_bytes= blastna_made= rows_made=` line gives the subjects held packed, the
  bytes of their packed forms, the times the traceback made blastna letters and the report rows made from them;
  `LOSAT/tests/blastn_subject_store.rs` compares reports and messages with and without the switch and reads it.
  Measured 2026-10-15 on the WSL machine (megablast against 28,036 16S subjects, 40.6 Mb; max RSS, one run each
  in one series, the same output in every run): one thread, one query 68.0 -> 30.5 MiB (NCBI BLAST+ 2.17.0 `-db`
  42.0), a query without hits 67.2 -> 29.4 (34.6), 100 queries 135.3 -> 97.7 (95.8); 8 threads, one query 81.1 ->
  44.7, 34.8 with `LOSAT_X_BNSTREAM` (43.8), 100 queries 206.7 -> 170.8, 108.9 with `LOSAT_X_BNSTREAM` (145.0 with it
  alone; NCBI 115.0). One 5.5 Mb subject (Sakai x MG1655, 8 threads) keeps its peak (196.2 -> 194.3; NCBI 152.8):
  there the query's copies dominate, not the subject.
  Time: 100 queries at one thread, 4 alternating pairs (load 6-12), pair ratio median 1.06 (0.87-1.11).
  Stays opt-in (2026-10-17, perf/perfmb-default-bundle). Workstation pairs of 2026-10-18 (perf-mb
  `ws/defbundle_plan2.txt`, the `+BNSTORE` arm against that bundle's new default, 16 physical cores of one node,
  rotated rounds; NCBI BLAST+ 2.17.0 `-db` in brackets): megablast of 16S, one query, at 1 / 8 threads time 0.998 /
  1.000 and max RSS 64 -> 27 / 64 -> 26 MiB [86 / 82]; 100 queries 1.000 / 1.000, 86 -> 47 / 90 -> 54 MiB [140 /
  146]; Sakai x MG1655 at 8 threads 0.997, 145 -> 141 MiB [194]. So with the bundle's other switches (RPTSTREAM's
  streamed reports) the 16S cells are below NCBI `-db` without it too, and it lowers them by a further 37-39 MiB at
  no time cost. On dc-megablast of 1,775 cDNAs against 3 Mb of chromosome 2L it lost time: 1.033 (1.010-1.057 per
  round) at 4 threads and 1.187 (1.098-1.355) at 8, searches of 0.6 s (1.000 at 1 thread, 0.998 at 16); on the other
  dc-megablast and blastn cells its ratios lay on both sides of 1 (0.78-1.16), inside those cells' round spread. So
  it stayed opt-in then (2-3 rounds per cell; no rule separating the cells was tried).
  Default ON since 2026-10-18 (perf/bn-memory-rules; workstation, 16 physical cores of one node with their siblings
  idle, a warm-up and 6-8 rounds per cell with the arms rotated, max RSS of both arms in the same window, cycles:u by
  `perf stat`; v4-D `results/d10_w1_summary.md`, `d10_w2_summary.md`): the dc-megablast loss did not reproduce, 1.006
  (0.935-1.057) at 4 threads and 0.997 (0.917-1.047) at 8 in 8 rounds with equal cycles, so the old figures were the
  cell's round spread (the same old window's 3,492-cDNA cell read 0.78 / 0.88); blastn of the 1,775 cDNAs at 16
  threads 1.014 (0.961-1.041) and, on the candidate, 1.001 (0.835-1.244), cycles within 3 %; megablast of 16S, one
  query at one thread 0.999 at -38.0 MiB (26.0 MiB), 100 queries at 8 threads 0.999 at -38.2 MiB. It saves memory
  wherever subjects are held (-3 to -38 MiB) and never adds any, and no measured cell loses time, so it has no rule.
- SEEDRELEASE (graduated 2026-10-22, in the set since 2026-10-13; no variable; always on) (TBLASTX, memory)
  the tables that only the scan reads are freed before the HSP list of the batch's last subject is merged
  and linked, which is where the heap of a large search peaks
  (`tblastx/x_seed_release.rs`, exactness argument in its header). A query batch keeps its tables in two
  holders, each reached only through one stage per subject: the unit path (PAIRPAR) with its group
  tables, whose stage ends when the subject's units have run, before their lists are merged (the unit path
  runs as two calls for this, `run_subject_units` and `merge_subject` in `x_pair_par.rs`, the same statements
  in the same order); and the batch's lookup table, whose stage ends after the seed stage. The end of the
  batch's last stage of a holder frees it. The worker's diagonal table, which only the reference frame loop
  reads, is not created when the unit path takes the batch, and is otherwise dropped by
  the worker that ended the last lookup stage. The earlier off state, in which the holders kept the tables
  until the batch had been searched, went with the variable. Measured once on AP027131 x AP027133 (genetic
  code 4, one thread, bytes requested under massif), before -> after: the end of the units 110.4 -> 102.4
  MiB, the merge of their lists 155.8 -> 81.0 MiB, the first linking call 167.6 -> 107.3 MiB, which is the
  peak of the run in both; max RSS in one run each 167 -> 110 MiB at one thread and 221 -> 163 MiB at eight.
  Known limits: the gain is for the last subject of a batch (the only one of a pairwise search), since
  earlier subjects are merged and linked while later ones still need the tables; and with
  where PAIRLUT's fallback builds the groups from the full table the batch's full lookup table,
  which the frame loop names, stays alive through the merge, so that path keeps the higher seed-stage
  peak. `LOSAT_X_STATS=1` prints
  `[X_STATS] SEEDRELEASE batches= tables_released= unit_tables_released= batches_without_worker_table= worker_tables_released=`
  (batches searched; batches whose lookup table a stage end released; batches whose unit path a stage end
  released; batches whose workers held no diagonal table; worker tables dropped at the release). There is no
  shadow switch: the release changes no value and no order of computation, only when three allocations are
  freed, so there is no second result to compare. A read after a release is excluded by ownership instead:
  a stage is the only way to its tables and ending it consumes it (the compiler rejects a later use), the
  merge after the units stage is a function that cannot name a group table, an open stage keeps its tables
  alive, a subject that arrives after a release stops the run, and an assertion stops a frame loop that
  would run without the worker's table (`LOSAT/tests/tblastx_seed_release.rs` reads the counters).
- `LOSAT_X_SEEDCOMPACT=1` (TBLASTX, memory; in the strict set since 2026-10-17, on by default) the seed
  stage of the unit path (PAIRPAR) and its merge hold fewer and exactly sized lists
  (`tblastx/x_pair_par.rs`, exactness argument in the section above
  `XUnitHsp`). Without it each unit keeps its saved init HSPs as `InitHSP`s (72 bytes on 64-bit targets) in a
  list grown by doubling, the lists of all units of a subject live until its merge, and the merge converts
  each chunk's list into a list grown by doubling and copies the first chunk's list into the empty frame
  list. With it a unit keeps each init HSP as a 28-byte record, NCBI's own record (`BlastUngappedData` and
  the seed offsets: query and subject starts, length, score, the two seed offsets) and the context index,
  built in a per-thread vector that is kept across units and copied out at its exact length; the merge
  expands the records of a chunk's units, in unit order, into the same `InitHSP`s (the other fields come from
  the unit's subject frame, the context and the subject), converts the chunk's list into a list reserved at
  its length (`get_ungapped_hsp_list_into`, the body of `get_ungapped_hsp_list`), and moves the first
  chunk's list into the empty frame list (what `merge_tblastx_subject_chunk_hits` computes for an empty
  list, without its copy). No value, order or sort changes. With the switch off the statements of the unit
  path run as before, after one read of the mode per subject (`x_hit` and the unit scan are generic over
  the list they save into and compiled once for each). It acts only where the unit path runs (the
  scan loop is untouched). There is no shadow variable of its own: the unit test `compact_check` in
  `x_pair_par.rs` compares the unit path's `combined_ungapped_hits`, made from the compact lists, with
  the one made from the full lists (both are the unit path), and the per-unit hits with
  `masked_reference`; unit path against frame loop is checked only through the frozen NCBI bytes at
  1 thread (frame loop) and 4 threads (units). `LOSAT_X_STATS=1` prints, whenever the
  switch is on, `[X_STATS] SEEDCOMPACT subjects= units= records= bytes_saved=` (subjects and units merged
  from compact lists, records expanded, and records x (size of `InitHSP` - 28): the bytes the unit
  outputs no longer request; all 0 where no unit path runs). `LOSAT/tests/tblastx_seed_compact.rs` reads
  the counters. Measured once on 2026-10-17 (WSL machine under load, one binary, max RSS from `wait4`, means
  of three alternating runs without -> with the switch, every report the frozen one): AP027131 x AP027133
  (genetic code 4) 111.0 -> 88.7 MiB at one thread and 164.0 -> 154.5 at eight; NZ_CP006932 against itself
  132.2 -> 123.9 and 226.3 -> 215.1 (at eight threads the peak is then the first linking call). The medians
  of `units`, `merge` (`[TIMING] x_pairpar`) and `search_total` moved within the run-to-run spread of that
  machine; the probe build measured the same memory on the workstation (-22.6, -7.5, -8.1 and -9.4 MiB).
- `LOSAT_X_SEEDUNMAP=1` (TBLASTX, memory; Linux with glibc only, at two or more threads, and only where
  transparent huge pages are not "always"; in the strict set since 2026-10-17, on by default from the change
  that added it, under the Owner's rule of 2026-10-17 that a speed-up or memory reduction lands on by default)
  three blocks of the unit path (PAIRPAR) on the search pool get a capacity that glibc always maps,
  so that they are returned to the system when they are dropped (`tblastx/x_pair_par.rs`, the section above
  `X_MAPPED_BYTES`): the unit diagonal table and the record vector of `LOSAT_X_SEEDCOMPACT` that each rayon split
  of a subject's units keeps, and the copy of a unit's records that lives until the subject's merge. glibc
  raises its mmap threshold to the size of each mapped block that is freed (up to 32 MiB; SEG's buffer and the
  pair lookup's group tables are such blocks) and its trim threshold to twice that, so without the switch these
  blocks come from the worker threads' arenas, which keep their freed tops through the linking call, where the
  search peaks (10-14 % of the max RSS at 8 threads). With it, for a unit whose query-group span and subject
  chunk both have at least 65,536 letters (`x_unit_mapped`), the split's table and record vector are replaced
  by empty ones with a capacity of 32 MiB and 4 KiB (above glibc's largest dynamic threshold) unless they have
  that much room, and a copy of at least 128 KiB gets that capacity; the pages beyond their length are never
  touched and are not resident. Three rules keep today's capacities where mapping costs more than it returns.
  (1) Both sides long: a mapped table is new pages at every split, about one page fault per 1,024 cells of its
  width (span + chunk), while the unit's scan grows with span x chunk, so the cost per work is about
  1 / span + 1 / chunk. With an earlier rule (the chunk at least a quarter of the table), 200 queries of about
  3 kb from AP027131 against AP027133 (genetic code 4, 8 threads; 50 query batches, group spans of about 4,000
  letters, chunks of 202,000) mapped 1,748 tables: sys time 0.05-0.2 -> 0.4-0.7 s, wall x1.05-1.07; MelaMJNV
  against AP027133 in six pieces (span about 96,000, chunks about 33,000) mapped 159 (sys 0.034 -> 0.070 s); and
  AP027131 against 200 contigs of about 3 kb took x1.07 of the time with every table mapped. Under the two-side
  rule none of these three maps; AP027131 x AP027133 (span about 220,000, chunk 202,000), NZ_CP006932 against
  itself (219,000 / 219,000) and MelaMJNV x PemoMJNVA (96,000 / 98,000) do. (2) Transparent huge pages: where
  the kernel gives anonymous mappings huge pages by default (`/sys/kernel/mm/transparent_hugepage/enabled`
  "[always]": WSL2's kernel and the RHEL family's default), the mapped 32 MiB regions get 2 MiB pages on first
  touch, so blocks of 0.1-0.5 MiB are resident as 2 MiB (WSL2, 8 threads, 5 alternating pairs: the 200 queries
  above 53-56 -> 60-67 MiB of max RSS, AP027131 x AP027133 137-140 -> 124-145, higher in 4 of 5 pairs; the same
  binary with THP disabled for the process by prctl PR_SET_THP_DISABLE: 52 -> 44, 135 -> 122 and NZ_CP006932
  185 -> 168 MiB, lower in 5 of 5). So the switch acts only where that file selects "[madvise]" or "[never]" or
  `/proc/self/status` says `THP_enabled: 0`, read once per process; on machines with THP "always" (WSL2, the
  RHEL family's default) or where the file cannot be read it does nothing. (3) One thread: the units then run
  on the calling thread, whose blocks come from glibc's main arena, which glibc trims, so nothing is retained
  there and mapping only costs (AP027131 x AP027133: 77 -> 78 MiB and 2,800 more page faults when mapped); that
  serial path keeps today's capacities. No value, order or sort changes: the table is still cleared and resized
  to the unit's width, the record vector is still cleared for every unit and filled by the same pushes, and the
  copy holds the records in their order; the merge reads the lengths and the elements of the copies, never
  their capacity. With `LOSAT_X_SEEDCOMPACT=0` only the table acts; with the switch off, at one thread, with THP
  "always" and off Linux with glibc, every capacity is as before. There is no shadow variable: only capacities
  change, so there is no second result to compare. `LOSAT_X_STATS=1` prints, whenever the switch is on,
  `[X_STATS] SEEDUNMAP active= diag= records= copies=` (`active` 1 where the switch acts in the process: Linux
  with glibc and the THP condition; then the blocks given the mapped capacity at the three places, all 0 at one
  thread, where no unit path runs, where the rule never holds, and where the switch does not act).
  `LOSAT/tests/tblastx_seed_unmap.rs` reads the line. On the workstation (THP "[madvise]", so the switch acts;
  8 rounds alternating with `LOSAT_X_SEEDUNMAP=0`, one warm-up per block, every report the frozen one or the
  case's first run), the max RSS at 8 threads went from 128 to 118 MiB on AP027131 x AP027133, from 176 to 163 MiB
  on NZ_CP006932 against itself and from 41 to 37 MiB on MelaMJNV x PemoMJNVA (lower in 8 of 8 rounds each), with
  the time x1.006, x1.000 and x1.001 (median of the per-round ratios); MelaMJNV against AP027133 in six pieces
  (44 MiB), the 200 queries (46 MiB) and AP027131 against 200 contigs (84 MiB) keep the arena and did not change;
  nothing changed at one thread.

- TBXFRAMES (graduated 2026-10-26, in the set since 2026-10-17; no variable; always on) (TBLASTX, memory) a
  subject's frames are made when they are read and only as long
  as they are read (`Subjects::x_preliminary_frames` / `x_letter_frames` in `tblastx/subjects.rs`,
  `x_generate_frames_in_windows` in `tblastx/translation.rs`). Before, `process_subject` made, before
  the seed stage, the six frames of the subject's letters and, for a subject with ambiguous bases, the six
  frames of the bases the preliminary stage translates, each from a whole copy of its sequence and of its
  reverse complement, and kept both sets to the end of the subject (a stored sequence with ambiguity codes:
  about 7 bytes per base at once); that code went with the variable, and `Subjects::with_search_bases` and
  `preliminary_subject_bases`, which made those whole sequences, stay as the reference of the unit tests
  (`cfg(test)`). Now only the preliminary frames are made before the seed stage, in
  one pass over windows of 196,608 letters (a record's ambiguous letters resolved window by window from the
  one `CRandom` of its length, in letter order, which draws what `preliminary_subject_bases` draws; a stored
  sequence's 2-bit codes read directly) through a table of the residues of the 4,096 codons of ncbi4na masks
  (`GeneticCode::get` reads a codon only through the masks of its letters); they are freed after the seed
  stage, before the first linking call; and a subject with ambiguous bases gets the frames of its letters
  (the re-evaluation's) only when HSPs are left for the re-evaluation, as NCBI translates the ncbi4na
  subject there (blast_hits.c:2650-2653). The seed stage's HSP lists and the first linking call read only
  the number and lengths of the frames, which both sets share, so they take the preliminary frames. Same
  frames by construction (unit tests compare every frame with `generate_frames`, and the methods with
  `with_search_bases`, for records and stored sequences with runs); `LOSAT/tests/tblastx_subject_next.rs`
  compares whole reports. Measured on the workstation (node 0, CPUs 8-15, one thread, medians of 3
  alternating pairs, reports identical): 3 fruit fly transcripts x the fruit fly genome database (1,870
  sequences, N runs) 11.59 -> 9.75 s (NCBI `-db` 10.05 s), max RSS 257 -> 107 MiB (NCBI 162 MiB); the same
  with the genome FASTA as `-subject` 12.96 -> 11.39 s, 328 -> 217 MiB; 10 transcripts x the zebrafish cDNA
  database 62.74 -> 61.62 s, 46 -> 45 MiB.
- `LOSAT_X_TBXRESOLVE=1` (TBLASTX, speed; in the strict set from the pull request that added it, 2026-10-24;
  `=0` turns it off) the preliminary frames of a `-subject` (or `run_local`) record with a letter other than A,
  C, G and T (`Subjects::x_preliminary_frames`) read the base of each letter of one
  base (`A`, `C`, `G`, `T`, `U`, either case) from a table (`tblastx/subjects.rs`, `ONE_BASE`), and only the other
  letters go through `ncbi4na_of` and `resolve_ncbi4na_base`. Before, every letter of such a record went through
  both. NCBI draws from `CRandom` only for an ncbi4na value outside `ctable` (not 1, 2, 4 or 8;
  blast_objmgr_tools.cpp:441-450), so the draws and the bases are those of the letter-by-letter walk (a lib test
  compares every cell of the table with the resolve and every frame of both walks). Conditional: a record of A, C,
  G and T only copies its letters, and a `-db` sequence is read from its stored codes (`x_reach!` per record).
  Found by the gate-gap audit's `-subject` packing item (line G, 2026-10-24): with the fruit fly genome FASTA as
  `-subject`, `translate` took 2.31 s against 0.91 s for the same genome as `-db` (`LOSAT_TIMING`, 3 transcripts,
  one thread), and `resolve_ncbi4na_base` plus its caller were 15 % of the cycles. Measured on the workstation
  (node 0 CPUs 16-23, HT siblings idle, cores warmed 3 s, one warm-up round + 3 rounds, arm order rotating,
  medians of the per-round ratios against the dev binary before the switch, bf92a834; reports and max RSS
  identical; load 14-30): 3 fruit fly transcripts x the fruit fly genome FASTA (`-subject`, 1,870 records) 1 thread
  9.72 -> 8.31 s (cycles x0.855, instructions x0.948), 8 threads 1.52 -> 1.28 s (cycles x0.835); the genome as
  `-db` 8.22 -> 8.15 s (cycles x0.992, instructions x1.000). With `=0`: cycles x1.000-1.001 of the dev binary.
  WebAssembly (serial WASI command, one Node process, CPU 20, 3 rounds, against `=0`): the fruit fly genome FASTA
  as `-subject` 13.17 -> 11.38 s (x0.864), linear memory equal (206 MiB). Since `LOSAT_X_TBXSTORE`
  (2026-10-25) the command line searches its records packed: the switch acts on `run_local`'s records (LOSAT Web,
  the library) and with `LOSAT_X_TBXSTORE=0`.
- `LOSAT_X_TBXSTORE=1` (TBLASTX, memory; in the strict set from the pull request that added it, 2026-10-25;
  `=0` turns it off) the records of `-subject` are packed as they are read (`blastinput/x_packed_subjects.rs`,
  `x_read_packed_subjects`, shared with `LOSAT_X_TBNSTORE`: a record is cut to `-subject_loc` first; at two or
  more threads a second thread packs while the reader reads, natively), the records keep their deflines and no
  letters, and the search reads the stores as it reads a database's sequences (`Subjects::Stored`): the
  preliminary frames from the packed bases (NCBI's compressed subject: an ambiguous base resolved by
  `CRandom(length)`, the bases `x_resolved_record_frames` gives a record), the letter frames and the report rows
  from the stored letters (`report::SubjectLetters::Stored`, ordinals `0..n`); the numbers of the subject set, the
  warnings of records without letters and the outfmt 0 prolog read the stores' lengths. Before, the letters of
  the whole set (one byte per base) were held for the whole run. Conditional: `run_local` (LOSAT Web, the
  library) keeps its records, and a set with a byte without an ncbi4na code keeps every record's letters (the
  CLI reader gives none, it reads `U` as `T`: the reference side is a gap of the reach table). The same pull
  request made the packing and the stored windows faster for every caller, with no switch (the switch alone
  cost time, below, and a fix that costs time lands only with its recovery): `PackedNucleotide::
  from_subject_letters` packs blocks of 32 letters of one case at a time (the check vectorized, four letters by
  one multiplication), and `bases_into` / `letters_to` make a window's letters with one table read per packed
  byte. The lib tests `every_form_is_the_one_made_from_the_letters` (the packed bytes against
  `encode_subject_ncbi2na_packed`, `bases_into` against `codes_into`, letters and intervals, with block-shaped
  cases: one case per block, case changes and other letters inside and between blocks, N runs across blocks,
  every length around a block) pin them. Measured on the workstation (node 0 CPUs 16-23, HT siblings idle, cores warmed 3 s, one warm-up round +
  3 rounds, arm order rotating, medians of the per-round ratios against `=0` of the same binary; reports
  identical; load 21-35): 3 fruit fly transcripts x the fruit fly genome FASTA (`-subject`, 1,870 records) 1
  thread max RSS 218 -> 114 MiB (the genome as `-db`: 107 MiB), cycles x0.990, wall x0.984; 8 threads 409 ->
  316 MiB, cycles x0.992, wall x0.992; x the fruit fly cDNA FASTA (31,168 records) 1 thread 105 -> 43 MiB,
  cycles x0.985; 8 threads 106 -> 46 MiB, cycles x0.992; the genome as `-db` against the dev binary cycles
  x0.985 (the stored windows). Without the faster packing the switch cost cycles x1.016 (1 thread) and x1.019
  (8 threads, wall up to x1.13: the packing thread was slower than the reader) for the same memory. With
  `LOSAT_TIMING` the packing is in `read_fasta`. WebAssembly (serial WASI command, one Node process, a fresh
  instance per search, CPU 20, one warm-up round + 3 rounds against `=0`; output identical): the fruit fly
  genome FASTA as `-subject` 11.01 -> 10.82 s (x0.981), linear memory 206 -> 128 MiB.

Scheduling (each value computed by one thread with the original function; order of output kept):

- AHEAD, GREEDYSPEC (graduated 2026-10-22; no variables, no `LOSAT_X_AHEADSHADOW`; class
  conditional, the names stay in `GRADUATED` for the reach table) blastn/megablast ordered
  gapped-extension and traceback loops: helpers compute the next elements' alignments ahead, the
  tree is read and written in the original order by one thread (`utils/xahead.rs`; the window is
  the constant `WINDOW`, 16). The rule stays: a parallel search (2+ threads) without
  `LOSAT_DEBUG_COORDS` whose aligner is the dynamic-programming one (`config.use_dp`) or, for
  GREEDYSPEC, megablast's greedy one; and the ordered loop only where `LOSAT_X_QGROUPPAR` leaves
  it. GREEDYSPEC is the megablast greedy traceback in that schedule; `LOSAT_X_SPECBATCH` tunes the batch.
- PAIRPAR, blastn/megablast half (graduated 2026-10-22; no variable, no `LOSAT_X_PAIRPARSHADOW`;
  class conditional, the name stays in `GRADUATED` for the reach table; it was in the strict set; the
  TBLASTX half is below) one query-subject pair on a pool of 2+ threads. (1) The megablast lookup
  table is filled by the pool threads: every thread reads the whole query in order and applies the
  update only to the words whose cell lies in its blocks of `2^max(11, pv_array_bts)` cells, so each
  cell gets its words in the original order (`blastn/x_mb_lookup_par.rs`; not with database word
  counts or `BLEMIR_DEBUG`; the serial `in_pool == false` branch of the fill stays for the unit
  tests). (2) The one-hit seed stage of a subject chunk with the diagonal hash (query block over 8000
  letters, `window_size` 0, `scan_range` 0, contiguous megablast table, unmasked subject, one
  subject at a time on the owner (several one at a time under `LOSAT_X_BNFEWSUBJ`), no small-query
  word extension, no debug or trace): the subject scan grid is cut into strips that the pool threads
  scan (lookup chains and word extension), and the owner thread folds the strips in scan order
  through a copy of the diagonal-hash code (`blastn/blast_engine/x_pair_seed.rs`). Anything else
  runs the original code. The shadow (lookup tables and,
  per chunk, the hit list and the whole diagonal hash compared with the original code, with its
  `[X_PAIRPARSHADOW]` summary line) went with the variable, with the second build of the megablast
  table that was the only run-time comparison of MBLUTONCE and MBLUTPREFAULT, the one-hit
  comparison of PAIRPARFOLD and the `max_level` check of DIAGLV for this path; their `=0`
  comparisons in the tests remain. With
  `LOSAT_TIMING` a `[TIMING] x_pairpar_seed` line gives the strips, seeds and the owner's fold, scan
  and wait time.
- Private copies of the megablast fill's counters (`blastn/x_mb_lookup_par.rs`; used by the read-once fill of
  `LOSAT_X_MBLUTONCE` below; the fill of PAIRPAR itself updates the shared arrays). A switch that gave the
  PAIRPAR fill the same copies with `LOSAT_X_MBLUTONCE=0` was removed because it never acted in a default
  run; its code and tests went with it, and the output of every path that remains is unchanged. The counters of
  the longest chain (`helper_array`, one per 2048 cells) and the presence vector (`pv_array`) are counted and set by
  each thread in private copies, a whole number of 64-byte lines per thread in one buffer of the calling thread (each
  thread zeroes its own lines first), and each thread adds its own blocks' counters (`saturating_add`) and ORs its own
  blocks' presence-vector words into the tables (`flush_private`). The blocks of PAIRPAR's plan are dealt
  round robin, so a line of `helper_array` held the counters of up to 16 threads and a line of `pv_array` the words of
  up to 8 (`pv_array_bts` 10): false sharing. A counter or presence-vector word lies inside one block, so only its
  owner ever changes it, from 0, with the loop's increments and bits; those values, `longest_chain`, `hashtable` and
  `next_pos` are the loop's. Memory: per thread `ceil(hashsize / 32768)` + `ceil(pv_size / 16)` lines of 64 bytes (160
  or 96 KiB at `lut_word_length` 12, 136 or 72 KiB at 11, by `pv_array_bts`), all written, so a fill of P threads
  holds P times that: at `lut_word_length` 12 2.5 / 1.5 MiB at 16 threads, 20 / 12 MiB at 128 and 40 / 24 MiB at 255
  (the parts are at most 255), at 11 2.1 / 1.1, 17 / 9 and 34 / 18 MiB. On Linux with glibc the buffer is kept for the
  next fill and never freed: freeing it after each fill made glibc raise its mmap threshold to its size, after which
  the search's other allocations of up to that size stayed in malloc's arenas (megablast of 31,168 cDNAs at 64
  threads: max RSS 263-305 MiB against 221-225). Each buffer is made with a capacity of at least 32 MiB + 4 KiB, above
  glibc's largest dynamic threshold (the rule of `LOSAT_X_TBNSTORE`'s frames), so glibc maps it, freeing one never
  moves the threshold, and only the lines the fills wrote are resident; a kept buffer too small for a fill (only above
  32 MiB of copies, 205 threads or more at 160 KiB) is freed and a new one made, never grown (std reallocates a
  64-byte-aligned block by allocating, copying and freeing, which moves the threshold). The lock of the kept buffers
  is only tried (`try_lock`), never waited for: a fill that finds it taken uses a new buffer, which it keeps only if
  the lock is free at its end. On other targets (macOS, Windows, musl, wasm32) the buffer has the fill's size and is
  freed at the end of the fill (their allocators return large blocks themselves). With `LOSAT_TIMING` `main` prints
  `[TIMING] x_mb_lookup_build: <s>s (builds=<n> parallel=<n> once=<n>)` when a megablast table was built (the summed
  wall time of the builds, not the shadow's second build;
  `parallel=` counts the builds filled on a pool) and `[TIMING] x_mb_lookup_split: alloc= loop= pool_fill= chain=`
  (the allocation, the fill by the loop or on the pool with the re-threading of `ascending_cells`, and the longest
  chain; seconds, summed).
- `LOSAT_X_MBLUTONCE=1` (BLASTN megablast and blastn wherever PAIRPAR fills the contiguous megablast lookup
  table on a pool of two threads or more: `-subject` and `-db`, not with `-limit_lookup` or `BLEMIR_DEBUG`; in the
  strict set since 2026-10-18 (perf/bn-memory-rules, evidence at the end of this entry, "Default ON"): on unless it
  is `0`) each
  query position of that fill is read once (`blastn/x_mb_lookup_once.rs`, the first
  dispatch of `x_mb_lookup_par::fill`): the positions the loop reads form one stream cut into chunks of
  `min(2^19, P << 15)` positions (P the fill's threads), each chunk cut into one span per thread; per chunk each thread reads the words of its span with the
  loop's statements (entering a range up to `lut_word_length - 1` positions early) into a region of its own, copies
  them stably into a second region grouped by the owner of their cell (the blocks of PAIRPAR's plan), each
  owner then applies the words of its cells in reading order into the private copies (above) and keeps
  each cell's previous head, and each reader writes `next_pos` of its words from the kept heads (in the same job as
  its read of the next chunk; two `rayon::scope`s per chunk, a last one that links the last chunk and adds the copies
  in). Each cell gets the loop's words in the loop's order, so `hashtable`, `next_pos`, `pv_array`, `helper_array` and
  `longest_chain` are the loop's; for blastn (`ascending_cells`) the chunks, readers and words are taken in reverse,
  the order of the loop's re-threading, so its per-thread word lists are no longer needed. Instead of every thread
  reading every word (0.65 s of each thread's time at 16 threads on megablast of 31,168 cDNAs), each word is read once
  and written to `next_pos` by the thread that read it. Memory: one arena per fill, P regions of `2 x ceil(ceil(c / P)
  / 8)` lines of 64 bytes and P rows of `ceil((P + 1) / 16)` lines, `c` the smaller of the chunk and the positions
  read: at most `16 x (c + 8 P)` bytes plus the rows, about P x 0.5 MiB up to 16 threads (2.0 MiB at 4, 8.0 MiB at
  16) and 8.0 MiB at 16 to 64 threads (8.3 MiB at 255) for a query of 2^19 positions or more, less for smaller
  queries. That is per fill, not per process: fills that run at the same time each hold their own arena. Until
  2026-10-18 the chunk was 2^19 positions at every thread count and k concurrent fills (the slots of
  `LOSAT_X_QBATCHPAR`, removed on 2026-10-28) held about k x 8 MiB (three E. coli genomes against LC738868 at 16
  threads, 3 slots: +25 MiB against the loop's fill). On Linux with glibc the arena's capacity is at least 32 MiB + 4 KiB (the rule of `LOSAT_X_TBNSTORE`'s
  frames), so glibc maps it and unmaps it at the end of the fill without moving its mmap threshold and only the
  written part is resident; elsewhere (wasm32 included) it has its exact size and is freed at the end of the fill. The
  ranges of the stream are a list of one 48-byte entry per unmasked range read (0.1-0.3 MB for the megablast tables of
  31,168 cDNAs), freed with the fill; the private copies are those above, one set per thread of each fill
  (kept from fill to fill on glibc). With the `wasm-threads` feature the fill runs on the search pool (it builds no
  pool); without it, or without `parallel`, the loop runs. Until 2026-10-22 the shadow of PAIRPAR
  (`LOSAT_X_PAIRPARSHADOW=1`) made the second build of its comparison with the read-once fill and `main` printed `[X_SHADOW]
  MBLUTONCE tables= cells= diffs=0 identical`; it went with that variable, and the read-once fill is now compared with the
  loop's by the `=0` runs of `LOSAT/tests/blastn_mb_lookup_once.rs` and by its unit tests. With `LOSAT_TIMING` `[TIMING] x_mb_lookup_build` counts `once=<n>` (the builds the
  read-once fill made), and `[TIMING] x_mb_lookup_once:
  fills= chunks= read_wall= read_max= read_sum= read_late= update_... last_... prefault_... arena= prefault_rounds=
  prefault_touches=` gives, per phase summed over the chunks of the search's fills, the calling
  worker's wall time, the longest job, the jobs' summed time and the latest start (`read`: the read and the links of
  the chunk before; `update`: the owners' updates; `last`: the last links and the copies' flush; `prefault`: the
  prefault round of `LOSAT_X_MBLUTPREFAULT`), the largest arena (MiB), and the fills that ran that round and the
  zeros it wrote.
  Measured 2026-10-18 on the workstation (da08ced7, the same release code; 16 physical cores of node 1 with their
  siblings idle for 2-16 threads, 32 with their siblings for 64; a warm-up, then 3-5 rounds of wall runs and 1-2 of
  `LOSAT_TIMING` runs per cell, arms rotated: base without the switch, the switch, a since-removed switch alone
  (private copies in the PAIRPAR fill), and the switch with chunks of 2^18; every output had its expected sha256; perf-mb `logs/mbonce_ws_table_mbonce2_a.md` and
  `_b.md`), switch / base: megablast of 31,168 cDNAs (37 tables) fill 0.64 / 0.54 / 0.52 / 0.39 per round at 2 / 4 /
  16 / 64 threads (2.99 -> 1.92, 2.23 -> 1.21, 1.50 -> 0.79 and 1.71 -> 0.66 s), search 0.838 / 0.816 / 0.849 / 0.826;
  with the slots of `LOSAT_X_QBATCHPAR` (removed on 2026-10-28; the dsw candidate set) 0.843 at 16 threads, and 0.837 / 0.888 with 4 / 8
  slots at 64; megablast of 1,775 cDNAs fill 0.49, search 0.882 (16 threads); blastn of 31,168 cDNAs fill 0.43, search
  0.967; blastn of 1,775 cDNAs fill 0.59, search 1.064 (0.77-1.30 per round, the spread of its 1.9 s searches;
  the removed switch alone 0.988, 0.79-1.19); the three E. coli genomes against LC738868 fill 0.55, search 1.021
  (0.84-1.04). Max RSS medians base -> switch (MiB; NCBI BLAST+ 2.17.0 at one thread in brackets): 31,168 cDNAs 250.1
  -> 259.0, 249.2 -> 255.6, 244.5 -> 254.6 and 222.4 -> 253.2 at 2 / 4 / 16 / 64 threads [358.5], with 4 / 8 slots at
  64 threads 580 -> 622 and 1,026 -> 1,050 (wide ranges); 1,775 cDNAs by megablast 126.5 -> 126.5 [198.7] and by
  blastn 141.1 -> 145.3 [143.5]; blastn of 31,168 cDNAs 571.5 -> 576.0; E. coli 143.4 -> 153.3 [185.1]. The fill's
  system time rose (31,168 cDNAs at 64 threads 22.8 -> 29.0 s summed). Chunks of 2^18 (an arena of 4 MiB): fill 0.60
  against 0.52 at 16 threads and 0.36 against 0.39 at 64, max RSS 1-8 MiB below the switch's. So the switch wins time
  where megablast fills large tables (search 0.82-0.89 at 2-64 threads), ties on blastn and on the E. coli input, and
  costs max RSS: 0 to +10 MiB at 2-16 threads and +31 MiB at 64, above NCBI's for blastn of 1,775 cDNAs at 16 threads
  (145.3 against 143.5). It stayed off by default then (no rule separating the RSS cost was tried).
  Default ON since 2026-10-18 (perf/bn-memory-rules; workstation, physical cores of one node with their siblings
  idle (32 cores and their siblings for 64 threads), a warm-up and 4-8 rounds per cell with the arms rotated, max RSS
  of both arms in the same window, cycles:u by `perf stat`; v4-D `results/d10_w1_summary.md`, `d10_w2_summary.md`,
  `d10_w3c_summary.md`; switch / without, wall per round, and max RSS difference): with the 2^19 chunk megablast of
  1,775 cDNAs 0.783 at 16 threads and 0.764 at 4, of 31,168 cDNAs 0.818 (16), E. coli 0.823 (16) but +25.2 MiB
  (three slots of 8 MiB arenas), blastn of 1,775 cDNAs 0.973 (16) and 0.958 (4), of 31,168 cDNAs 0.996 (16); with the
  chunk by workers (`min(2^19, P << 15)`): E. coli 0.836 at +4.1 MiB (16), 0.827 at +17.8 MiB (32) and 0.823 at +24.4
  MiB (64), megablast of 31,168 cDNAs 0.870 at +4.1 MiB (16) and 0.876 at +22.7 MiB (64), of 1,775 cDNAs 0.835 at
  +1.2 MiB (4), blastn of 1,775 cDNAs 0.954 at +1.1 MiB (4). No measured cell loses time (64 threads included, where
  the removed private-copies switch lost on few-table inputs), and every cell of these windows with an NCBI figure stays below
  NCBI's max RSS (blastn of 1,775 cDNAs at 16 threads 113.0 MiB with the switch, NCBI 143.5 at one thread); the
  chunk size changes no table element (the unit tests fill with chunks of 1 to 2^19 positions).
- `LOSAT_X_MBLUTPREFAULT=1` (BLASTN wherever the read-once fill of `LOSAT_X_MBLUTONCE` builds the contiguous
  megablast lookup table; in the strict set since 2026-10-19 (perf/mb-fill-prefault, evidence at the end of this
  entry): on unless it is `0`) a prefault
  round before the fill's first update (`blastn/x_mb_lookup_once.rs`, the form of `LOSAT_X_DISCLUTPAR`'s round below,
  page by page): when the first chunk of the processing order has at least 8 positions per 1,024 cells of the table
  (131,072 for 4^12 cells, the first chunk of 4 threads), its readers mark the 4 KiB page of each word's cell, and the
  fill's threads write 0 to the middle element of each marked page (one thread per page, the owner of that element's
  block) before any cell is updated. The update reads a cell before it writes it, so without the round a fresh page is
  first mapped to the kernel's zero page and then copied on write, which flushes its entry from the TLBs of every other
  CPU running the search; written first, it takes one fault and no flush. Each element written holds the 0 it was
  allocated with (the table is zeroed, NCBI's `calloc`), and only pages the fill writes are marked, so every array
  and every report is what it is without the round, and at 4 KiB pages it maps no page the fill does not map (no
  `next_pos` page: the links are written without a read). The saving rests on Linux's zero page and copy-on-write;
  the timing evidence is Linux only. Memory: the readers' marks, one bit per page per reader on whole cache lines
  (2,112 bytes per reader for a 4^12 table, 33 KiB at 16 threads), and their OR, freed with the fill. With
  `LOSAT_TIMING` `[TIMING] x_mb_lookup_once` adds `prefault_wall= prefault_max= prefault_sum= prefault_late=` (the
  round, as the other phases) and `prefault_rounds= prefault_touches=` (fills that ran it, zeros written). (Until
  2026-10-22 the shadow of `LOSAT_X_PAIRPARSHADOW` built its read-once table with the round and compared every element.)
  Evidence. Line R's R5 (workstation, Linux 6.8, THP `madvise`, `v4-R/REPORT_R5.md`): the fills of three E. coli
  genomes (65 records) against LC738868 at 64 threads made 102,000 such copies-on-write, 6.4 million TLB shootdowns
  to the search's CPUs; a probe that wrote every page before the fill removed them (2) and took whole runs to
  0.73-0.94 of the read-once fill's time at 16-64 threads. This round, on the i9 (2026-10-19, ws unreachable; WSL2,
  8 alternating pairs per cell in one run, load 4-6, every report of a cell identical, the cDNAs' equal to NCBI
  BLAST+ 2.17.0's; v4-D `results/i9_ab*_summary.md`), round / without, wall per round, and max RSS (MiB): E. coli x
  LC738868 0.69 at 4 threads (141 -> 141), 0.64 at 8 (228 -> 241), 0.61 at 16 (330 -> 349), 0.63 at 32 (332 -> 354),
  0.73 at 16 with `LOSAT_X_QBATCHPAR=0` (the batch loop, now the only form; 148 -> 148); MG1655 x LC738868 (one table) 0.88 at 16 (129 -> 129) and
  0.93 at 32 (130 -> 130); megablast of 31,168 D. melanogaster cDNAs against 3 Mb of chromosome 2L 0.79 at 8 (382 ->
  383), 0.78 at 16 (620 -> 604), 0.78 at 32 (856 -> 830), 0.86 at 16 and 0.92 at 32 with `LOSAT_X_QBATCHPAR=0` (258 -> 257, 262 -> 264);
  1.00 at one thread (no pool fill). The search's system time fell in every cell (0.25 -> 0.08 s on MG1655 at 16 threads, 6.29 -> 2.43 on the cDNAs at 32). Max RSS moved only where
  `LOSAT_X_QBATCHPAR`'s slots (removed on 2026-10-28) filled tables side by side (+13 to +21 MiB on E. coli at 8-32 threads, -26 to +1 on the
  cDNAs, spreads of 20-40 MiB): the round brings a slot's table to full residence at the start of its fill instead
  of during it, within the batch memory the slots already admit. Those cells are above NCBI's max RSS without the
  round too (NCBI with `-subject` does not split by queries and stays at 183 MiB for E. coli and 356 for the cDNAs at
  1-16 threads), because of the slots, not of the round. On this WSL machine (THP `always`) megablast of 5 kb to 4.6
  Mb of MG1655 against LC738868 at 8 threads kept its max RSS within 0.4 MiB of the fill without the round, at 0.57-0.86
  of the time from 150 kb on.
- `LOSAT_X_MB24=1` (BLASTN megablast, the contiguous megablast table of `build_mb_lookup` with 4^12 cells; memory; in the
  strict set from the pull request that added it, 2026-10-26 (perf/mb-narrow); `=0` builds the 4-byte table) keeps the
  table's head table (`hashtable`) and chains (`next_pos`) with elements of 3 bytes instead of 4 (`blastn/x_mb_narrow.rs`,
  argument in its header). NCBI's `Int4` arrays (`blast_nalookup.c:1262`, `:979`) hold 1-based query offsets, all below
  2^24 when `next_pos` has at most 2^24 elements, so element `i` is the low 24 bits of the little-endian `u32` at byte
  `3 i` of a `Vec<u32>` (`(3 n + 4) / 4` words for `n` elements, so the 4-byte read of the last one stays inside): the
  heads of a 4^12 table take 48 MiB instead of 64 and the chains 3 bytes a query position instead of 4 (Sakai's both
  strands, 11.0 M positions: 31.5 MiB instead of 42). Where: a table of 4^12 cells (megablast of a query batch of about
  300,000 table entries or more), newest first (not the blastn task's re-threaded order), without database word counts or
  `BLEMIR_DEBUG`, `next_pos` of at most 2^24 elements (a batch of up to about 16.7 M positions of both strands), queries
  whose offsets ascend without overlap and no inverted mask; the smaller tables, the blastn task and larger batches keep
  the 4-byte table. Every fill applies the same updates in the same order in an instance of its own (`H24`): the serial
  fills (`LOSAT_X_MBFILL`'s prefetch and batch queues, `x_fill_mb_queued`) read a head as 4 bytes and write it back with
  the next cell's byte, and write a link as one 4-byte store whose fourth byte, the next element's first, is not written
  yet (the words come in increasing order of offset); the pool fills (PAIRPAR's and `LOSAT_X_MBLUTONCE`'s) read and write
  every element as exactly its 3 bytes, since the next cell or offset may be another worker's; the prefault round of
  `LOSAT_X_MBLUTPREFAULT` marks the pages of a cell's first and last byte (12,289 pages of a 4^12 table, at least 98,304
  positions in the first chunk). The heads are the batch loop's kept table (`LOSAT_X_MBLUTKEEP`) read as bytes: a build
  takes its first 12,582,913 words; after a 4-byte 4^12 build in the same loop the keep holds the 64 MiB table, and the
  3-byte builds use its start. The readers (`for_each_hit`, `contains_hit`, `count_hits`) find a 3-byte table behind the
  bounds check of the head that already sends the CSR tables of `LOSAT_X_BNCSR` elsewhere, so the scans of 4-byte tables
  run the same instructions (`=0` against `dev`: x1.0001-1.0010 instructions; the blastn task, 4-byte tables only,
  x1.0000). Evidence (workstation; one node 0 core with its HT sibling idle, or node 1
  CPUs 32-47 with theirs, membind to the node, cores warmed 3 s before each run; arms rotated every round; `dev`
  ce281c9a, outputs identical in every arm; peak `RssAnon` and max RSS = max(`%M`, a 10 ms `VmHWM` sample)): megablast of
  E. coli Sakai against MG1655 at one thread `RssAnon` 126.5 -> 99.9 MiB (max RSS 131.6 -> 105.4), wall x1.004 of `=0`
  in 16 rounds (x0.954 in 8 rounds under `LOSAT_TIMING`: the build x0.92, the search x0.96), `cycles:u` x1.01-1.035
  (instructions x1.021: the 3-byte stores and reads, about 9 a word), system time 0.08 -> 0.06 s; at 16 threads
  `-subject` 129.4 -> 102.4 MiB (132.7 -> 102.2), wall x0.995, `-db` 129.6 -> 102.3 (131.3 -> 99.7), wall x0.987; megablast of 31,168 D. melanogaster cDNAs against dmel_dna_3M (37 batches, 18 of
  them 4^12 tables) at one thread 156.2 -> 130.6 MiB (163.0 -> 137.1), wall x0.993, `cycles:u` x1.037, system time 0.45 ->
  0.13 s (minor faults 232 k -> 64 k), and with `-db` at 16 threads 174.7 -> 169.6 MiB (177.2 -> 173.1), wall x0.993;
  blastn of the cDNAs against MG1655 (no 3-byte table) wall x0.988. Wasm (the serial WASI command in Node 26.8.2, all arms'
  modules in one process, a fresh instance per search, 8 rounds; two node 0 cores): the linear memory of Sakai against
  MG1655 133.1 -> 106.4 MiB, x0.995 of the time, of the cDNAs 234.3 -> 208.8 MiB, x0.992.
- `LOSAT_X_MBCSR=1` (BLASTN megablast, the contiguous tables of 4^12 cells with a packed head, `-subject` and `-db`, any
  thread count, not on wasm32; speed and memory; in the strict set from the pull request that added it, 2026-10-30
  (perf/mb-csr12): on unless it is `0`) lays those tables out as CSR with a packed head (`blastn/x_bn_csr.rs`, `x_in_scope`), as
  `LOSAT_X_BNCSR` and `LOSAT_X_BNCSRPACK` do for 4^8 to 4^11 cells, when the subjects' letters are at least 8 times the
  positions of the word loop plus the cells (`X_MB12_LETTERS_PER_POSITION`). The scan then reads a word's presence bits and its run
  through the prefetch ring instead of one miss per head and per chain element on `next_pos`, which took 45 % of the
  cycles of draft200 x one volume of ref_prok_rep_genomes (4 Gnt) at one thread. The rule keeps the chains' cheaper
  fill for small subjects: without it Sakai x MG1655 was ×1.34 and draft200 x an 8.7 MB database ×1.69 (the CSR's build
  costs 13-22 ns more per position and about 40 ms per build, its scan saves about 3.2 ns per subject letter). v4-R R11 measured only a 5 Mb
  subject and found no gain. Evidence (ws, one core of node 0, three rounds, the same reports): draft200 x that volume
  ×0.72 of `=0` with max RSS ×0.90, 235 contigs x 1.5 Gnt of plasmids ×0.74 with ×0.86, Meles x 129 Mb of nt ×0.81 with
  ×0.92; below the rule the same as `=0`. `LOSAT/tests/blastn_mb_csr.rs` checks both sides of the rule on whole
  searches.
- `LOSAT_X_MBHITPF=1` (BLASTN, every task, x86_64, the loop over the offset pairs of `search_query_batch`, not the seed
  stage of `LOSAT_X_PAIRPAR`; speed; in the strict set from the pull request that added it, 2026-10-30 (perf/mb-csr12): on unless it is `0`) prefetches the query line of
  the hit 16 places on (`x_mb_hit_prefetch`, `X_HIT_AHEAD`). A hit's first extension reads the query at a place its
  lookup cell chose, which took about 20 % of the cycles of draft200 x one volume of ref_prok_rep_genomes at one thread.
  Evidence (ws, one core of node 0, four rounds, the same reports, max RSS equal): that case ×0.922, 235 contigs x
  1.5 Gnt of plasmids ×0.942, Meles x 129 Mb of nt ×0.945, Sakai x MG1655 ×0.982. Prefetching the hit's diagonal cell as
  well gave ×1.000 and was dropped.
- `LOSAT_X_MBFILL=1` (BLASTN megablast and blastn, the one-thread fill of `build_mb_lookup`: the contiguous megablast
  table and the two-stage blastn table; speed; in the strict set from the pull request that added it, 2026-10-23; `=0`
  runs the queued loop `x_fill_mb_queued`) reads the words of the fill run by run
  (`blastn/x_mb_fill.rs`, argument in its header): the first `lut_word_length - 1` bases of each run of unambiguous bases
  build the code, every later base of the run ends a word, and the words are NCBI's words in NCBI's order; each update
  (`blast_nalookup.c:1081-1091`, the same statements) is applied 64 words after a prefetch of its cell (the queued loop's
  `Delay` queue was 16 deep). The loop over a run is a function of its own, so the tables and the queue stay in
  registers; the indexes are checked once per build instead of per word, and the presence bit comes from a table, so the
  shift by `pv_array_bts` keeps `cl` (x86-64 without BMI2). Where: x86_64 and wasm32, a
  table built off a pool of two threads or more (PAIRPAR's pool fill takes the others), without database word counts
  (`limit_lookup`, which the command line does not offer). wasm32 has no prefetch instruction, so there the queue is a
  batch (since 2026-10-24, perf/mb-fill-wasm): 64 words are collected, their cells read in one loop (the cache misses
  overlap), then the words applied in order, as the queued loop's `Batch` queue does; other targets (aarch64) keep the
  queued loop, where the run loop has not been timed. R10 (line R, megablast Sakai x MG1655 at one thread): the queued loop spent 0.190 s, about 87
  instructions and 42 cycles a word, on 10.9 M words (NCBI 45.5 cycles); the loop shape of R10's probe (a closure and a
  queue 32 deep) gained 11 ms, and a branch-free update of the presence bit and helper counter lost here as in R10
  (megablast of the cDNAs below x1.004: about 21 more instructions a word for few saved misses).
  Measured 2026-10-23 on the workstation (node 1, one pinned core, HT sibling idle, warmed 3 s; 5 rounds, the arm order
  reversed every other round; `perf stat` cycles against the binary before, `dev` ceaef53e; outputs and max RSS
  identical in every arm): megablast of Sakai against MG1655 x0.906 (0.683 -> 0.623 s; the fill 0.189 -> 0.131 s under
  `LOSAT_TIMING`, 9 % fewer instructions), of 31,168 D. melanogaster cDNAs against dmel_dna_3M x0.874 (5.73 -> 5.06 s,
  37 table builds), of Sakai with every 20th letter N (runs of 19 letters, 8 words a run) x0.907, blastn of the cDNAs
  against MG1655 x0.985 (91 two-stage builds); megablast of Sakai at 8 threads (pool fill) x1.001. `=0` x0.976-0.999
  of the binary before (0.2-2.4 % fewer instructions: the queued loop sits behind the dispatch now).
  Measured on wasm32 2026-10-24 (the serial WASI command in Node 26.8.2, every arm's module in one process, a fresh
  instance per search, 8 rounds, the arm order rotating; `LOSAT_TIMING`; against `dev` b1f09848; outputs identical):
  on the workstation (node 0, two pinned cores, HT siblings idle) megablast of Sakai against MG1655 x0.957 (1.403 ->
  1.340 s; the fill 0.237 -> 0.182 s), of the 31,168 cDNAs against dmel_dna_3M x0.891 (8.37 -> 7.46 s; the fill 3.15
  -> 2.22 s), of Sakai with every 20th letter N x0.924, of 21,994 150-letter pieces of Sakai against MG1655 x0.975;
  `=0` x0.987-0.994. On this WSL i9 (a client CPU, loaded) the fill of the same loop read x0.74-0.84 of the queued
  loop's in every case. Applying each word at once, or batching only the first 64 words of a run, was faster on the
  workstation for long runs (the cDNAs x0.82) but slower than the queued loop on the i9 (short runs x1.24, Sakai x1.07),
  so wasm32 takes the batch. The machine code of the x86_64 loops is unchanged.
- `LOSAT_X_MBLUTKEEP=1` (BLASTN megablast, blastn and dc-megablast: the array of one element per cell of a batch's
  lookup table, the contiguous megablast table's `hashtable` of `build_mb_lookup` and since perf/bn-keep-cells a CSR
  table's `offsets` and a one-template discontiguous table's chain heads; in the strict set since 2026-10-19
  (perf/mb-lut-keep, evidence at the end of this entry): on unless it is `0`)
  the batch loop keeps the table's `hashtable` from one query batch to the next
  (`blastn/x_mb_lut_keep.rs`). NCBI allocates each batch's table zeroed (`calloc`, `blast_nalookup.c:1262`) and frees it
  with the batch (`sfree`, `:1376`); a 4^12-cell table is 64 MiB, which glibc maps per batch and unmaps at the batch's
  end, so each batch faulted the table's pages in again (or in the prefault round of `LOSAT_X_MBLUTPREFAULT`) and each
  free unmapped them and flushed their TLB entries on every CPU of the search. With the switch the loop holds the
  longest table built so far; the next build takes it when it is at least as long as the build's table, writes 0 to the
  cells of its table that earlier builds wrote (on the search pool when the build runs on a pool of two threads or more;
  the other cells are 0 already) and uses its first cells as its table. NCBI's `CBatchSizeMixer` alternates full batches
  with small probe batches (megablast of 31,168 cDNAs: 18 tables of 4^12 cells, 19 of 4^8 or 4^11), so a probe batch's
  table is the start of the kept one, and no build allocates a second table beside it (a longer build frees the kept one
  first and allocates). The keep records the longest start of the table that builds wrote since it was last zeroed, so a
  probe batch zeroes its own cells only and the rest of a full batch's cells wait for the next build that uses them
  (megablast of 1,775 cDNAs builds tables of 4^11, 4^8 and 4^12 cells; zeroing the whole 4^11 start for the 4^8 build
  cost 1.7 % of the run at one thread). The prefault round does not run on a kept table (its pages are mapped; the round
  writes zeros over zeros): the zeroing writes every cell of the build, also those past the start that builds wrote,
  which a packed CSR of `LOSAT_X_BNCSRPACK` leaves unwritten (a fill that read such a page first mapped the zero page
  and copied it on its write, one TLB flush of every CPU of the search per page). The fill finds every cell of its
  table at 0, as in a fresh one, and every reader takes cells below the table's size, so every array and every report
  is what it is without the switch; `next_pos`, `pv_array` and
  the counters are allocated as before. Where: only the batch loop of `run_in_pool`, for the batches it searches itself
  and the query chunks of a split batch. Not the discontiguous tables of two templates, the PAIRPAR sub-lookups or the tests.
  Memory: between batches the loop holds the longest table instead of freeing it; a probe batch's peak with it stays
  below the full batch's. Since perf/bn-keep-cells the keep also serves the tables that have no `hashtable`: a build
  laid out as CSR (`LOSAT_X_BNCSR`, `LOSAT_X_DCCSR`) builds its `offsets` (one element per cell and one or two more) in
  the kept table (`Keep::take_buffer`; the build writes 0 to every element it uses before counting into it) and gives
  it back when the batch drops the CSR, where it freed the kept table before; a one-template discontiguous table takes
  its chain heads from the keep as the megablast table does. These arrays depend on the table's cells only (4^11 or
  4^12), so one kept table serves the CSR builds of `CBatchSizeMixer`'s full batches and the chain builds of its probe
  batches (dc-megablast of 31,168 cDNAs: 18 CSR and 19 chain builds of 4^11 cells); a fresh table has two elements more
  than its cells for a later CSR. Why: glibc raises its mmap threshold to the size of a mapped block it frees (up to
  32 MiB); a 16 MiB array freed after each batch put it up after the first batch, the heap then took the batches' other
  buffers and the free space in its arena stayed resident (the dc-megablast above at one thread: anonymous RSS 84-92
  MiB depending on the build, against 80.8 for the memory floor of plan v4 line MF). With `LOSAT_TIMING` `[TIMING]
  x_mb_lut_keep: kept= fresh= csr= zero=` gives the loop's builds of a `hashtable` or chain heads on the kept table, its
  builds (of either layout) on a fresh one, its CSR builds on the kept table and the time spent zeroing (also in
  `x_mb_lookup_split`'s `alloc=`). Evidence. Premise probe (v4-D
  `results/d10d2_probe1_summary.md`, workstation node 1, `LOSAT_X_QBATCHPAR=0`): megablast of 31,168 D. melanogaster
  cDNAs against 3 Mb of chromosome 2L at 16 threads spent 0.141 s unmapping the batches' tables at their ends (serial,
  outside `[TIMING] total`) and 0.060 s in prefault rounds; keeping the table took 0.038 s of zeroing instead, with
  instructions:u -0.4 %, system time 2.24 -> 1.30 s, page faults 558K -> 268K and TLB shootdowns 152K -> 53K. Speed
  (workstation node 1, windows v4d9-w2 and v4d9-w3 against dev 80405895 and 9bc7c3cb, `results/d10d2_w2_summary.md` and
  `d10d2_w3_summary.md`, medians of rounds in alternating order on warmed CPUs): in the batch loop
  (`LOSAT_X_QBATCHPAR=0`, against dev's loop) megablast of 100-250 kb E. coli pieces against LC738868 took 0.84-0.85 of
  the time at 8-32 threads and 0.80-0.82 at one thread, megablast of the 31,168 cDNAs 0.85-0.88 at 8-32 threads and
  0.90-0.91 at one, blastn and dc-megablast of cDNAs 0.97-1.02; by default (`LOSAT_X_QBATCHPAR`'s slots, removed on
  2026-10-28, which the switch left alone, at 8-32 threads) 0.97-1.02; megablast of 1,775 cDNAs (the loop by default) 1.00 at one thread and
  1.02-1.03 at 8 and 16 threads, within the rounds' spread (0.97-1.09) with the same page faults, TLB shootdowns, user
  and system time as dev. Max RSS: the same (the 31,168 cDNAs in the loop 236 MiB against 249). Evidence of
  perf/bn-keep-cells (CSR `offsets` and discontiguous chain heads; workstation, 2026-10-27, against dev c780af86, plan
  v4 line D sessions 20-22; the cDNAs against the first 3 Mb of chromosome 2L `-subject`, one sha per cell): one thread
  on a node-0 core, 6 rounds in rotated order, median of the per-round ratios: dc-megablast of the 31,168 cDNAs peak
  `RssAnon` 91.0 -> 79.4 MiB (the floor 80.8), time x0.997 (0.958-1.026), minor faults 192 K -> 183 K; blastn of them
  62.1 -> 60.2 MiB, x1.003 (0.977-1.028), system time 0.41 -> 0.55 s and minor faults 146 K -> 231 K: glibc's threshold
  now stops at 7.7 MB (a batch's CSR positions, freed once), its trim threshold at twice that, and the batches' scratch
  freed at the top of the heap (the CSR positions, hit vectors, interval trees, sort buffers) goes back to the system
  17 times and faults in again (`brk` 52 -> 188 calls); blastn of 5,000 cDNAs 54.5 -> 57.9 MiB (dev's own builds
  52-60), x0.985; megablast of the 31,168 cDNAs and of E. coli Sakai against MG1655 unchanged (130.6 and 100.0 MiB,
  x1.002 and x1.001). At 16 and 32 threads (node 1, 6 rounds): with `-subject` QBATCHPAR's slots (removed on 2026-10-28) ran and the keep did
  not act (dc-megablast, blastn and megablast of the 31,168 cDNAs, blastn of 5,000: x1.001-1.023, cycles x1.001-1.002,
  the same peaks); with `-db` the loop keeps (dc-megablast of the 31,168 cDNAs at 16 threads `RssAnon` 138 -> 131 MiB,
  x0.993; blastn of 5,000 88 -> 87, x0.994).
- `LOSAT_X_BNHITGROUP=1` (BLASTN, speed, every task, `-subject` and `-db`; in the strict set from the pull request that
  added it (feature/manyq-reads): on unless it is `0`) the preliminary stage's collector (`blastn/blast_engine/run.rs`
  `PrelimHitCollector::x_add_grouped`) keeps one list per query, made at the first subject with an HSP, and visits only
  the queries of a subject's HSPs. `=0` runs `add_per_query`, which made those lists anew for every subject with an HSP
  and visited them all, as NCBI's `hsp_list_array` (`hspfilter_collector.c:116-161`): a 24-byte cell and a visit per
  query x subject, on the thread that takes the subjects in OID order. With reads against a database of many subjects
  (1e6 x 150 nt x 16S, 28,036 subjects, about 33,000 queries per batch) that was 53 % of the cycles at 8 threads and held
  the search's wall time. A query's hit list does not read another's and the HSPs of a query keep their order, so that the
  lists are the same (lib test `grouped_hit_lists_equal_the_lists_per_query`; `tests/blastn_hit_group.rs` compares reports
  with `=0`: `-subject` and `-db`, 1 and 4 threads, blastn-short and megablast).
  Evidence 2026-10-30 (workstation; dev 5ce40a9e0 / 13f17367e / that head with `=0` in one window, sampled VmHWM, one sha
  per case): megablast of ENA DRR1003646 reads against 16S `-db` at 16 threads, 1e6 reads 37.47 -> 5.34 s (x0.146; NCBI
  12.76 s in the same window) with VmHWM x0.9997 (537.8 -> 536.4 MiB; NCBI 5,132), 1e5 x0.174 and x1.008, 1e4 x0.520 and
  x0.999; at one thread 1e6 x0.466, instructions:u x0.359. No gain and no loss where a batch holds few queries per subject:
  QIIME-form amplicons 1e5 at 16 threads x0.996 (518.7 -> 516.6 s, one pair), VmHWM x0.9995; cDNAs against 3 Mb
  (blastn, megablast, dc-megablast, `-subject`) and one query against a genome at one thread instructions:u x0.99997-1.00003.
  The one control with a VmHWM above dev past 0.5 %: 5,000 cDNAs, blastn at 16 threads x1.010 (`=0` x1.007, the binary's placement);
  16 pairs with `MALLOC_MMAP_THRESHOLD_=131072` x0.9997.
- `LOSAT_X_BNHLGROW=1` (BLASTN, memory, every task, `-subject` and `-db`; in the strict set from the pull request that
  added it (perf/bn-hl-grow): on unless it is `0`) a hit list's array (`blastn/hsp.rs` `HitList`, the preliminary
  stage's lists per query and the final ones) is no longer reserved to NCBI's array size when a list comes
  (`kStartValue` = 100 lists, then doubled, `blast_hits.c:3219-3240`): it doubles from 4 as lists come and never holds
  more room than that size (`HitList::x_room_for_one`, `reserve_exact`). NCBI's array holds pointers (800 bytes for
  100), LOSAT's holds the lists themselves (40 bytes each), so a search of many short queries with few subjects each
  kept about 4 KB per query and per hit list (blastn-short of 30,000 reads: 114 MiB at one thread; v4-MF's arm
  MF_HLCAP). `hsplist_current`, the number `Blast_HitListUpdate` compares with the count, is NCBI's as before; only
  the capacity changes, and nothing but `push` reads it, so reports do not change (lib test
  `hit_list_array_grows_as_lists_come_within_ncbi_array_size`; `tests/blastn_hl_grow.rs` compares reports with `=0`:
  `-subject` and `-db`, 1 and 4 threads, lists past 100 subjects and lists cut by `-max_target_seqs 50`). An array
  that its owner fills without `update` is left to `push`, as before. Evidence 2026-10-26 (v4-D12 s15, workstation,
  head bfbc83c6 on dev ce281c9a; 30,000 reads of 24 letters against dmel_dna_3M, `-task blastn-short -evalue 0.001`,
  every report identical and the same as NCBI's): peak RssAnon, MiB, switch / `=0`: one thread (node 0, 2 runs) 99.2,
  99.4 / 185.4, 187.3 (NCBI 328.7); max RSS (largest VmHWM sample, node 1, 3 runs, median) `-subject` 16 threads 198.8
  / 304.6, 32 threads 334.4 / 436.8, 48 threads 433.5 / 536.4 (NCBI `-subject` 389.3, one thread: NCBI does not thread
  a subject search); `-db` (2 runs) 32 threads 296.6, 301.3 / 408.9, 391.5 (NCBI 444.3), 48 threads 396.4, 397.0 /
  493.4, 502.0 (NCBI 493.7). Time: one thread, `perf stat` `:u` against `=0`, 4 rounds rotated with an A/A arm (`=0`
  again): instructions x0.99998, cycles x0.9993 (0.986-1.012) where the A/A arm read x0.9970 (0.984-1.010), minor
  faults 94.1k -> 40.6k, system time 0.15 -> 0.07 s; 16 threads (node 1, CPUs warmed 3 s before each run), 20 pairs in
  alternating order: wall x0.961 (median; geometric mean x0.962; 20 of 20 pairs below 1), minor faults 115k -> 60k.
- BNLUTDROP (graduated 2026-10-27, in the set since 2026-10-24 (perf/mb-lut-drop, evidence at the end of this
  entry); no variable; always on) (BLASTN, every task, `-subject` and `-db`) a query batch drops its lookup tables
  when its preliminary stage ends, before its traceback (`blastn/blast_engine/x_lut_drop.rs`), where NCBI's traceback
  stage drops its table (`m_LookupTable.Reset(NULL)`, `traceback_stage.cpp:208-211`); before it graduated the tables
  lived until the batch returned. The traceback reads the HSPs that the preliminary stage kept and does not scan, and the
  tables are passed to `process_subject` by the preliminary stage's calls only (`None` in the traceback's), so the
  compiler shows that no traceback reads them; the values taken from them before the scan (the offset array size, the
  shape of the subject search) are plain numbers. A `LOSAT_X_MBLUTKEEP` table goes back to the batch loop's keep when
  the tables drop, as at the batch's end before (zeroed and taken by the next build as before); what the traceback no
  longer holds is the rest: the chains (`next_pos`, 4 bytes a position of both strands), the PV array, and every table
  that is not kept (the discontiguous tables). A batch split into query
  chunks has no table of its own (its chunks drop theirs when their calls return after their preliminary stage). With
  `LOSAT_TIMING` `[TIMING] x_bn_lut_drop: batches=` gives the number of batches that dropped their tables before their
  traceback. Evidence (workstation node 1, 16 cores with their HT siblings idle, membind 1, cores warmed 3 s before each run; dev
  a1fb63f9 against 82176dc4, 4-6 pairs in alternating order; outputs identical; peak `RssAnon` and max RSS = max(`%M`, a
  10 ms `VmHWM` sample)): megablast of E. coli Sakai against MG1655 `-subject` 135.7 -> 131.3 MiB at one thread and
  142.5 -> 136.7 at 16, `-db` 135.0 -> 130.6 and 142.6 -> 133.7 (the peak moves from the traceback to the scan; NCBI's
  anonymous memory 155.3 and 133.9, plan v4 line MF); megablast of 31,168 D. melanogaster cDNAs against dmel_dna_3M (37
  batches) at 16 threads 579.2 -> 552.1, blastn of the cDNAs 741.2 -> 705.7 and dc-megablast 392.4 -> 390.9 (16 threads,
  `LOSAT_X_QBATCHPAR`'s slots, removed on 2026-10-28, whose tables were not kept). One case reads higher: the cDNAs' megablast at one thread
  217.5 -> 219.7 (the batch loop with the keep). It is glibc's dynamic mmap threshold (freeing a mapped block of up to 32
  MiB raises it): the chains of a batch are freed before its traceback instead of after, so the threshold rises one
  traceback earlier and that traceback's medium buffers stay in the heap; with the threshold fixed
  (`MALLOC_MMAP_THRESHOLD_=131072`) the run reads 217.2 -> 217.0, and `LOSAT_X_BNLUTDROP=0` reads 217.5 as dev. Time:
  wall x0.997-1.001, cycles:u x0.991-1.000 and instructions:u x0.999-1.000 in every case.
- `LOSAT_X_BNQPACK=1` (BLASTN, every task, `-subject` and `-db`; in the strict set since 2026-10-25
  (perf/mb-query-pack, evidence at the end of this entry): on unless it is `0`) the query records read before the
  search go into a store that the batch loop takes each batch out of (`blastn/blast_engine/x_query_drop.rs`,
  `x_query_stream::Queries::Taken`), their letters packed as the reader reads them (`sequence/packed_letters.rs`: 2-bit
  codes, the runs of other bytes, the lowercase mask as intervals), made again when the loop takes the batch and
  dropped once the batch has encoded them and its query chunks have cut theirs (`record.sequence = Vec::new()`, as
  `LOSAT_X_BNSTORE` does with the subjects'); with `=0` there is no store: every record keeps its letters until the
  search ends. Since 2026-11-01 (line SF S3, family bn-query-read) this one switch holds what two did: `LOSAT_X_BNQDROP`
  (the store and the drop; in the strict set from 2026-10-25) and this switch (the packing, only for a query file of at
  least two of the task's largest batches and not for standard input); below two batches packing reads the same time
  with less memory, so small files and `-query -` are packed too. Evidence of the fold (workstation, 2026-11-01; dev
  afe173946 against the fold, 8 pairs in alternating order at one thread on node 0 and at 16 on node 1, max RSS = 2 ms
  `VmHWM`; outputs identical in every arm and equal to NCBI's `-num_threads 1`): of the files packed now, megablast of
  9.5 MB of cDNAs max RSS x0.955 at one thread (x0.979 at 16), 20 MB of them from standard input x0.922 (x0.926),
  dc-megablast of 4 and 9.5 MB x0.976 and x0.949 (x0.989, x0.954), blastn of 1 and 1.9 MB x0.988 and x0.967 (x0.991,
  x1.008), one record (E. coli Sakai) and 100 queries against 16S x0.999-1.001 (x1.001, x1.022: 70-76 MiB in every arm
  at 16 threads, the store of those queries is 0.15 MB); time x0.994-1.004 at one thread and x0.981-1.018 at 16; user
  instructions at one thread x1.0003-1.0054 (after the merge of dev 803b3c23e). Earlier evidence of the drop alone
  (perf/mb-query-letters, 2026-10-25): megablast of E. coli Sakai against MG1655 `RssAnon` 131.6 -> 126.4 MiB at one
  thread and 135.2 -> 129.9 at 16, 31,168 cDNAs at one thread 219.7 -> 215.3 MiB, time x0.98-1.00, outputs identical.
  The rule (class conditional, `x_reach!("LOSAT_X_BNQPACK", store_rule_acts)` recorded in `search_cli`): every output
  format shows no query letter (`BlastnOutputFormat::reads_query_letters`: outfmt 6, 7, 8, 9, 12 and SAM without both
  `SQ` and `SR`; BLASTN takes no custom tabular column), no `-query_loc`, no caller takes the hits' records, and the
  records are read before the search (`LOSAT_X_QSTREAM` not acting; not `run_local`, whose caller holds its records);
  the formats that show a query's letters and `-query_loc` keep the records whole. A batch split into query chunks cuts
  them from the records' letters, so it drops them after its chunks, before its traceback. The records keep their IDs,
  titles and messages, and every length read later comes from the lengths taken before (`query_lengths` in the
  traceback, the searched lengths for an invalid query's message and the Seq-aligns). With `LOSAT_TIMING`
  `[TIMING] x_bn_q_drop: batches= records= letters=` gives the batches, records and letters dropped. The packing: NCBI
  reads its batches one at a time (`GetNextSeqBatch`); LOSAT reads all records first (it needs every length before
  the first batch), so their letters, touched while the file was read, stayed with the process to the end of the
  search: dropping each batch's letters after its encoding (the drop) freed pages the heap kept for later
  allocations. Memory only; no NCBI counterpart. The reader packs each nucleotide record in place of the copy of its
  buffer (`FastaInputSource::x_pack_letters`; `U` read as `T` first, the mask kept as intervals), so no record's
  letters are held after its reading; a record keeps its letters where the store and its place in the list would take
  no fewer heap bytes than they do (a short record, many other runs; glibc's chunk sizes), and records with gap lines
  or of proteins are never packed. The batch loop makes a batch's letters again when it takes the batch
  (`Store::take`, by the search pool's threads with two threads or more), before any reader of them (WindowMasker,
  DUST, the encoding): the letters are the reader's, byte for byte. Packing reads 32 letters at a time with AVX2
  where the CPU has it (checked at run time) and 8 in a `u64` elsewhere; the letters are made 32 at a time with AVX2,
  else 4 from a table per packed byte. With `LOSAT_TIMING`
  `[TIMING] x_bn_q_pack: records= letters= bytes= unpack_s=` gives the records packed, their letters, the bytes held
  for them and the time spent making their letters again. Evidence (workstation, outputs identical; peak `RssAnon` of
  a 10 ms sample, max RSS = max(%M, VmHWM)): megablast of 31,168 D. melanogaster cDNAs (89.3 Mb, 37 batches) against
  3 Mb, `-subject` at one thread (node 0 core 28 with its sibling idle, 8 pairs against dev 83b9673b, alternating
  order, load 9-11): `RssAnon` 215.3 -> 156.7 MiB, max RSS 221.9 -> 163.5, wall x0.997 (7 of 8 faster), `cycles:u`
  x1.002, `instructions:u` x1.005, minor faults -15.0k (the reading takes 0.152 s against 0.163); again after the
  review's change (records that save no heap bytes keep their letters), under a load of 93-106: 215.3 -> 156.2 MiB,
  wall x1.002 (3 of 8), `cycles:u` x1.003, `instructions:u` x1.005; `-db` at 16 threads (node 1, CPUs 32-47, 8 pairs,
  load 20): `RssAnon` 234.4 -> 175.3 MiB, max RSS 235.5 -> 177.8, wall median x1.001 (geometric mean x0.995, 4 of 8
  faster, every pair within 0.96-1.03), `instructions:u` x1.004, `cycles:u` x1.009 (the pool's threads woken to make
  each batch's letters again: 13 ms on the loop's thread at one thread); megablast of the 5.5 Mb Sakai genome (one
  batch, not packed then): `instructions:u` x1.000. The 22.3 MB of packed letters replace 89.3 MB. Before the reader
  packed in place of its copy, packing after it with `PackedNucleotide` cost x1.020 (one thread) and x1.098 (`-db`,
  16 threads) of the time (plan v4 line A10 session 19).
- `LOSAT_X_BNQSLOTS` (REMOVED on 2026-10-28 with `LOSAT_X_QBATCHPAR`, whose slots it served; see the `LOSAT_X_QBATCHPAR`
  entry below). The `LOSAT_X_BNQPACK` store (once `LOSAT_X_BNQDROP`) no longer lends a batch to anything.
- `LOSAT_X_GREEDYROWFUSE=1` (BLASTN megablast, the distances of the non-affine greedy alignment that `x_greedy_row` runs:
  every traceback, and the score-only calls the fused loop GREEDYFUSE declines; speed; in the strict set from the pull
  request that added it, 2026-10-24; `=0` runs `x_greedy_row_avx2`, the three passes) runs one distance in two passes per
  block of 256 diagonals instead of three on x86_64 CPUs with AVX2, BMI1, BMI2 and LZCNT
  (`blastn/alignment/x_greedy_row_fuse.rs`, argument in its header). Pass 1 also keeps each diagonal's eight-byte load
  offsets and whether either lies outside its sequence; the eight-byte slide test, which pass 2 ran as a vector loop over
  the block whose lanes went back and forth between vector and general registers (line F s5, probe P6), runs as scalar
  code inside the pass over the diagonals, for the surviving diagonals only. The same cells, state and
  `x_greedy_cell_slow` calls in the same order; the unit tests compare the two builds and the plain build on random single
  distances and on whole traceback alignments. Other CPUs and targets (wasm32, aarch64, x86_64 without AVX2) run the three
  passes (`x_greedy_row_any`) as before; their machine code is unchanged, as is that of `x_greedy_row_avx2` (only the
  dispatcher `x_greedy_row` changed). With `LOSAT_X_STATS`, once a greedy row ran, `[X_STATS] GREEDYROWFUSE
  kernel=fused|avx2|plain fused_rows=N` gives the build every row of the process took and the rows of the fused build
  (counted only with the variable). Measured 2026-10-24 on the workstation against `dev` 375afb0a (outputs and max RSS
  identical in every arm; `perf stat` counts `:u`): at one thread (node 0, one pinned core, HT sibling idle, warmed 3 s; 3
  rounds, the arm order reversed every other round; Sakai x MG1655 also 16 rounds on node 1) megablast of Sakai against
  MG1655 cycles x0.900-0.908 (pairs 0.894-0.913), instructions x0.964, wall 0.623 -> 0.567 s; of Sakai with every 100th
  letter N x0.908 (instructions x0.966); of 100 16S rRNAs against the 28,036 of `16s_all` (`-db`) x0.988 (instructions
  x0.995); of 31,168 D. melanogaster cDNAs against dmel_dna_3M equal to `=0` (instructions x0.998; greedy rows are a small
  part of that search). At 8 threads (node 1, 16 rounds) Sakai against MG1655 cycles x0.915 (0.909-0.918), wall 0.282 ->
  0.274 s; at 16 threads cycles x0.920 (0.907-0.941), wall 0.259 -> 0.252 s. `=0` read x1.000-1.009 of `dev` (cycles;
  instructions x1.001-1.002: the dispatcher). Kernel (P6 on 7e73ec9c): 17.5 -> 12.3 cycles per cell of the band.
- `LOSAT_X_DISCLUTPAR=1` (BLASTN `-task dc-megablast`, and megablast given `-template_type`,
  `-template_length` and `-word_size` 11 or 12; `-subject` and `-db`; in the strict set since 2026-10-17,
  default ON: `LOSAT_X_DISCLUTPAR=0` turns it off; evidence below, "Default ON") the
  discontiguous megablast lookup tables
  (`lookup.rs` `build_disc_mb_lookup`, NCBI `s_FillDiscMBTable`: `hashtable`, `next_pos` and the counters
  of the longest chain, and with two templates `hashtable2`, `next_pos2` and their counters, both
  templates sharing one presence vector) are filled by the threads of the search pool
  (`blastn/x_disc_lookup_par.rs`). The words are read in chunks of 524,288 positions, each in three
  phases on the pool: (A) the chunk's positions are cut into one span per thread, and each thread
  reads the words of its span with the loop's statements (entering a range up to `template_length - 1`
  positions early, so the words and their template indices are the loop's) and sorts them into
  buckets by the owner of their cell (blocks of `2^max(11, pv_array_bts)` cells dealt round robin, the
  plan of PAIRPAR's megablast table); (B) each owner applies the loop's update to the words of
  its cells span by span, so each cell gets its words in the original order, and keeps each word's
  previous head; (C) each thread writes `next_pos` (and `next_pos2`) of the words of its span from the
  kept heads. Each presence-vector word is written by one owner (a presence bit is only ever set);
  each owner counts the chains of its cells (the counters of the longest chain) in a private copy of
  the counters on cache lines of its own, which the calling thread adds into the counters after the
  last chunk (the same wrapping increments of each counter by its one owner, and nothing reads the
  counters before the longest chain); and `longest_chain` is computed afterwards by the original code.
  Each bucket's header sits on a 64-byte line of its own (so the counter lines and the header lines
  are not shared between threads; the presence-vector words of two neighbouring blocks still share a
  line, which costs nothing measurable). The buckets hold one
  8-byte entry per word of the current chunk and template; before each chunk a bucket gives back the
  room it grew above four times its share of a chunk, so the buckets of one template reserve at most
  24 MiB plus 608 bytes x threads^2 (544 of entries and the 64-byte header) whatever the input; the
  private counters take threads x 8 KiB per template at word size 11 (threads x 32 KiB at 12). Measured (unit test, chunks of 524,288
  positions): 8 MiB per template on random letters at 4, 16 and 64 threads; 16 / 15 / 12 MiB on
  blocks of unmasked short repeats (each chunk's words on the cells of one or two owners, changing
  from chunk to chunk), where buckets that kept their room reserved 16 / 32 / 48 MiB (up to 2 x
  threads x 4 MiB); on (a) below 6.0 MiB at 4 threads and 7.8 MiB at 16. It acts only where a LOSAT search pool
  of two threads or more is current: the search's pool (never rayon's global pool); with one thread, a
  build without the `parallel` feature or WebAssembly without `wasm-threads` the original loop runs.
  With `LOSAT_TIMING` (any switch state) `main` prints `[TIMING] disc_lookup_build: <s>s
  (builds=<n> parallel=<n>)` when a table was built: the summed wall time of the builds (allocation,
  fill and longest chain) and how many were filled on a pool. A second line, `[TIMING]
  disc_lookup_split: alloc= loop= pool_fill= chain= chunks= read_wall= read_max= read_sum= read_late=
  update_... link_... buckets_peak= x_prefault_wall= x_prefault_max= x_prefault_sum= x_prefault_late=
  x_prefault_rounds= x_prefault_touches=`, splits the same builds: the allocation, the fill by the
  original loop, the fill on the pool (its prefault rounds included) and the longest chain (seconds,
  summed); for the pool fills their chunks and, per phase (read = A, update = B,
  link = C) summed over the chunks, the calling worker's wall time, the longest job, the jobs' summed
  time and the latest job start; the most the buckets of one fill reserved; and the same four times of
  the prefault rounds (below), their number and the zeros they wrote (one per page). The jobs' times
  are taken only with `LOSAT_TIMING`. Measured with #252's fill (before
  the false-sharing fix below) on the workstation (2x Xeon
  Gold 6338, powersave, 16 idle physical cores of one node with the siblings idle, `numactl --membind`;
  `-outfmt 6`, subject 3 Mb of D. melanogaster 2L, `LOSAT_X_PAIRPARDISC=1`; off runs carry an unread
  `LOSAT_X_PADDISCLUTPAR=1`; build times from `LOSAT_TIMING` runs, walls from runs without it),
  medians of alternating off/on pairs: (a) 1,775 cDNAs (5.0 Mb; one big query batch and two small
  ones), 16 threads: build 0.206 -> 0.089 s, wall 2.173 -> 2.050 s (on/off 0.943, 8 pairs); (b) 3,492
  cDNAs (two big batches), 16 threads: build 0.378 -> 0.158 s, wall 4.640 -> 4.433 s (0.955, 6 pairs);
  (a) at 4 threads: build 0.21 s either way, wall 2.184 -> 2.180 s (1.00, 6 pairs; in the machine's fast
  phase the build took 0.140 s without and 0.204 s with the switch). A first version in which every
  thread read every word and wrote `next_pos` itself built (a) in 0.127 s at 16 threads. Max RSS within
  the pairs 91-96 MiB on / 91-93 off for (a) at 16 threads, 104-107 / 105-107 for (b). At 4 threads
  the build was not faster. The split line on the workstation (fix round 1, node 0 CPUs 16-31, (a),
  three tables in 21 chunks, medians of 4 `LOSAT_TIMING` runs per variant; the 4-thread runs were all
  in the fast phase, 1.75-1.80 GHz, the 16-thread runs all in the slow phase, 1.28-1.38 GHz) shows
  where the time goes. 4 threads: the loop 0.204 s; the pool fill 0.224 s = phases A / B / C 0.081 /
  0.133 / 0.014 s of wall, longest jobs 0.077 / 0.130 / 0.009 s, jobs' summed time 0.238 / 0.497 /
  0.033 s. 16 threads: the loop 0.219 s; the pool fill 0.086 s = 0.017 / 0.061 / 0.009 s, longest
  jobs 0.011 / 0.057 / 0.005 s, summed 0.159 / 0.853 / 0.061 s. So the joins cost little: the latest
  job start is 0.0024-0.0033 s per phase summed over the 21 chunks, and the walls exceed the longest
  jobs by 0.011 s (4 threads) and 0.013 s (16 threads) in all, about 0.2 ms per join. The time is in
  the jobs' own run time (mostly stalls, as the cause below shows), whose sum is 3.8 times (4
  threads) and 4.9 times (16 threads) the loop's whole fill; phase B (the owners' updates) alone sums to 2.4 and 3.9 times the loop, with balanced jobs
  (longest / mean 1.04 and 1.08); at 4 threads phase A is also costly (the longest reader 1.29 times
  the mean). The runs with the switch also have more system time (0.22 against 0.15 s at 4 threads,
  0.90 against 0.60 s at 16). Fewer joins would save at most the 0.01 s above. Fix round
  1's bucket bound costs nothing measurable: (a) at 16 threads, 4 rounds old (24c714fd) / new
  (b1e3df69) / off, all slow phase: walls 2.056 / 2.055 / 2.174 s (new/old 1.000, 0.999-1.001; new/off
  0.945), max RSS 91.3-95.3 / 91.3-94.3 / 88.9-92.9 MiB; (b) at 16 threads one round in a common phase
  4.419 / 4.419 s (the other mixed the phases). Every output equals the sha256 of the run without the
  switch.
  The cause was false
  sharing (perf-mb-20261014 `logs/disclutcause_report.md`: hardware counters inside the fill and
  per-thread `getrusage` per phase, on probes of a scratch branch): in phase B the 4-byte counters of
  `helper_array`, one per block of 2048 cells with the blocks dealt round robin, put 16 owners (16
  threads) or 4 (4 threads) on each 64-byte line, and every word on a non-empty cell wrote one; in
  phase A the 24-byte bucket headers of neighbouring readers shared lines, and every push stores the
  length. That gave 3.2-4.5 million HITM loads (loads served by a line another core had modified)
  and as many RFOs (stores that had to take a line from another core) per fill of (a), 6.4-9.0
  million of each for (b), against at most a few thousand in the loop, about 80 % of them in phase B.
  The 3.8 and 4.9 times are the jobs' summed time, mostly stalls, not work: the instructions were
  1.46-1.6 times the loop's, and with all jobs on one thread the fill took 0.215 s against the
  loop's 0.213 s with the same cycles. Remote memory was no cause and the first-touch page faults
  (the loop's count, dearer with many threads at once) a secondary one (pre-faulting the tables is a
  separate, smaller follow-up).
  The fix (padded bucket headers, private counters, added once after the last chunk; the counter
  lines and the bucket-header lines are no longer shared, while the presence-vector words of two
  neighbouring blocks still share a line, which the cause report measured to cost nothing: 1.1-2.4
  million HITM loads left, no time): on the
  workstation (node 1 CPUs 32-47 with the siblings idle, `numactl --membind 1`, `LOSAT_X_PAIRPARDISC=1`;
  off = the fixed binary without the switch, #252 = dev 0b7335ce with it, fix = c33d2e2a with it; per
  input and thread count a warm-up, then rounds of the three in rotated order, 4 at 4 and 16 threads, 2
  at 2 and 8 for the fill and 6 for the wall; fill = the split line's `loop` or `pool_fill` of
  `LOSAT_TIMING` runs, wall from runs without it; every output has its sha256, (a) 749453559dfae63b,
  (b) 6c12337933ef66fe) the fill took, off / #252 / fix, (a) 0.141 / 0.196 / 0.090 s at 2 threads,
  0.171 / 0.197 / 0.083 s at 4, 0.203 / 0.127 / 0.060 s at 8, 0.201 / 0.082 / 0.049 s at 16; (b)
  0.259 / 0.366 / 0.155, 0.368 / 0.375 / 0.151, 0.374 / 0.231 / 0.100, 0.310 / 0.154 / 0.083 s (fix /
  #252 0.40-0.59, fix / off 0.24-0.64; medians over all runs). The machine's one-thread speed flips
  between runs (powersave; walls of one configuration 1.4 or 2.2 s for (a)), and the loop's fill
  follows it while the pool fills hardly do: (a) at 4 threads, off 0.139 s in the fast mode and 0.204 s in
  the slow one, fix 0.083 and 0.085 s; (b) at 16 threads, off 0.249 and 0.373 s, fix 0.083 s (slow).
  Within one speed mode (`logs/disclutfix_ws_table.md`) fix / off is 0.22-0.64 and fix / #252
  0.39-0.59 for every input and thread count, so the switch beats the loop from 2 threads on in every
  mode, and it needs no thread threshold. Walls are compared within one speed mode as well (runs
  within 10 % of the group's slowest; at 2 threads, where (a) has a fast and a slow mode, the fast
  one, and for (b) all runs): (a) at 2 / 4 / 8 / 16 threads 1.446 / 1.508 / 1.386,
  2.203 / 2.200 / 2.070, 2.190 / 2.112 / 2.044, 2.190 / 2.062 / 2.027 s; (b) 3.037 / 3.130 / 2.925,
  4.030 / 4.045 / 3.854, 4.244 / 4.173 / 4.044, 4.486 / 4.393 / 4.309 s: fix / #252 0.92-0.98, fix / off
  0.93-0.96, about the fills' differences. Max RSS medians of the arms at most 2.5 MiB above the loop's
  and within 1 MiB of #252's (88-109 MiB).
  The prefault round (2026-10-17, perf-mb-20261014 `logs/disclutfault_report.md`): phase B reads a
  cell (and a presence-vector word) before it writes it, so on a page the process had not touched the
  kernel mapped its zero page and then copied it at the write: two faults, and a TLB shootdown to the
  other CPUs of the search (on the workstation one thread touching 4,096 fresh pages took 8,193 faults
  and 13-18 ms reading first, 4,063 faults and 7 ms writing first). Now, between phases A and B of a
  table's first chunk, the readers having marked the pages of their words' cells and presence-vector
  words, one job per page (the owner of the block of the page's middle element) writes 0 to the
  page's middle element, for the pages of the whole 2 MiB regions (512 pages by address) that lie
  within the table and that the chunk writes in full; the later chunks run without marks or round
  (and the readers then run a copy of the loop without the marking code). It runs only when the first
  chunk has at least 8 positions per page of the table (`PREFAULT_RULE`). The writes store 0 over the
  0 the element was allocated with (nothing has been written before the first chunk's update), only
  on pages the fill writes anyway, so the tables and reports do not change, and neither does max RSS,
  with 4 KiB pages or with transparent huge pages: a huge page the round maps holds only pages the
  fill writes. That is the reason for whole regions: a first version touched every page the chunk
  wrote, and on a host with THP `enabled=always` (this WSL machine, with free 2 MiB blocks; the
  workstation runs `madvise`) its first writes mapped whole 2 MiB pages where the fill's first reads
  map the huge zero page and then only the 4 KiB pages written: dc-megablast of 1,000 letters at 4
  threads took 23 MiB against 17 MiB for the fill without the round at word size 11 and 134 against
  18 MiB at word size 12 with two templates, and 2.8 MiB more remained at word size 12 even with 160
  positions per page (some pages of every region stay unwritten on real sequence), so no size
  threshold alone bounded it. With whole regions, max RSS on that host equalled #263's fill (and the
  loop's, apart from the buckets) at every size measured, 1 kb to 1.3 Mb, word sizes 11 and 12
  (perf-mb `logs/disclutfault_rss_*.md`). On those test genomes no region was written in full, so the
  round touched nothing there; on the workstation's cDNA batches the first version touched every page
  of the tables, and whole regions keep part of its gains there (measured 2026-10-18, "Default ON"
  below). `next_pos` (written without a read) is not touched, which measured
  faster (1,775 cDNAs, fill over #263's: 0.79 at 16 threads and 0.91 at 4 with the cells alone, 0.89
  and 1.01 with `next_pos` as well). The timing evidence is Linux only: the saving rests on Linux's
  shared zero page and copy-on-write; on other systems the round is exact and maps no extra memory,
  but its gain is unmeasured. The unit tests check, for every partitioned build, that each zero went
  over 0, that there was one round at most, and that the pages touched are exactly the first chunk's
  pages of the accepted regions, array by array.
  Default ON (2026-10-17, Owner rule of 21:50: speed-ups land in the strict set in their PR). On the
  workstation (node 1, 16 idle physical cores with the siblings idle, `numactl --membind 1`; the rest of
  the dsw `cand` env, `LOSAT_X_PAIRPARDISC`, `QBATCHPAR` (removed on 2026-10-28), `DIAGLV`, `PAIRPARFOLD` and `QDUSTPAR`, in
  every arm; off = the head with `LOSAT_X_DISCLUTPAR=0`, dev = 7a5a7f87 with the switch (#263's fill),
  head = 1a05a8a2 by default, the first version of the round; rotated rounds after a warm-up, 4 at 2 and 4 threads, 3 at 8 and 16, 2
  at 1; every output had its sha256) the fill took, off / dev / head: (a) 1,775 cDNAs 0.140 / 0.089
  / 0.080 s at 2 threads, 0.170 / 0.081 / 0.071 at 4, 0.205 / 0.066 / 0.055 at 8, 0.204 / 0.051 /
  0.043 at 16; (b) 3,492 cDNAs 0.264 / 0.151 / 0.149, 0.287 / 0.168 / 0.152, 0.354 / 0.149 / 0.142,
  0.375 / 0.091 / 0.080; 31,168 cDNAs (37 batches, summed over them) 2.19 / 1.21 / 1.24 s at 4 threads
  and 2.72 / 1.13 / 1.11 at 16; 3 cDNAs (5 kB) 0.011 / 0.010 / 0.004 at 4 and 0.011 / 0.008 / 0.003
  at 16. The search took head / off 0.90-0.93 (2 threads), 0.89 / 0.93 (4), 0.76 / 0.86 (8), 0.68 /
  0.69 (16) for (a) / (b), 0.95 and 0.85-0.88 for 31,168 cDNAs at 4 and 16 threads, 0.86-0.90 and
  0.87 for 3 cDNAs, and 1.00 at 1 thread, where the loop runs either way (head / off 1.001 and
  1.005); against dev 0.96-0.99 at 8 and 16 threads and 0.97-1.06 at 2 and 4 (inside the rounds'
  spread, 0.99-1.12). So the switch wins wherever it runs and needs no
  thread or size rule. The prefault cut the TLB shootdowns received by the search's CPUs from 108,000
  / 110,000 / 301,000 to 100 / 1,000 / 5,300 per run at 16 threads ((a) / (b) / 31,168 cDNAs), the
  minor faults by 7,400-33,000 and the system time by 0.15-0.23 s. Max RSS of head within the run's
  spread of dev's and off's (-7.6 to +6.6 MiB at 31,168 cDNAs, at most 3.1 MiB elsewhere) and below
  NCBI BLAST+ 2.17.0 `-db` of an earlier window (31,168 cDNAs 305 against 437 MiB at 4 threads and 497
  against 700 at 16; (b) 176 against 199 at 4). At 31,168 cDNAs and 4 threads the fill of 1a05a8a2 was
  1.02 of dev's, its open loss: its reader loop with the mark tests took 7-11 % longer in phase A in
  every cell, also where nothing was marked; the current code runs the marking loop only in a
  table's first chunk and a copy without the marking code elsewhere, which removed that loss. The final
  form (whole regions, the size gate and first-chunk marking; 2026-10-18, the same setup and arms with
  head = that code; perf-mb `logs/disclutfault_ws_table_e7143dfc.md`) took, fill off / dev / head: (a)
  0.136 / 0.086 / 0.084 s at 2 threads, 0.138 / 0.068 / 0.066 at 4, 0.203 / 0.061 / 0.059 at 8, 0.202 /
  0.050 / 0.050 at 16; (b) 0.257 / 0.152 / 0.148, 0.283 / 0.168 / 0.153, 0.290 / 0.128 / 0.122, 0.364 /
  0.085 / 0.084; 31,168 cDNAs 2.18 / 1.22 / 1.20 s at 4 threads and 2.70 / 1.15 / 1.08 at 16; 3 cDNAs
  0.011 / 0.010 / 0.011 at 4 and 0.011 / 0.008 / 0.007 at 16 (the size gate keeps the round from
  running there). Head over dev per round for (a), (b) and 31,168 cDNAs: fill 0.91-0.99 at 2-8 threads
  and 0.93-1.00 at 16 ((a) at 16 threads 1.00, rounds 0.92-1.04), search 0.97-1.02 (inside the spread),
  max RSS -2.9 to +2.3 MiB. TLB shootdowns per run at 16 threads, dev / head: 108,000 / 68,000 ((a)),
  111,000 / 73,000 ((b)), 261,000 / 112,000 (31,168 cDNAs). Whole regions leave the pages of partly
  written regions to the fill's own faults, so the final round removes fewer shootdowns than the first
  version (100 / 1,000 / 5,300) and gains less over #263's fill (the first version: 0.83-0.99, 1.02 at
  31,168 cDNAs and 4 threads). Apart from 3 cDNAs at 4 threads, fill 1.07 per round (0.7 ms; the round
  does not run there; search 0.94), it loses nowhere measured, so it needs no rule. Default ON for the switch itself
  rests on its margin over the loop (the previous default), which no version of the round changes much. Tests: the
  unit tests of `x_disc_lookup_par.rs` (partitioned and pool fills equal the original loop for every template type,
  length and word size, one and two templates, 1-20 parts, chunks of 1 to 5,000 positions and random
  layouts of queries and masks, pools of 1, 2, 3, 5, 8 and 16 threads, a 2.6 M-position query in the
  default chunks on pools; the template indices read only the last `template_length` bases; the owners
  run on at least two workers of a 3-thread pool; a test names a change in any one array or the
  longest chain; the buckets' bound; the private counters give the loop's counters on pools of 2, 3, 5
  and 16 threads also from start values that the increments wrap past u32::MAX, with blocks whose
  counters wrap and blocks whose counters stay zero, and at least two owners' private copies hold
  counts before the flush (a fill that counted into one shared copy would be just as exact but bring
  the false sharing back); the bucket headers and the counter copies are
  whole 64-byte lines) and `LOSAT/tests/blastn_disc_lookup_par.rs`
  (off / on equal at 1 and 4 threads, `-subject` and `-db`, one and two templates, with
  `LOSAT_X_PAIRPARDISC` and split query chunks, megablast with and without a
  template, the timing and split lines, and AP027131 against NZ_CP006932 with one and two templates,
  whose tables span two chunks).
- `LOSAT_X_DCCSR=1` (BLASTN `-task dc-megablast` with one template, `coding` or `optimal`, and megablast
  given `-template_type`, `-template_length` and `-word_size` 11 or 12; `-subject` and `-db`, any thread count;
  in the strict set from the pull request that added it, 2026-10-21: `LOSAT_X_DCCSR=0`
  turns it off) lays out the first template's chains of a discontiguous megablast table as CSR
  (`blastn/x_dc_csr.rs`, argument in its header): `offsets[cell] .. offsets[cell + 1]` of `qs` holds the
  cell's 1-based query offsets in the chain's order (newest first), built from the fill's own words by a
  count pass, a running sum (which also makes the fill's presence vector and longest-chain counters from
  the counts) and a place pass, instead of NCBI's `hashtable` / `next_pos` (`s_FillDiscMBTable`), which
  are then not allocated: `qs` has one element per word where `next_pos` had one per query position, so
  the table takes no more memory. Only tables whose word loop reads at least a quarter as many positions
  as the table has cells (1,048,576 at word size 11) get it; smaller ones, such as NCBI's one-query
  batches, keep the chains, which fit the caches there. Two-template tables (`coding_and_optimal`) keep
  their chains. With a search pool of two threads or more (where `LOSAT_X_DISCLUTPAR`'s fill would run)
  the table is built on the pool (`x_disc_lookup_par::csr`): the fill's readers put the words into the
  owners' buckets, once to count and once to place, and each thread owns a contiguous run of cell blocks;
  the whole build is one `rayon::broadcast` whose threads pass a spin barrier between the reading and the
  owners' work of each chunk (one `rayon::scope` per phase and chunk started its last job about 0.15 ms
  late each time; v4-R R16). The scans read a word's run instead of walking its chain: the strip scan of
  `LOSAT_X_PAIRPARDISC` and of its diagonal-array stage, and the reference scan (one thread with `LOSAT_X_BNDCONE=0`, `-db` with
  subjects searched by several threads, and the cases outside those stages) put the words through a FIFO
  of 32 words that prefetches each word's `offsets` element on entry and the first and last cache line of
  its run 16 words later (x86_64; elsewhere the FIFO only delays the words); the words, pairs and their
  order are the chain walk's. The masked-word test of `s_TypeOfWord` (`contains_hit`) searches the run.
  The `[TIMING] disc_lookup_build` interval includes the CSR builds (in `loop` or `pool_fill` of
  `disc_lookup_split`; its phase times count the chain fills only), and `[TIMING] disc_lookup_csr:
  builds= parallel= words=` counts them. The unit test `the_runs_are_the_chains_with_the_fills_presence_vector_and_counters`
  compares every cell's run with its chain (and that it decreases), the presence vector, the counters and
  the longest chain, and a pool build with the serial one.
  Evidence (ws, node 1 timed, dev 7bea6267 against perf/dc-csr 8b8b3f51, medians; v4-A15 ab2 / ab3): 31,168
  dmel cDNAs against 3 Mb of dmel genome with `-subject`, 16 threads 2.619 -> 1.847 s (x0.705, cycles:u 101.7 ->
  69.1 G), 8 threads x0.733, 1 thread x0.788; the same queries against that genome as `-db`, 16 threads x0.796;
  1,775 cDNAs x0.86 (16 threads) and x0.80 (1 thread); 236 cDNAs (one table at 0.31 of the cells, just above the
  size rule) x0.92-0.96 at 1, 8 and 16 threads; 100 16S queries against a 16S `-db` (small tables, chains in
  both arms) x1.003. `=0` is dev's (x0.998); instructions of megablast and blastn at one thread move by at most
  0.04 %. Peak RSS with `MALLOC_ARENA_MAX=1`: 374 against 386 MiB (`-subject`, 16 threads) and 113 against
  115 MiB (`-db`); with the default arenas both arms vary between about 130 and 144 MiB on `-db`.
- `LOSAT_X_DCSEEDLEAN=1` (BLASTN `-task dc-megablast`, and megablast given template options, at two threads or more;
  in the strict set from the pull request that added it, 2026-10-28: `LOSAT_X_DCSEEDLEAN=0` turns it off) changes
  three LOSAT-only parts of dc-megablast's seed stage; each gives the same table, seeds and hits in the same order:
  - the CSR of a one-template table built on a pool (`LOSAT_X_DCCSR` above; `x_disc_lookup_par::csr_on_pool_steps`):
    two sets of the owners' buckets, chunk k's in set k % 2. Step k of each pass (count, then place) reads chunk k into
    its set with the lean reader (`Fill::read`'s word loop with the template a const parameter, so a word's index is
    computed without the dispatch per word of `compute_discontiguous_index`, and a range read through one slice), then
    does the owner's work on chunk k - 1 from the other set; one spin barrier ends the step. `csr_on_pool` waits twice
    per chunk (after the reading and after the owners' work), so each chunk costs the slowest reading plus the slowest
    owner's work; a step costs the slowest thread's reading and owner's work together (about as many barriers per pass:
    chunks of half the size, one barrier each). Chunks are 2^18 positions (`STEP_POSITIONS`, half of
    `CHUNK_POSITIONS`, so the two sets keep room for about as many entries as one set did); the place pass prefetches the
    `offsets` element of the entry 16 ahead, as the count pass does. A cell gets its words in the same order (chunks
    in order, readers' spans in stream order, a bucket in reading order), so the counts, sums and runs are
    `csr_on_pool`'s.
  - the partitioned stage of `LOSAT_X_PAIRPARDISC` (`LOSAT_X_PAIRPARFOLD`) in `LOSAT_X_PAIRPARHOME`'s home loop: a
    strip's scan pushes each seed into its partition's list (`x_pair_seed_part::PartSink`; the partition of the seed's
    diagonal bucket from a table of the 512 buckets made at the stage's start), numbered in scan order, instead of into
    one list that a stable counting sort then cuts by partition (two more passes over the strip's seeds and a copy of
    each 16-byte seed and its 4-byte position). A partition gets the same seeds in the same order with the same
    positions, so its folds and the merge's key order are the same.
  - the folds of that stage (home loop or not) do not record the buckets their inserts change (`insert_t` without
    `TRACK`), and the merge writes back every bucket of each partition: a bucket no fold changed has the stage's chain,
    cells below the stage's occupancy, so writing it back changes nothing. That is one copy of the table's cells per
    stage, where `rebuild` copies the table into each partition at the stage's start; the stages of dc-megablast begin
    with a nearly empty table (at most 256 cells in the 27 stages of the 31,168 cDNAs below).
  The scan into lists and the untracked folds are dc-megablast's only: the partitioned stages of `-task blastn` and
  `blastn-short` (`LOSAT_X_PAIRPAR`, `LOSAT_X_PAIRPARNA`) keep the cut and the tracked folds. `x_reach!` at each pool
  build and each partitioned stage, always the switch's side while it is on. Tests: `csr_on_a_pool_lists_the_chains`
  (`csr_on_pool` and `csr_on_pool_steps` on pools of 2, 3 and 5 threads, chunks of 1, 97, 5,000 and `CHUNK_POSITIONS`
  positions, against the chains of the serial loop), `seed_stage_part_equals_the_reference_fold` (the stage with and
  without the lists and untracked folds, in every pairing with the home loop for each level setting, one scratch
  passing its strip buffers from one form to the other, against the reference fold's hits and table cell for cell),
  the dc-megablast stage test of `x_pair_seed_disc.rs` (both forms in both `LOSAT_X_FOLDLEAN` instances against the
  reference fold) and `LOSAT/tests/blastn_dc_seed_lean.rs` (the reports at 2 and 4 threads, and at 4 with
  `LOSAT_X_PAIRPARHOME=0`, against `=0` at 2 and 4: Sakai against MG1655, a CSR table built on the pool, and 9,800
  letters of LC738874 against LC738868, a table of chains; with the sides of `LOSAT_X_REACH`, of `LOSAT_X_DCCSR`,
  `LOSAT_X_DISCLUTPAR` and `LOSAT_X_PAIRPARFOLD` too).
  Probe (v4-A15 s15-s18, 2026-10-26/28, exp/a15-csrpool; workstation 2x Xeon Gold 6338, node 1 booked through its FIFO,
  16 / 8 physical cores with the siblings idle, `numactl --membind 1`, the arms in rotated order, 16 rounds at 16
  threads; dc-megablast of 31,168 D. melanogaster cDNAs against 3 Mb of chromosome 2L, `-subject`, every output NCBI's
  5c357f81): dev's form 1.441 s, 147.07 G instructions, 58.39 G cycles, max RSS 237 MiB; the lean reader, the scan into
  lists and the untracked folds x0.957 (cycles x0.946, 133.04 G instructions, 221 MiB); with the steps as well (this
  switch) x0.951 (median x0.952, no round slower, cycles x0.939, 137.12 G instructions, 220 MiB), at 8 threads x0.946.
  Of the parts (s16, the same window): the untracked folds -2.7 G instructions and -0.5 G cycles; the scan into lists
  -4 G instructions and -15 MiB (the strip buffers' copies) but about the same cycles (its stores go to buffers other
  cores read in their folds); the steps against the lean reader alone -0.45 G cycles of the build (the waits of two
  barriers per chunk, each waiting for the slowest thread, not the overlap of reading and owner's work: s15's form
  taking one entry after each word read gained no more). Chunks of 2^19 positions took 55.42 against 54.75 G cycles
  (screening on node 0). One thread takes none of these paths (no pool, no partitioned stage): 117.68 G instructions as
  dev's.
  Pull request A/B (v4-A15 s19, 2026-10-28; 7466dbc1b against dev 1ad262297, a byte copy of dev's binary (A/A) and the
  pull request's binary with `=0`; node 1 as above at 16 threads unless noted, every arm's outputs the same, NCBI's
  5c357f81 for the 31,168 cDNAs; median of the per-round ratios, pull request / dev [dev's copy / dev]): dc-megablast
  of the 31,168 cDNAs x0.951 (24 rounds, none slower; 1.418 -> 1.348 s, CPU x0.947, max RSS 235 -> 220 MiB) [0.997], at
  8 threads x0.944 (12 rounds, none slower; 2.518 -> 2.377 s, 212 -> 200 MiB) [0.998], against the 3 Mb as `-db` x0.953
  (12, none slower; 129 -> 120 MiB) [1.000]; of 1,775 cDNAs x0.951 (40, none slower; 98 -> 93 MiB) [1.000]. The small
  stages, inside the A/A spread: 236 cDNAs (two stages, 1.37 M seeds, 53 hits; 45 ms) 0.997 [1.001] and 7 cDNAs (no
  hit) 0.995 [0.992]; with `LOSAT_TIMING` (node 0, 8 runs per arm) the merge of the 236 cDNAs takes under 1 ms with and
  without the switch, and the folds 0.033-0.040 s against 0.040-0.043 s. Other tasks (the switch off their paths):
  megablast of the 1,775 / 31,168 cDNAs 0.998 / 0.973 [0.999 / 0.968], the 31,168 with `LOSAT_X_QBATCHPAR=0` 1.006
  [0.986]; blastn of the 1,775 / 5,000 cDNAs 1.001 / 0.999 [0.996 / 1.002], the 31,168 with `LOSAT_X_QBATCHPAR=0` 1.005
  [1.008]; blastn-short of the 236 cDNAs 1.001 [0.997]; max RSS of these within the A/A spread. `=0` against dev:
  0.993-1.005 in every case. One thread (node 0, `perf stat`, outputs the same): instructions of dev, the switch and
  `=0` within 0.01 % (31,168 cDNAs 117.750 / 117.751 / 117.749 G; 1,775 cDNAs dc-megablast, megablast and blastn too).
- `LOSAT_X_BNSEEDLISTS=1` (BLASTN `-task blastn`, `megablast` and `blastn-short`: the partitioned stages of
  `LOSAT_X_PAIRPAR` and `LOSAT_X_PAIRPARNA` with `LOSAT_X_PAIRPARFOLD`, two threads or more, one unmasked subject, a
  query block above 8,000 letters; in the strict set from the pull request that added it, 2026-10-28:
  `LOSAT_X_BNSEEDLISTS=0` turns it off) gives those stages the form `LOSAT_X_DCSEEDLEAN` gives dc-megablast's
  (`x_pair_seed_part::seed_stage_part_lists`); each partition gets the same seeds in the same order, so the hits and
  the table are the same:
  - in the home loop (`LOSAT_X_PAIRPARHOME`) a strip's scan pushes each seed into its partition's list (`PartSink`)
    instead of into one list that a stable counting sort then cuts by partition (`produce`: a count pass, the
    partition of each seed's bucket, and a copy of each 16-byte seed and its 4-byte position);
  - the folds (`fold_t`, `fold_na_t` without `TRACK`) do not record the buckets their inserts change, and the merge
    writes back every bucket of each partition (an unchanged one as it was).
  The scans (`scan_strip`, `scan_strip_csr`, `scan_strip_na`) take the receiver as a type parameter (`SeedOut`), so
  the lists are one more instance of each; the stages without partitions and one thread are unchanged. `x_reach!` at
  each partitioned stage, always the switch's side while it is on. Tests: `seed_stage_on_strips_equals_one_scan` of
  `x_pair_seed.rs` (the partitioned stage in both forms on 1 to 8 threads and 2 to 16 partitions, on tables with and
  without level bounds, two chunks in a row on one scratch, against the serial fold) and of `x_pair_seed_na.rs` (word
  sizes 6, 7 and 8, both forms on 2 to 6 threads), and `LOSAT/tests/blastn_bn_seed_lists.rs` (blastn and megablast of
  LC738874 against LC738868 and blastn-short of 700 24-mers of LC738868: the reports at 2 and 4 threads and at 4 with
  `LOSAT_X_PAIRPARHOME=0` against `=0` at 2 and 4, with the sides of `LOSAT_X_REACH`, `LOSAT_X_PAIRPARFOLD` and the
  stage's switch).
  Probe (v4-A15 s20, 2026-10-28, exp/a15-bnlists; workstation as above, 16 threads unless noted, against the 3 Mb of
  chromosome 2L): the cut took 7-8 % of the cycles of blastn-short of 236 D. melanogaster cDNAs (three stages, 313 M
  seeds) and about 2 % of blastn of 31,168 and 5,000 cDNAs (216 M and 35 M seeds; `perf record`, node 0). Time A/B on
  node 1 (dev's form, a byte copy of its binary, this form; every output the same; median of the per-round ratios
  [A/A]): blastn-short x0.972 (16 rounds, none slower; max RSS 408 -> 385 MiB) [1.002], at 8 threads x0.960 [1.000];
  blastn of the 31,168 cDNAs x0.988 [0.999], against the 3 Mb as `-db` x0.991 [1.001] (104 -> 95 MiB), with
  `LOSAT_X_QBATCHPAR=0` x0.993 [0.996] (101 -> 93 MiB); of 5,000 cDNAs x0.988 [1.002]; of 1,775 x0.998 [0.999].
  Pull request A/B (v4-A15 s21, node 1: dev 2087618e6, a byte copy of its binary, d6e69cd1f, and d6e69cd1f with `=0`;
  every output the same; median [95 %], rounds slower, [A/A]): blastn-short of the 236 cDNAs x0.969 [0.966-0.974]
  0/16 [0.997] (2.359 -> 2.285 s, max RSS 416 -> 386 MiB), at 8 threads x0.961 0/8 [1.000] (383 -> 373 MiB); blastn
  of the 31,168 cDNAs x0.986 1/16 [0.999] (392 -> 383 MiB), as `-db` x0.996 [1.002] (103 -> 94 MiB), with
  `LOSAT_X_QBATCHPAR=0` x1.003 8 rounds [1.005] and x0.998 16 rounds [1.000] (102 -> 92-96 MiB); of 5,000 x0.996
  [1.001]; of 1,775 x0.994 [0.998]; megablast of 1,775 x0.994 40 rounds [0.989], of the 31,168 x0.999 [0.985];
  dc-megablast of the 31,168 (not reached) x0.997 [0.999]. Max RSS rises in no case. `=0` against dev 0.996-1.005.
  The one case it does not help: megablast of the 31,168 cDNAs with `LOSAT_X_QBATCHPAR=0` at 16 threads, x1.020 in 12
  rounds and x1.012 [0.998-1.024] in 40 (25 slower; A/A 0.999; CPU x1.013 in both), at 8 threads x1.004 [0.980-1.021]
  (A/A 0.992). Its partitioned stages are small: 25 stages, 131 k seeds, 15 seeds a strip, so about one per partition,
  which the lists put in 16 vectors per strip where the cut wrote one; `[TIMING]` fold and merge the same,
  instructions +0.1-0.2 % at 16 threads and +0.02 % at 8, page faults the same (both arms 70-150 k, run to run), sys
  time +0.03-0.05 s (node 0, `perf stat`). One thread (node 0, `perf stat`, outputs the same): instructions of dev,
  the switch and `=0` within 0.01 % (blastn-short of the 236 cDNAs 260.834 / 260.821 / 260.821 G; 1,775 cDNAs
  blastn, megablast and dc-megablast too).
- `LOSAT_X_BNCSR=1` (BLASTN `-task blastn` and every task with word size 9, 10 or 11 that builds the contiguous
  two-stage table of `build_mb_lookup`, `-subject` and `-db`, any thread count, speed; in the strict set from the pull
  request that added it, 2026-10-25 (perf/bn-csr, evidence at the end of this entry): on unless it is `0`) lays out the
  chains of that table as CSR, as `LOSAT_X_DCCSR` does for a discontiguous table (`blastn/x_bn_csr.rs`, argument in
  its header): `offsets[cell] .. offsets[cell + 1]` of `qs` holds the cell's 1-based query offsets in the chain's
  order, newest first for the megablast-type table (word size 11's) and in reading order for the table that stands
  for NCBI's small table (`ascending_cells`, its count kept two elements later so that the place pass moves each
  cell's cursor up from the run's start), built from the fill's own words by a count pass, a running sum (which also
  makes the fill's presence vector and longest-chain counters) and a place pass, instead of NCBI's `hashtable` /
  `next_pos`, which are then not allocated (`qs` has one element per word where `next_pos` had one per query
  position). Only tables of `lut_word_length` 8 to 11 without database word counts and without `BLEMIR_DEBUG` whose
  word loop reads at least a quarter as many positions as the table has cells get it: with NCBI's choice of table, a
  query batch of 1,048,576 positions or more (about 524 kb, both strands) at word size 11, 262,144 at 10, 65,536
  (and 16,384 for the reading-order table of 4^8 cells) at 9; word sizes 12 and above build 4^11-cell tables only
  below that size and keep their chains, and so do megablast's 4^12-cell tables below `LOSAT_X_MBCSR`'s rule. Such a table takes nothing from the
  batch loop's keep (`LOSAT_X_MBLUTKEEP`): its `offsets` is allocated per build and freed with the table (before the
  traceback, BNLUTDROP), and a build laid out as CSR frees the table the keep holds. With a search pool
  of two threads or more (where PAIRPAR's fill would run) the table is built on the pool
  (`x_mb_lookup_once::csr`): the readers of `LOSAT_X_MBLUTONCE` put each chunk's words into the owners' sub-ranges of
  its arena, once to count and once to place (so the build holds one chunk's words, the fill's arena, not the whole
  query's), each thread owns a contiguous run of cell blocks, and the whole build is one `rayon::broadcast` whose
  threads pass a spin barrier between the reading and the owners' work of each chunk. The scans read a word's run
  instead of walking its chain: the strip scan of PAIRPAR (`scan_strip_csr`) and the reference scan of
  `search_query_batch` put the words with hits (after the presence test) through the FIFO of `LOSAT_X_DCCSR`; the
  words, pairs and their order are the chain walk's. `for_each_hit`, `count_hits` and the masked-word test of
  `s_TypeOfWord` (`contains_hit`) read the run where the table has no chains. The `[TIMING] x_mb_lookup_build`
  interval includes the CSR builds (in `loop` or `pool_fill`), and `[TIMING] x_bn_csr: builds= parallel= words=`
  counts them. The unit tests compare every cell's run with the loop's chain, in both orders and for `lut_word_length`
  8 to 11, the presence vector, the counters and the longest chain, a pool build in both orders (chunks of 1, 97
  and 5,000 positions and the default size, 2 to 40 threads) with the serial one, and a table above the rule built by
  `build_two_stage_lookup` with its chains; `LOSAT/tests/blastn_bn_csr.rs` compares whole reports of the default with
  `=0`. Evidence (workstation, 2026-10-25, tree c14f064f on dev 105e4a60, medians of per-round ratios to `=0`, every
  output equal and AB's equal to NCBI's): 16 threads on node 1, 31,168 dmel cDNAs × dmel 3 Mb `-task blastn`
  `-subject` ×0.854 (0/10 rounds slower, CPU ×0.844), 5,000 cDNAs ×0.877 (CPU ×0.845), 1,775 cDNAs ×0.966 (3/30),
  at `-word_size 9` ×0.768 and 10 ×0.921, the 31,168 against a `-db` ×0.958 and with `LOSAT_X_QBATCHPAR=0` ×0.952,
  megablast (no CSR) ×0.979 (7/16, noise); one thread on one node 0 core (perf stat, user cycles) ×0.80 to ×0.88 at
  1.005 to 1.03 times the instructions; max RSS (time %M with a 2 ms sampler, three runs per arm) not above `=0`'s
  except the 1,775 cDNAs' 82 → 85 MiB (dev 84): 31,168 cDNAs 721 → 712, word size 9 329 → 322 MiB.
- `LOSAT_X_BNCSRPACK=1` (BLASTN, the CSR tables of `LOSAT_X_BNCSR` whose runs list their words newest first: every
  such table but the reading-order one of 4^8 cells, `-subject` and `-db`, any thread count, memory; in the strict set
  from the pull request that added it, 2026-10-29 (perf/bn-csr-pack2, evidence at the end of this entry): on unless it
  is `0`) packs the head of the table (`blastn/x_bn_csr.rs`, `x_pack_on`; `x_dc_csr::Csr::x_bits`): two words per 64
  cells (the cells' presence bits, then the number of non-empty cells before them) and one `offsets` element per
  non-empty cell and one more, indexed by that rank, instead of one element per cell. A table of 4^11 cells has
  4,194,305 offsets (16 MiB) per build, of which a query batch of about 1 Mb of `-task blastn` fills 1.2-1.6 M; the
  packed head is 1 MiB of bits and ranks and 4.6-6.1 MiB of offsets. A packed table has no size rule: every
  newest-first table of `lut_word_length` 8 to 11 is laid out as CSR and packed, also below `LOSAT_X_BNCSR`'s rule
  (the first and last query batch of a run, a run of one small batch, megablast's tables of 4^11 cells or fewer).
  Below the rule those batches built chains on the batch loop's kept table (`LOSAT_X_MBLUTKEEP`), whose pages past the
  start a packed head writes were then first mapped by the chains' zeroing at the end of the run (16 threads, 5,000
  cDNAs: max RSS 2 MiB above dev at the run's last 1-3 % while the live data was 12 MiB below it, read with the mmap
  threshold pinned in a diagnostic run, s27 w21). The table's form is chosen once per table (`Csr::x_packed`):
  `Csr::run` and `Ring` are the code of a head of one element per cell, used by every discontiguous table and by a
  contiguous one with the switch off, and the contiguous scans (`scan_strip_csr`, the reference scan of
  `search_query_batch`) are generic over `x_dc_csr::XHead`, so each form's scan is compiled on its own (the
  dc-megablast and `=0` scans have dev's instructions). The packed ring (`XRingPacked`) drops a word of an empty cell
  (its run is empty, so the emit would read nothing) and prefetches the rank's two offsets; the one-word answers of
  the lookup (`for_each_hit`, `count_hits`, `contains_hit`) read either form (`Csr::x_run_any`). The build counts each
  word into one byte per cell (`x_count8`: up to 255, the words past 255 of a cell in a map of the thread that owns
  the cell; 4 MiB for 4^11 cells), takes each block's presence bits from its 64 bytes (`x_presence8`, eight bytes at a
  time without a branch), writes the blocks' ranks and the counts by rank (`x_rank_counts8`), then the sum pass
  (`x_sum_packed`: the running sums of the non-empty cells, and the presence vector and helper counters that
  `sum_cells` gives over the same cells) and the place pass by rank. The bytes are the first `cells` bytes of the run
  list's capacity (`qs`, which has room for every position the word loop reads and for the bytes), zeroed by the
  cells' owners; the place pass writes the run list only after the last count is read, so the build allocates nothing
  more and frees nothing in the middle of a batch. The serial build (`x_serial_packed`) and the pool build
  (`x_mb_lookup_once::csr_on_pool`: each owner counts, marks and ranks its blocks, one barrier to publish its number
  of non-empty cells and one after the ranks) read the words twice, as the build of one offset per cell. With
  `LOSAT_X_MBLUTKEEP` the kept table serves as the packed offsets (a fresh packed build hands the keep a zeroed table
  of `cells + 2` elements, which a later chains build of the loop takes), and the keep's zeroing writes every cell of
  a chains build, also those no packed build wrote (`Keep::take`; its entry). The runs, their order, the presence
  vector and the counters are those of the unpacked head, so every value read from the table is the same. Not packed:
  the reading-order table (4^8, `ascending_cells`; it keeps the size rule), the discontiguous tables of
  `LOSAT_X_DCCSR` and megablast's 4^12-cell chains (below `LOSAT_X_MBCSR`'s rule). The CI fixtures do not reach the reading-order table (it needs
  `-word_size 9` with a query batch of 16,384-21,000 positions), so its reach row is a gap. The tests of the chains'
  fill (`LOSAT/tests/blastn_mb_fill.rs`, `blastn_mb_lookup_once.rs`, `blastn_mb_lut_keep.rs`) run `-task blastn` with
  the switch off, which keeps the small tables' chains. Tests:
  `x_dc_csr::tests::{a_packed_head_gives_the_same_runs_and_the_ring_skips_empty_cells,
  byte_counts_give_the_presence_ranks_and_exact_counts}`, `x_bn_csr::tests::the_packed_head_lists_the_chains_too`
  (packed and unpacked against the loop's chains, presence vector, counters and longest chain, `lut_word_length` 8-11,
  with cells of 255, 256 and 1,200 words), `x_mb_lookup_once::tests::csr_on_a_pool_is_the_serial_build` (the pool
  build in chunks of 1 to 2^19 positions on 2 to 40 threads against the serial one, with those cells in half the
  rounds) and `LOSAT/tests/blastn_bn_csr.rs` (the whole search with `LOSAT_X_BNCSRPACK=0`: the same reports, and both
  reach sides: word size 11 packed, word size 9 the reading-order table). Evidence (workstation 2x Xeon Gold 6338,
  2026-10-29, perf/bn-csr-pack2 bebbad34 on dev ee39301d: three arms in one window, dev, the branch and the branch
  with `=0`, order rotated per round, medians of per-round differences and ratios, one output per cell; D.
  melanogaster cDNAs (AB 31,168, C5 5,000, C4 4,000, Q300 / Q500 the first 107 / 180 of them, one query batch below
  `LOSAT_X_BNCSR`'s rule) against the first 3 Mb of chromosome 2L, DB the same 5,000 cDNAs against those 3 Mb as
  `-db`, MB megablast of the 31,168; node 1 through its FIFO, `numactl --membind 1`, peak RssAnon and VmHWM sampled
  every 2 ms): 16 threads (8 rounds; C5, DB, Q300, Q500 10), peak RssAnon dev -> branch in MiB with the median of
  per-round differences and the rounds above dev, then wall: AB 100.5 -> 89.8 (-10.7, 0/8) x0.991, AB at 8 threads
  69.6 -> 55.7 (-13.1, 0/8) x0.990, at 32 125.6 -> 116.6 (-9.0, 0/8) x0.993, C5 84.6 -> 72.7 (-11.7, 0/10) x0.985, at
  8 60.2 -> 48.4 (-11.9, 0/8) x0.987, at 32 97.1 -> 86.0 (-10.6, 0/8) x0.997, C4 80.8 -> 71.7 (-8.3, 0/8) x0.978, DB
  82.0 -> 73.4 (-9.0, 1/10) x0.983, Q300 32.8 -> 17.9 (-14.9, 0/10) x0.955, Q500 37.6 -> 21.8 (-15.8, 0/10) x0.977, MB
  169.5 -> 167.6 (-1.0, 4/8) x0.983; `=0` against dev -3.7 to +2.9 MiB and x0.997-1.023; control dc-megablast (AD, the
  31,168 cDNAs) -4.6 MiB, x0.994, instructions x1.0001. One thread (node 0 CPU 22, sibling idle, same three arms, 8
  rounds, perf stat): AB 57.7 -> 47.0 MiB (-10.7, 0/8) x0.993, C5 53.5 -> 40.8 (-12.7, 0/8) x0.990, Q500 30.5 -> 17.0
  (-13.4, 0/8) x0.900; MB 132.4 in every arm (the branch 4 KiB above dev: its peak is in a batch of 4^12-cell chains,
  which the switch does not change) x0.999; AD instructions x1.0000 of dev (`=0` x1.0000), wall x0.993. Earlier forms:
  three readings of the words (a presence pass first): the same RssAnon, wall x0.973-1.006 but cycles x1.010-1.020 at
  16 / 32 threads (`[TIMING] x_mb_lookup_build` 2.755 -> 3.539 s at AB n16, the scan 8.711 -> 8.469 s), hence the byte
  counts; the byte counts in a vector of their own, freed after each build: cycles x0.977-0.987, but TLB shootdowns 6
  k -> 121 k at AB n16 and DB n16 x1.021 in wall (each trim interrupted every thread), hence the bytes in the run
  list's capacity; without the keep's change, ~3,000 flushes per run (49 k shootdown interrupts at 16 threads) from
  the chains build that followed packed builds; one `Csr::run` / `Ring` for both forms (a test of the form per word
  and per call): dc-megablast x1.009 of dev at one thread with the switch on or off, hence a scan per form; the size
  rule kept for packed heads: c516 / db16 max RSS +2.4 / +1.9 MiB over dev (7/10 rounds), hence no rule.
- `LOSAT_X_PAIRPARDISC=1` (in the strict set since 2026-10-17 (perf/perfmb-default-bundle, evidence at the end of this entry): on unless it is `0`; it does not depend on PAIRPAR) blastn `-task dc-megablast`, one query-subject pair
  on a pool of 2+ threads (or on one thread: `LOSAT_X_BNDCONE`): the two-hit seed stage of a subject chunk with a discontiguous megablast
  table (one or two templates) and the diagonal hash (query block over 8000 letters, `scan_range`
  0, unmasked subject, a subject the thread searches by itself, no debug or trace). The pool threads scan strips of the subject
  with the discontiguous scan kernel (both templates' chains in NCBI's order, the query context of
  each pair), and the owner thread folds the pairs in scan order through a copy of the two-hit
  diagonal-hash code (window test, `s_TypeOfWord`, ungapped extension, every insert;
  `blastn/blast_engine/x_pair_seed_disc.rs`). The helpers scan at most two strips per thread
  ahead of the owner and reuse the owner's consumed buffers. The owner's fold is the slower side
  for dc-megablast (most pairs reach the diagonal hash), so the stage stops gaining at about 4
  threads. Anything else (contiguous tables, which PAIRPAR takes in one-hit mode) runs
  the original code. With `LOSAT_TIMING` a
  `[TIMING] x_pairpardisc_seed` line gives the strips, seeds and the owner's fold, scan and wait
  time. Test: `LOSAT/tests/blastn_pairpar_disc.rs`.
  Default ON since 2026-10-17 (perf/perfmb-default-bundle; workstation (2x Xeon Gold 6338, 16 physical cores of one
  node with their siblings idle, alternating rounds, every output the expected one; perf-mb
  `logs/defaultset_ws_report.md` unless named); a reduction is 1 - with / without the switch, the other switches of
  the candidate set on): dc-megablast of 1,775 / 3,492 / 31,168 D. melanogaster cDNAs against 3 Mb of chromosome 2L
  (`-subject`) at 16 threads -73.8 / -71.0 / -62.3 % wall, at 4 threads -62.4 / -61.6 %; megablast and blastn: it
  does not act. Max RSS +10-72 MiB, below NCBI `-db` in every cell (494 vs 700 MiB at 16 threads on the 31,168
  cDNAs). Its helpers only yield on wasm32 (`std::thread::sleep` panics under the browser's WASI shim).
- `LOSAT_X_PAIRPARDISCARR=1` (in the strict set since 2026-10-20 (perf/dc-seed-batch, evidence at the end of this
  entry): on unless it is `0`; it does not depend on `LOSAT_X_PAIRPARDISC`)
  blastn `-task dc-megablast`, one query-subject pair on a pool of 2+ threads (or on one thread: `LOSAT_X_BNDCONE`):
  the two-hit seed stage of a subject chunk whose query block has NCBI's diagonal array instead of the hash (`kQueryLenForHashTable`: at most 8000
  letters with both strands, so a query of at most about 4000 bases; `LOSAT_X_PAIRPARDISC` takes the hash only).
  Same scope otherwise (`scan_range` 0, unmasked subject, a subject the thread searches by itself, no debug or trace). The pool threads scan
  the strips of `LOSAT_X_PAIRPARDISC` (same kernel, same strips), and the owner thread folds the pairs in scan
  order through a statement-for-statement copy of the two-hit diagonal-array branch of `process_offset_pairs`
  (`s_BlastnDiagTableExtendInitialHit`: the diagonal masked by `diag_mask`, `last_hit` and `flag` in
  `hit_level_array`, the hit length as a byte in `hit_len_array`; `fold_two_hit_array` in
  `blastn/blast_engine/x_pair_seed_disc.rs`). The array of such a block has at most 8192 + window diagonals, and
  its chunks carry few pairs (31,168 cDNAs against 3 Mb: 10 chunks, 42,107 seeds in all, against 188.7 M in the
  27 hash chunks); the stage replaces the one-thread scan of the whole subject, which NCBI's `CBatchSizeMixer`
  brings once per small batch (a batch of one short cDNA after a large batch).
  Anything else runs the original code. With `LOSAT_TIMING` a `[TIMING] x_pairpardiscarr_seed` line
  gives the chunks, strips, seeds and the owner's fold, scan and wait time. Test:
  `LOSAT/tests/blastn_pairpar_disc_arr.rs`.
  Default ON since 2026-10-20 (perf/dc-seed-batch; workstation (2x Xeon Gold 6338, node 1 booked through its FIFO,
  8 / 16 / 32 physical cores with the siblings idle, `numactl --membind 1`, rotated rounds, every output the
  expected one; v4-D `results/v4d10_probe1`); a ratio is the median of the per-round with / without, 65a2283b
  against origin/dev 7e73ec9c, both built on that machine): dc-megablast of 31,168 D. melanogaster cDNAs against
  3 Mb of chromosome 2L (`-subject`) in the batch loop (`LOSAT_X_QBATCHPAR=0`) x0.950 / 0.918 / 0.900 at 8 / 16 / 32
  threads (its 10 array chunks took 0.024 s each on one thread; `[TIMING] scan_lookup` summed over the batches
  1.57 -> 1.35 s at 16 threads), by default (QBATCHPAR slots, removed on 2026-10-28) x0.974 / 0.996 / 1.006; 1,775 cDNAs (the batch loop
  by default) x0.919 / 0.879 at 8 / 16 threads, x0.996 at one thread (8 rounds; the stage needed a pool then: see
  `LOSAT_X_BNDCONE`); 31,168
  cDNAs at one thread x1.004 (2 rounds); megablast and blastn of 31,168 cDNAs at 16 threads x0.991 / 0.999 (it
  does not act there). Max RSS unchanged in every cell.
- `LOSAT_X_PAIRPARNA=1` (in the strict set since 2026-10-20 (perf/bn-short-strip, evidence at the end of this entry): on
  unless it is `0`; it does not depend on PAIRPAR) blastn with a lookup
  table that is not two-stage (the lookup width is the word size: `-task blastn-short`'s word 7, and word sizes 4 to 8
  when the table has enough entries), one subject of up to 5 Mb (a longer one is searched in chunks in parallel already)
  and a query block on a pool of 2+ threads: the one-hit seed stage of the chunk with the diagonal hash (query block over 8000
  letters, both strands counted, so about 4000 query letters; `window_size` and `scan_range` 0, unmasked subject, no
  debug or trace). The pool threads scan strips (512 x max(threads, 4) of them, at least 512 subject offsets each: 512
  offsets on 3 Mb from 12 threads on) with the scan kernel of the reference
  (`run.rs`, the branch after `if two_stage_lookup_ref.is_none()`), each lookup hit becoming a seed with the context of
  its query offset, and the pairs are folded in scan order through a copy of the reference's one-hit diagonal-hash code
  (`last_hit` test, `s_TypeOfWord`, the ungapped extension on the context's query slice, the cutoff, every insert;
  `blastn/blast_engine/x_pair_seed_na.rs`); with `LOSAT_X_PAIRPARFOLD` (default) the fold is split by diagonal-hash
  bucket over the threads, without it the owner folds. Every lookup hit is a seed here (about 100 per subject offset for
  30,000 24-mers against 3 Mb), so the helpers scan at most two strips per thread ahead of the slowest fold. Anything
  else (the diagonal array of a query block up to 8000 letters, such as a few primers; one thread; several subjects or
  subject chunks searched in parallel) runs the original code. Before it, blastn-short against one subject ran on one core at any thread count. With `LOSAT_TIMING` a
  `[TIMING] x_pairparna_seed` line gives the strips, seeds and the fold, scan and wait time. Tests:
  `LOSAT/tests/blastn_pairpar_na.rs`, `x_pair_seed_na.rs` `seed_stage_on_strips_equals_one_scan`.
  Default ON since 2026-10-20 (perf/bn-short-strip; workstation (2x Xeon Gold 6338, node 1 booked through its FIFO,
  16 physical cores with the siblings idle, `numactl --membind 1`, rotated rounds, every output the same bytes in every
  arm; v4-B `results/b13ab`); a ratio is the median of the per-round with / without, a718932f (the switch opt-in, same
  code paths) built on that machine): blastn-short of 30,000 24-mers (taken at even steps from the subject) against
  3 Mb of D. melanogaster chromosome 2L (`-subject`, `-evalue 0.001`, outfmt 6) x0.510 / 0.284 / 0.189 / 0.094 at 2 /
  4 / 8 / 16 threads (30.3 -> 2.97 s at 16 threads, where NCBI 2.17.0 `-db` takes 29.4-29.8 s on one core), x1.001 at
  one thread (the stage needs a pool); 30,000 queries of 18-30 nt x0.110 at 16 threads (34.3 -> 3.78 s); outfmt 0
  x0.140 at 16 threads. Max RSS 191-197 -> 238 / 270 / 275 / 304-311 MiB at 2 / 4 / 8 / 16 threads (the seeds of the
  strips in flight and one table per partition), below NCBI's 392 / 414 MiB at 1 / 16 threads. `-task blastn` (the
  megablast table, where the switch does not act) at 2 threads: same instructions (7.50 / 7.51 G) and cycles. At 8 and
  16 threads the threads spend 1.2-1.8x the user time of the serial scan, mostly waiting for the slowest partition
  (cycles up at equal instructions, sys 2-5 s).
- `LOSAT_X_TBXQFRAMES=1` (TBLASTX, speed; in the strict set since 2026-10-18, on by default, in the pull
  request that added it; evidence at the end of this entry) the SEG of a query batch's frames is made on the search
  pool, one task per frame, instead of one frame after the other before the lookup table is built
  (`tblastx/x_query_frames.rs`, called from `search_query_batch` for every query batch, `-subject` and `-db`). NCBI
  filters the contexts of a batch one after the other (`blast_filter.c:1280-1290`); each task reads and writes only
  its own frame (letters, intervals, unmasked copy) and only reads the `SegMasker`, and SEG's two caches (the memo
  of SEG's entropy and the histogram memo of the fast SEG path) are `thread_local!` tables of a pure function of
  the window's letter counts, so every frame gets the serial loop's intervals. It acts only with a pool of two or
  more threads and two or more frames in the batch; with one thread, a build without the `parallel` feature,
  WebAssembly without `wasm-threads`, the switch off, or `-seg no`, the serial loop runs. The translation of the
  frames stays serial (at most a few ms per batch). BLASTX needs no such switch: its query SEG already runs on the
  pool (at 2+ threads). Tests: the unit tests of `x_query_frames.rs` (pool == serial at 1, 2,
  3 and 7 threads for three SEG parameter sets, with batches of one query and none; a probe called once per frame
  shows that the frames are masked by at least two workers of a 3-thread pool, and by none with the switch off or
  one thread) and `LOSAT/tests/tblastx_query_frames.rs` (off / default equal at 1
  and 4 threads, with several query batches, other SEG parameters and `-seg no`; the reach rows show the pool
  path at 4 threads and not at 1). On by default since 2026-10-18
  (perf/tbx-query-seg-par), on this evidence: the shadow over the TBLASTX fixtures at 2 and 4 threads
  (`check_losat.py`, `ci_fast_regressions.py --programs tblastx --all-cases`) compared 2,910 frames (58,542 intervals) in 242 batches of 200 runs,
  all identical, and every report equals the one without the switch; i9 (8 alternating off/on pairs after a
  warm-up, load 5-9, `LOSAT_TIMING`, every report its frozen sha256), 8 threads, medians: AP027131 x AP027133 (code
  4) SEG 0.223 -> 0.044 s, outside the search 0.381 -> 0.200 s, wall 1.864 -> 1.722 s (on/off 0.92); MelaMJNV x
  PemoMJNVA 0.048 -> 0.009 s, 0.247 -> 0.210 s (0.86); AvCLPV x PsCLPV 0.122 -> 0.024 s, 1.130 -> 1.028 s (0.91);
  NZ_CP006932 against itself 0.208 -> 0.042 s, 4.28 -> 4.22 s (0.96); the first 3 kb of MelaMJNV x PemoMJNVA 0.012
  -> 0.011 s; at one thread both arms run the serial loop (SEG equal; wall ratios 0.96-1.03 within the spread of
  the P and E cores). Max RSS at 8 threads within 2 MiB (AP027131 x AP027133 134-137 MiB, MelaMJNV x PemoMJNVA
  46-48 MiB, 3 pairs).
- `LOSAT_X_QPREPPAR=1` (BLASTN, all tasks; in the strict set since 2026-10-17 (perf/perfmb-default-bundle, evidence at the end of this entry): on unless it is `0`) the BLASTNA encoding of a query batch and the ungapped Karlin blocks of its contexts are made on the search pool
  instead of one query / context after the other (`blastn/blast_engine/x_query_prep_par.rs`, called from
  `search_query_batch`, so for every query batch; `-db` searches too).
  NCBI encodes the queries of a batch into one buffer and computes each context's block in turn
  (`blast_setup_cxx.cpp:500-603`, `blast_stat.c:2754-2784`); both run before the lookup table is built, with no other
  thread working (about 0.03 s per 5 Mb batch; 5 % of a batch at 4 threads, 11-14 % at 8-16). With the switch on and a
  `SearchPool` of two or more threads current, both are ordered parallel maps run with `SearchPool::install`, so on the
  pool of the search that calls them (the whole search pool), never on rayon's global pool. Encoding: the buffer is cut with `split_at_mut` into the span of each
  query (first context start to last context end; the spans must be in query order and disjoint, else the serial function
  runs) and each query writes its own span with the two calls of the serial function, plus context then minus context.
  Blocks: `par_iter().map().collect()` over the contexts (rayon splits them into ranges; the blocks are per context, the
  jobs are not). Exact because the encoding of a query is a function of its letters, its two context offsets and the block
  length (constant tables, and the reverse-complement table, a `OnceLock` that every caller fills with the same 256
  bytes), and the block of a context is a function of its bytes and the reward and penalty (locals only; `compute_karlin_params_ungapped`
  has none); there is no cache, memo table, counter
  or random generator, results are collected in order, and no floating-point operation changes order (the sums are inside
  one context). It acts only with a pool of two or more threads and two or more queries (encoding) or contexts (blocks);
  with one thread, a build without the `parallel` feature, WebAssembly without `wasm-threads`, or
  the switch off, today's serial code runs. With `LOSAT_X_STATS` (any switch state) `main` prints
  `[X_STATS] QPREPPAR encode_s=<s> karlin_s=<s> fallbacks=<n>`, the summed wall times of the two steps over every batch and the
  fallbacks (counted when the switch is on: batches of two or more queries on a pool whose spans were not in query order
  and disjoint, which the serial function encoded; always 0 in the layout of `search_query_batch`).
  Since 2026-11-01 (line SF S3, family bn-query-read) the switch also makes the DUST masks of the batch on the pool
  (formerly `LOSAT_X_QDUSTPAR`, in the strict set from 2026-10-17, same rule: a pool of two or more threads and two or
  more queries; `blastn/x_query_dust_par.rs`, called from `coordination.rs` `apply_dust_masking`, so by `query_masks`
  for every batch and by `chunk_query_masks` for the parts of a split query; the per-query map is
  `par_iter().map().collect()` in query order with `SearchPool::install`, exact because a mask is a function of the
  query's letters and the three DUST parameters only; `-dust no` has no map at all). With `LOSAT_TIMING` `main` prints
  `[TIMING] query_dust: <s>s` and with `LOSAT_X_STATS` `[X_STATS] QDUSTPAR dust_s=<s> pool=<n>`, the summed wall times of
  the DUST calls and how many of them made their masks on the pool, only when at least one batch was masked. The numbers that justified the parallel masks (QDUSTPAR alone, 16
  threads, workstation): the DUST step of 1,775 cDNAs 0.211 -> 0.019 s, wall 2.19 -> 2.01 s (x0.92); blastn of them
  0.211 -> 0.022 s, x0.91; outputs identical.
  Measured on the workstation (2x Xeon Gold 6338, powersave, node 1 CPUs 32-47 with the siblings idle, `numactl
  --membind`; `-outfmt 6`, subject 3 Mb of D. melanogaster 2L; both arms run the same binary with `LOSAT_X_PAIRPARDISC`,
  `LOSAT_X_QBATCHPAR` (removed on 2026-10-28), `LOSAT_X_DIAGLV`, `LOSAT_X_PAIRPARFOLD`, `LOSAT_X_QDUSTPAR` and `LOSAT_X_DISCLUTPAR`, the off arm
  carries an unread `XOSAT_X_QPREPPAR=1` of the switch's length; alternating pairs, the order rotated per round, a
  discarded warm-up run before each pair; at most 1 other CPU busy at start and end of every run), dc-megablast, medians,
  on/off wall (per-round ratio median, range): (a) 1,775 cDNAs (5.0 Mb) 4 threads 0.653 -> 0.629 s (0.966, 3 rounds), 16
  threads 0.498 -> 0.475 (0.941, 4 rounds); (b) 3,492 cDNAs 4 threads 1.243 -> 1.208 (0.963), 16 threads 0.643 -> 0.596
  (0.922); Ad 31,168 cDNAs (37 batches) 4 threads 10.24 -> 9.73 (0.949), 16 threads 3.43 -> 3.34 (0.974). The two steps
  (`encode_s` + `karlin_s`, one stats run per arm): (a) 16 threads 0.034 -> 0.008 s, (b) 0.100 -> 0.020, Ad 0.55 -> 0.17
  (n16) and 0.53 -> 0.155 (n4). Cycles:u within 4 % (the cycles move with the GHz, 1.9-3.1); max RSS 99 / 99, 103 / 104, 109 /
  111, 194 / 195, 238 / 236, 507 / 510 MiB (off / on, in the order above). Every output equals the sha256 of the
  run without the switch. Tests: the unit tests of `x_query_prep_par.rs` (pool == serial at 1, 2, 3 and 7 threads on random
  queries with N, IUPAC letters, lower case, empty queries, batches of none / one / two, in the production layout of
  `search_query_batch` (which never takes the fallback; also with bytes 0..=255 on pools of 2, 3, 5 and 16 threads), in
  compact layouts with touching spans or sentinel gaps, and in layouts the parallel path refuses (the fallback); blocks bit
  for bit; a probe called once per item shows that both steps run on at least two workers of a 3-thread pool, and on none
  with the switch off, one thread or one item; a caller that is not a pool worker still runs on the given pool (fails
  without `SearchPool::install`); a pool other than the search's and nested calls) and
  `LOSAT/tests/blastn_query_prep_par.rs` (off / on equal at 1 and 4 threads for the four tasks, `-dust no`, other
  scores, `-lcase_masking`, split query chunks, the `[X_STATS]` line).
  Default ON since 2026-10-17 (perf/perfmb-default-bundle; workstation, perf-mb `logs/qpreppar_report.md`,
  alternating rounds): dc-megablast of 1,775 cDNAs against 3 Mb x0.966 / x0.941 at 4 / 16 threads, of 3,492 x0.963 /
  x0.922, of 31,168 x0.949 / x0.974; max RSS within 3 MiB.
- `LOSAT_X_QGROUPPAR=1` (in the strict set since 2026-10-19 (perf/bn-query-groups, evidence at the end of this entry): on unless it is `0`)
  BLASTN, every task, one subject at a time searched on one thread with two or more threads in the pool (one subject,
  or several one at a time under `LOSAT_X_BNFEWSUBJ`): the preliminary gapped
  loop of a subject chunk and the traceback of the subject run per query on the search pool instead of on the thread
  that owns the subject over all queries in one order (`blastn/blast_engine/x_qgroup_par.rs`, dispatched in
  `search_query_batch`). NCBI runs both loops of a subject on one thread with one interval tree
  (`blast_gapalign.c:3850,3917-3922`, `blast_traceback.c:375,404-405`); with one subject (`-subject` with one record,
  the loop's batches of a query file) the look-ahead helpers only computed alignments ahead of that loop, and the two
  loops with the purges and sorts around them were 57 % of a blastn batch's wall at 16 threads (5,000 and 31,168 cDNAs
  against 3 Mb of 2L; v4-D session 8 probe, `results/s8_probe_i9/REPORT.md`). Form: the chunk's ungapped hits are split
  by query in the order the seed stage saved them (a counting sort), each query task sorts its hits with the batch's
  comparator (stable), runs the same loop body (the loop is a closure, `x_prelim_list`, called once on the owner with
  the look-ahead helpers or once per query without them) and purges its list; the purged lists are concatenated in query
  order and sorted by score as before (rayon's stable `par_sort_by` on the search pool). The kept HSPs of the traceback
  are the queries' hit lists one after another, so each query's slice runs the traceback loop (`x_trace_list`, no
  helpers and no speculative batches in a task), the passes up to the second purge (`x_post_purge`) and the score sort,
  final tree and e-values (`x_final_list`); the final lists are concatenated in query order and sorted stably by score.
  The queries are handed out from an atomic counter, the query with the most HSPs first (ungapped hits per query before
  the preliminary loop, kept HSPs per query before the traceback), each thread with its own tree and alignment scratch.
  Why exact: the tree answers per context, the purges sort by context first and compare within a context (both the
  preliminary purge and the traceback's two: kept HSPs, then trimmed ones in reverse removal order,
  per context the same in both forms), the re-evaluation and identity test read one HSP and its context's cutoff and
  statistics, the sorts are stable and the serial purges' outputs are ordered by context first, so the concatenations in
  query order are the serial purges' outputs and the stable score sorts give the serial lists element for element
  (`chunk_prelim_hits` and the subject's final HSPs): the merge of a subject's chunks and `-subject_besthit` see the
  same lists. The owner loop stays (counted as `owner=`) where a stage has fewer than two queries with HSPs or one query
  holds more than a quarter of its HSPs (a long genome query: E. coli genomes against LC738868), for chunks searched in
  parallel (subjects above 5 Mb; their traceback still runs per query), with several subjects (the pool searches
  subjects), on one thread, and on the diagnostic paths (`-verbose`, `LOSAT_DEBUG_BLASTN`, `BLEMIR_DEBUG`,
  `BLEMIR_DEBUG_WINDOW`, `LOSAT_DEBUG_COORDS`, the blastn trace), which print in the loop's order. With the switch on,
  the look-ahead and the speculative traceback (also for megablast) act only where the owner loop
  stays. The chunk's ungapped hits are grouped by query in place (counting-sort destinations, 4
  bytes per hit, applied by swaps), so the stage holds no second copy. In the per-query form the `[TIMING]` steps inside
  the stages (`traceback_alignment`, the purge passes, `reevaluate_ambiguities`, the final tree) are sums over the
  threads, while `gapped_extend` and `traceback_prune` stay the stages' wall times. The one-thread path's instructions
  did not change (workstation, 1,775 cDNAs: blastn +0.003 %, megablast +0.001 %, dc-megablast -0.4 %).
  `LOSAT_X_STATS`: `[X_STATS] QGROUPPAR prelim= prelim_queries= traceback= traceback_queries=
  owner=` when a stage was planned (`owner=` counts the stages with HSPs that kept the owner loop); with
  `LOSAT_WASI_THREADS_DEBUG=1`, `[losat-thread-stage] program=blastn stage=prelim_queries|traceback_queries`. Tests: the
  unit tests of `x_qgroup_par.rs` (the rule and the plan's order, the in-place grouping, slices, `run` serially and on a
  4-thread pool, the pooled stable sort on ties) and `LOSAT/tests/blastn_qgroup_par.rs`
  (off / default equal at 1 and 4 threads for the four tasks, masking and other scores, split query
  chunks, the owner rule and the counters). Evidence (workstation ws, node 1, CPUs
  32-39 / 32-47 / 32-63 for 8 / 16 / 32 threads, governor powersave, THP madvise; dev 0797cc40 and this branch in one
  window, arms rotated per round, the run's CPUs warmed 3 s before every run; every report one sha per case, equal to
  NCBI's; v4-D `results/d10d_w1_summary.md`): in the batch loop (`LOSAT_X_QBATCHPAR=0`) blastn of 5,000 D. melanogaster
  cDNAs (15.4 MB) against 3 Mb of chromosome 2L took x0.79 / x0.62 / x0.51 of the loop's time at 8 / 16 / 32 threads
  (2.59 -> 2.06, 2.19 -> 1.35, 2.14 -> 1.09 s) and of 31,168 cDNAs (96 MB) x0.81 / x0.61 / x0.50 (14.92 -> 12.04, 12.60
  -> 7.70, 12.64 -> 6.32 s), with max RSS no higher (108-156 against 111-167 MiB) and fewer cycles (31,168 cDNAs at 16
  threads 234 against 344 Gcycles:u); megablast of 31,168 cDNAs and of three E. coli genomes against LC738868 and
  dc-megablast of 31,168 cDNAs x0.98-1.00 (their batches are bound by the lookup and the seed stage). By default
  (`LOSAT_X_QBATCHPAR`'s slots, removed on 2026-10-28, the stages per query only inside slots of two or more threads) the 15 cells took
  x0.95-1.03 of dev's time (medians of 4-8 rounds) at the same max RSS (G2's slots hold about 250 or 320 MiB from run to
  run in both arms), and 1,775 cDNAs (5.4 MB, a file the loop searches by default) x0.80 / x0.65 at 8 / 16 threads with
  blastn (0.93 -> 0.74, 0.80 -> 0.51 s), x1.00 with megablast, x0.97-0.99 with dc-megablast. With this switch the loop
  was still slower than the default slots of that date at 16-32 threads (31,168 cDNAs, blastn 7.70 against 5.89 s at 16 threads,
  megablast 1.23 against 0.96 s), at 121-250 instead of 287-1,484 MiB.
- `LOSAT_X_DIAGLV=1` (in the strict set since 2026-10-17 (perf/perfmb-default-bundle, evidence at the end of this entry): on unless it is `0`)
  blastn, every task, `-subject` and `-db`, with the diagonal hash (query block over 8000 letters): the
  table keeps per bucket an upper bound of the levels of its cells (the largest level an insert wrote since
  the table was made or cleared; 512 values), and a seed's retrieve of its own diagonal returns "no cell"
  without walking the bucket's chain when that bound lies below the seed's level bound: `s_off_pos + 1`
  and, with two hits, at most `s_end_pos - window_size` = `s_off_pos + word_length - window_size` (19 to
  24 below `s_off_pos` for dc-megablast's window 40 and templates 21 to 16), a constant per chunk added
  to `s_off_pos`. A cell below it makes the decisions of no cell (the `last_hit` test; the two-hit window
  test), which are the only readers of the retrieve; the insert still walks the chain itself, and the next
  seed reads only the table, so the bound needs no term for the insert, and the hit list and the table
  are the reference's (`blastn/blast_engine/x_diag_lv.rs`, exactness argument in its header). It acts in
  the reference fold (`process_offset_pairs`, not with debug, trace or window diagnostics, which print the
  values read) and in the strips folds of PAIRPAR and `LOSAT_X_PAIRPARDISC`; the off-diagonal
  retrieves (`scan_range`) stay walks. With the switch off the table has no bounds, an insert tests the
  length of an empty vector, and the retrieves walk as before (callgrind, release: +0.005 to +0.10 %
  instructions against the base).
  The unit tests of `x_diag_lv.rs` check the bound invariant (`the_bounds_hold_under_inserts_advances_and_clears`:
  every cell of a bucket lies at or below the bucket's bound) and its sufficiency
  (`retrieve_below_makes_the_decisions_of_retrieve`, `the_bound_is_tight_at_every_reachable_setting`: the cell the
  retrieve's walk returns gives the seed's decisions, the `last_hit` test and with two hits the window test, of no
  cell, where the walk is skipped); the tests of the strips folds compare the bounds of the two tables with the
  reference fold's (the one-hit strips fold had `LOSAT_X_PAIRPARSHADOW` until PAIRPAR graduated on 2026-10-22).
  Measured on the workstation (alternating pairs, node 1 CPUs 32-47, `LOSAT_X_PAIRPARDISC=1` in the dc-megablast
  runs; 1,775 or 3,492 D. melanogaster cDNAs against 3 Mb of chromosome 2L): dc-megablast one batch 2.237 ->
  1.806 s at 1 thread (0.808, 4 pairs) and 2.089 -> 1.304 s at 16 (0.624, 4 pairs), two batches 4.442 -> 2.663
  s at 16 (0.600, 2 pairs); blastn with the first query set 4.123 -> 4.136 s at 1 thread (1.003, 4 pairs;
  instructions +0.1 %: it walks about two cells per seed, so the skips save about what the bound per seed and
  its upkeep per insert cost) and 0.997 at 16 (2 pairs); peak memory unchanged. The shadow skipped 70 % of the
  retrieves on the dc-megablast inputs (62 % with the prototype's bound, which also had the insert's stale
  limit) and 69 % on blastn.
  Tests: the unit tests of `x_diag_lv.rs`, the strips tests of `x_pair_seed.rs` and `x_pair_seed_disc.rs`
  (with the bounds at 1, 2, 3 and 7 threads), and `LOSAT/tests/blastn_diag_lv.rs`.
  Default ON since 2026-10-17 (perf/perfmb-default-bundle; workstation (2x Xeon Gold 6338, 16 physical cores of one
  node with their siblings idle, alternating rounds, every output the expected one; perf-mb
  `logs/defaultset_ws_report.md` unless named); a reduction is 1 - with / without the switch, the other switches of
  the candidate set on): at 16 threads dc-megablast of 1,775 / 3,492 / 31,168 cDNAs against 3 Mb -8.6 / -8.1 / -12.0
  % wall, of 3,492 at 8 threads -26.5 %; plain blastn -2.0 % (its PR: x0.997 at 16 threads, x1.003 at one); max RSS
  unchanged.
- `LOSAT_X_DIAGLV2=1` (in the strict set since 2026-10-22 (perf/dc-fold-lean, evidence at the end of this entry): on unless it is `0`)
  blastn with a two-hit diagonal hash (window above 0: dc-megablast), `-subject` and `-db`, with `LOSAT_X_DIAGLV`
  on: the table also keeps per bucket the latest insert's `(diag, level, hit_len, hit_saved)` and `other`, an upper
  bound of the levels of the bucket's other cells (every insert raises `other` to the previous latest level, then
  becomes the latest). A retrieve that `LOSAT_X_DIAGLV`'s bound does not skip is answered without the walk when
  `other` lies below the seed's bound: the latest insert's values for its diagonal (the latest insert's cell is the
  first cell of its diagonal in the chain: an insert stops at the first cell of its diagonal or the first stale
  cell and writes there, or appends at the head, and no insert of the bucket came after it), no cell for any other
  diagonal (any cell it would find lies below the bound and makes the decisions of no cell, `LOSAT_X_DIAGLV`'s
  argument). The insert is unchanged apart from the memo, so the table and the hit list are the reference's
  (`blastn/blast_engine/x_diag_lv.rs`, (iv) to (vi) in its header). The table clear, `copy_hash`, the stage copy
  and the merge of `LOSAT_X_PAIRPARFOLD` keep the memo with the chains. Only the strips fold of `LOSAT_X_PAIRPARDISC`
  keeps and reads it: each stage makes it again from `LOSAT_X_DIAGLV`'s bounds at its start (`x_lv2_reset`: no latest
  insert, `other` the bucket's bound, which answers no cell exactly where the bound skips), so the reference fold
  (`process_offset_pairs`: one thread without `LOSAT_X_BNDCONE`, and every worker with many subjects) runs the code before the switch
  (`DiagHashTable::insert_plain`, `retrieve_below::<false>`) and may change the table between stages. (Keeping the memo
  in the reference fold made 100 16S queries against a 16S database at 16 threads x1.008 slower, where the walks it
  saves are one cell long.) One-hit tables (megablast, blastn, blastn-short) keep no memo, and the one-hit strips folds
  (PAIRPAR, `LOSAT_X_PAIRPARNA`) take the same memo-free calls. The dc strips stage runs an instance per state
  of the table's memo (`MEMO`, with `LOSAT_X_FOLDLEAN`'s `LEAN`: see there); with the switch off (or `LOSAT_X_DIAGLV`
  off) it is the instance without the memo.
  Tests: the unit tests of `x_diag_lv.rs` (`the_memo_makes_the_decisions_of_retrieve`: the bucket's walk
  finds the latest insert's values in the first cell of its diagonal and every other cell at or below
  `other`, and every answer makes the decisions of the walk), the dc strips stage test of `x_pair_seed_disc.rs` (the memo
  on 1 to 8 threads and 2 to 16 partitions), the partition tests of `x_pair_seed_part.rs`,
  `LOSAT/tests/blastn_diag_lv.rs` (`the_memo_and_the_lean_folds_give_the_same_report`: dc-megablast at 1 and 4
  threads with each switch, the memo answering; megablast keeps no memo).
  Evidence (ws, node 1 timed, dev f664ceee against perf/dc-fold-lean 33397e4a, medians of the per-pair ratios,
  pairs in rotated order, every output dev's; v4-A15 ab9 and the many-pair runs fin6; both switches on unless
  named): 31,168 dmel cDNAs
  against 3 Mb of dmel genome with `-subject`, 16 threads 1.842 -> 1.611 s (x0.872, 8 of 8 pairs; cycles:u 69.35
  -> 59.37 G, instructions 167.81 -> 150.18 G), 8 threads 3.083 -> 2.660 s (x0.862), 1 thread x0.998 (then the
  reference fold, no memo; since `LOSAT_X_BNDCONE` one thread runs the strips stage too: see there); the same queries against that genome as `-db`, 16 threads 2.065 -> 1.786 s (x0.864);
  1,775 cDNAs x0.890 (16 threads, 60 pairs) and x0.997 (1 thread); 236 cDNAs (a small table: the memo answers 7 %
  of the retrieves) x0.995 / x0.985 / x0.981 / x1.001 at 16 / 8 / 4 / 1 threads (300 / 150 / 80 / 60 pairs); 100
  16S queries against a 16S `-db` (every worker in the reference fold) x0.996 with dev's instructions; megablast
  and blastn (no memo) with dev's instructions: the 31,168 cDNAs with blastn at 16 threads x1.002 [0.999-1.004]
  (16 pairs, CPU time x0.999), megablast at 1 thread x0.998, 100 16S queries with megablast at 16 threads x0.999.
  Both switches off: x0.991 (31,168 cDNAs, 16 threads, instructions x0.9996) and x0.990 (236 cDNAs, 16 threads).
  Each switch alone (ab6, 85bf2fd2 against dev 249c7454, 31,168 cDNAs at 16 threads): this switch x0.927,
  `LOSAT_X_FOLDLEAN` x0.957. Before the partition states had their own lines, 236 cDNAs at 16 threads were x1.014
  (both off x1.025: the same instructions, memory-ordering machine clears x1.37). Memory (RssAnon at the peak, the
  same runs): `-subject` at 16 threads 417,916 -> 415,444 kB (`MALLOC_ARENA_MAX=1`: 382,040 -> 374,184 kB);
  `-db` at 16 threads with `MALLOC_ARENA_MAX=1` 118,408 -> 118,688 kB, within the runs' spread (with the default
  arenas it takes two values, by where the arenas lie).
- `LOSAT_X_FOLDLEAN=1` (in the strict set since 2026-10-22 (perf/dc-fold-lean, evidence at the end of this entry): on unless it is `0`)
  blastn dc-megablast, the strips folds of `LOSAT_X_PAIRPARDISC` (`fold_two_hit`, diagonal hash) and
  `LOSAT_X_PAIRPARDISCARR` (`fold_two_hit_array`, diagonal array): the stage runs one instance per state of the
  switch (`LEAN`) and of the table's memo (`MEMO`, `LOSAT_X_DIAGLV2`), `seed_stage_of::<LEAN, MEMO>` /
  `seed_stage_array_of::<LEAN>`, so with the switch off on a table without the memo it is the code before the two
  switches (the same instructions apart from the layout: the field offsets of the table, which gained the memo, and
  the partition states' indexing). The `LEAN`
  instance takes `s_TypeOfWord`'s early return
  itself (na_ungapped.c:531 `if (word_length == lut_word_length) return 1;`, always so for a discontiguous word: no call, no
  arguments, the offsets unchanged and nothing extended) and inlines the diagonal-hash insert
  (`DiagHashTable::insert_inl`: the statements of `insert_plain`, or of `insert_memo` with the memo). Same statements in the same order. The reference
  fold (`process_offset_pairs`) is not changed (one-thread searches reach the strips folds with `LOSAT_X_BNDCONE`).
  Tests: the dc strips stage tests of `x_pair_seed_disc.rs` (both instances on a table with bounds and memo),
  `LOSAT/tests/blastn_diag_lv.rs` (`the_memo_and_the_lean_folds_give_the_same_report`).
  Evidence: with `LOSAT_X_DIAGLV2`'s (above); this switch alone (ab6) x0.957 on 31,168 cDNAs at 16 threads
  (instructions 166.53 -> 152.44 G) and x0.983 on 236 cDNAs at 8 threads.
- `LOSAT_X_PAIRPARFOLD=1` (in the strict set since 2026-10-17 (perf/perfmb-default-bundle, evidence at the end of this entry): on unless it is `0`) blastn, where the strips paths of PAIRPAR (one-hit fold: megablast, blastn) and
  `LOSAT_X_PAIRPARDISC` (two-hit fold: dc-megablast) run, on a pool of 2+ threads: the diagonal-hash fold
  of the chunk is split by bucket (`bucket % P`, one partition per thread) over the owner and `threads - 1`
  jobs of an `in_place_scope` (`blastn/blast_engine/x_pair_seed_part.rs`). A chunk's scan is cut into 64 x
  max(threads, 4) strips (the one-hit fold since perf/bn-few-subj, before 8 x max(threads, 4): with a long chunk
  those strips of seeds held a quarter of its seeds; see `LOSAT_X_BNFEWSUBJ`). Each participant scans the
  next strip when it lies within `2 x threads` strips of the slowest partition, cuts its seeds by
  partition (a stable counting sort) and folds partition p of the next strip of any partition it can
  lock (each partition's state on cache lines of its own, 128-byte aligned: the participants poll every
  partition's lock while one of them folds into it, perf/dc-fold-lean; the same in the one-hit and NA strips stages), into p's own table, with the folds of the
  strips paths (`fold`, `fold_two_hit`, a const
  parameter and a sink of (strip, position) keys for every pushed hit and every insert that appends a
  cell). The merge puts the hits in key order and gives every appended cell the number `occ0 + rank` of
  its key (links and heads mapped; the capacity of the reference's doublings; the bounds of
  `LOSAT_X_DIAGLV` per bucket), and writes the changed buckets into the chunk's table. Exact because a
  seed's fold reads and writes only its own diagonal's bucket chain (retrieve, the `last_hit` and window
  tests, insert), `s_TypeOfWord` and the ungapped extension read only inputs, and with `scan_range == 0`
  there is no off-diagonal read; the argument is written out in the module's header. Every stage copies
  the chunk's table into every partition table; only their allocations stay in the subject scratch (until
  the query batch ends). A helper job returns at once on a participant's thread; a helper on any other
  thread leaves once the stage has stalled (no strip scanned, no partition folded) for 10 ms. Only that
  exit ends the wait of a helper that a thread runs while it waits inside a piece of a participant's
  nested parallel call that it stole (the pattern of #218). A stage also stalls for 10 ms when the OS
  deschedules the owner on a loaded machine or when one long fold unit holds back the window of strips;
  the idle helpers then leave for good and the stage ends on fewer threads (the same output, more time;
  `helpers_left=` counts them). No `thread::sleep` on wasm32. With LOSAT_TIMING the partitions'
  ungapped-extension times are added to the shared counters once per stage; with the switch the
  `[TIMING] x_pairpar_seed` / `x_pairpardisc_seed` lines end with `partitions= merges= merge=
  helpers_left=` (without it they are as before), their `fold=` is the CPU time summed over the
  participants (not wall time), and `[TIMING] scan_lookup` takes off only the owner's extension time, so
  it stays a figure of the owner's wall time. The one-hit
  fold's own shadow (`LOSAT_X_PAIRPARSHADOW`) went with PAIRPAR's graduation on 2026-10-22: the `=0` runs
  of `LOSAT/tests/blastn_pairpar_fold.rs` compare the reports, and the `[TIMING] x_pairpardisc_seed` fields
  (`merges=` equal to `chunks=`, `partitions=` the thread count) show that every chunk folded in partitions.
  RSS grows with the thread count: every partition table is a full copy of the chunk's table (P =
  threads, at most 512; about P x the table's capacity x 20 bytes), and up to `2 x threads` strips of
  seeds are in flight. On the inputs below the tables are small (thousands of cells per chunk), so the
  growth stayed below the spread of single runs, and max RSS with the switch did not rise from 16 to 64
  threads (one run each on the workstation's node 1, 64 threads on 60 CPUs): dc-megablast one batch 104 /
  98 / 96 MiB at 16 / 32 / 64 threads (off 93 / 90 / 89), two batches 119 / 121 / 111 (off 108 / 99 / 99),
  blastn 127 / 127 / 120 (off 142 / 158 / 133). With LOSAT_TIMING at 16 threads no helper left a stage
  of these inputs (`helpers_left=0`).
  Measured on the workstation (alternating pairs, node 0 CPUs 16-31, `LOSAT_X_PAIRPARDISC=1` and
  `LOSAT_X_DIAGLV=1` in every run, off runs with an unread placebo; 1,775 or 3,492 D. melanogaster cDNAs
  against 3 Mb of chromosome 2L): dc-megablast one batch 1.374 -> 0.797 s at 16 threads (0.580, 8 pairs,
  rounds 0.567-0.616), 0.687 at 8 (2 pairs), 0.835 at 4 (6 pairs); two batches 2.649 -> 1.527 s at 16
  (0.577, 4 pairs), with `LOSAT_X_QBATCHPAR=1` (slots, removed on 2026-10-28) 1.620 -> 0.959 s; blastn with the first query set 2.832 ->
  1.435 s at 16 (0.507, 4 pairs). Instructions +12 % (dc-megablast) and +4 % (blastn): the keys, the
  counting sort, the merge and the spinning of idle participants. Max RSS dc-megablast 90-92 -> 104-108 MiB
  (one batch) and 105-107 -> 119-123 MiB (two batches; the partition tables and up to `2 x threads`
  strips of seeds), blastn 137-151 -> 121-129 MiB. Against the prototype (`LOSAT_X_DCFOLDPART`,
  tmp/dcfold-instr d537d145) in the same windows: dc-megablast 0.99-1.01 of its time, blastn 1.14-1.29 (two
  windows, 4 rounds; instructions +0.8 %). Those numbers are of a build that still kept the partition tables
  from one stage to the next (no stage on the command line used them); the build that copies the table at
  every stage, against it with the switch on at 16 threads (4 alternating pairs each, node 1 CPUs 32-47, one
  speed level): dc-megablast one batch 0.720 -> 0.720 s (1.001), two batches 1.373 -> 1.370 s (0.997),
  blastn 1.357 -> 1.373 s (1.011, rounds 0.966-1.022; instructions equal), max RSS 3-6 MiB less for
  dc-megablast.
  Tests: the unit tests of `x_pair_seed_part.rs` (random op streams at 2, 3, 7 and 16 partitions on 1..8
  threads with and without bounds, growth past 256 cells, a stage after an offset move, the clear, a
  reference chunk in between, in every other case with the owner sleeping in its folds until a helper
  has folded so that helpers fold on any machine; a table changed between two stages; a nested parallel
  call in every participant; a panicking participant), the strips tests of `x_pair_seed.rs` and
  `x_pair_seed_disc.rs` (partitioned, two chunks in a row), and `LOSAT/tests/blastn_pairpar_fold.rs` (also
  the field names of the `[TIMING]` lines with and without the switch).
  Default ON since 2026-10-17 (perf/perfmb-default-bundle; workstation (2x Xeon Gold 6338, 16 physical cores of one
  node with their siblings idle, alternating rounds, every output the expected one; perf-mb
  `logs/defaultset_ws_report.md` unless named); a reduction is 1 - with / without the switch, the other switches of
  the candidate set on): dc-megablast of 1,775 / 3,492 cDNAs against 3 Mb at 16 threads -49.4 / -42.0 % wall, at 4
  threads -17.9 / -13.7 %; plain blastn of 1,775 cDNAs -11.6 % at 16 threads (its PR: x0.507 of the owner's fold). At
  4 and 8 threads plain blastn reads -1.0 % and -9.0 % (rounds at a matched clock -1.0 / +6.9 %), inside the spread
  of the same rounds (dc-megablast-only switches, which cannot act there, read -1.0..+7.1 %). Max RSS +7..+48 MiB in
  those rounds: dc-megablast of 31,168 cDNAs +39 MiB and plain blastn of 1,775 cDNAs +48 MiB at 16 threads (the latter
  beside the earlier budget's `LOSAT_X_QBATCHPAR` slots (removed on 2026-10-28), which the budget rule then kept off for that file), below
  NCBI `-db` in every cell (494 vs 700 MiB; plain blastn 122 vs 136 MiB under the budget rule).
- `LOSAT_X_PAIRPARHOME=1` (in the strict set since 2026-10-24 (perf/dc-fold-sync, evidence at the end of this entry): on
  unless it is `0`) blastn, in the partitioned stages of `LOSAT_X_PAIRPARFOLD`
  (above; the strips paths of PAIRPAR, `LOSAT_X_PAIRPARDISC` and `LOSAT_X_PAIRPARNA` with the diagonal hash): how the
  participants share the folds (`seed_stage_home` in `blastn/blast_engine/x_pair_seed_part.rs`). Partition p is the
  home of participant `p % threads`, which folds p's strips as they become ready, every ready strip from p's next one
  on in one lock (a run), and scans a strip within the window between runs. Another participant folds p only
  (`Run`): while p has no participant at home in the stage (its helper job has not started, returned at once on a
  participant's thread, or left a stalled stage), until that participant comes; one strip when p is the slowest
  partition and the participant has nothing else to do (the window is closed); and every ready strip when the
  participant has nothing else to do once every strip is handed out. Without the switch any participant folds the
  next strip of any partition it can lock (`try_lock` on every partition in every turn), one strip per lock: a stage
  of PAIRPAR or PAIRPARDISC has 64 x threads strips and so 64 x threads^2 fold units (65,536 at 32 threads; PAIRPARNA
  512 x max(threads, 4) strips), the partitions move from core to core unit by unit, and every unit takes the strip
  slot's lock and a shared pointer and adds to a per-strip and a stage-wide shared count. With the switch a strip's
  buffer goes back once the slowest partition has passed it (no per-strip count; the participant that sees the
  slowest partition move gives the passed strips back under the spare list's lock, and an idle participant only
  reads), the strips are read under a read lock, the helpers' stall test sums the cursor and the partitions' next
  strips (no stage-wide count), and the partitions' next strips and states lie on cache lines of their own. Each
  partition keeps its strips' order and one folder at a time (its lock), so the output is the same (the module
  comment's argument holds as it stands; only who folds which run changes); the owner alone can still scan every
  strip and fold every partition. With LOSAT_TIMING the `[TIMING] x_pairpar_seed` /
  `x_pairpardisc_seed` / `x_pairparna_seed` lines end with `runs= away= units=` (the runs folded, the strips of
  partitions folded away from home, and the stages' fold units, strips x partitions). Tests: the unit tests of
  `x_pair_seed_part.rs` run both ways (random op streams in alternating stages, with the strip buffers held within
  the window; the nested parallel call; the panicking participant; the table changed between stages; a stage that
  stalls past `HELPER_STALL`, so that the helpers leave), `home_partitions_are_folded_away_from_home` (a participant
  with slow scans, and more participants than the pool has threads), and `LOSAT/tests/blastn_pairpar_fold.rs` (the
  reports at 2, 4 and 8 threads also with `LOSAT_X_PAIRPARHOME=0`, and the field names of the `[TIMING]` lines with
  and without it).
  Default ON since 2026-10-24 (perf/dc-fold-sync; workstation (2x Xeon Gold 6338, node 1 booked through its FIFO,
  8 / 16 / 32 physical cores with the siblings idle, `numactl --membind 1`, rotated rounds, every output the
  expected one; v4-D `results/v4d11_w2`); a ratio is the median of the per-round with / without, 9381cc0a against
  origin/dev 7a2d300e, both built on that machine): dc-megablast of 31,168 D. melanogaster cDNAs against 3 Mb of
  chromosome 2L (`-subject`) in the batch loop (`LOSAT_X_QBATCHPAR=0`) x0.950 / 0.900 / 0.831 at 8 / 16 / 32 threads
  (4 / 6 / 6 rounds; cycles x0.950 / 0.878 / 0.740, instructions x0.97-0.98; the folds' CPU time summed over the
  participants 8.36 -> 6.73 s at 16 threads and 11.40 -> 8.27 s at 32, where 60,352 runs folded the 1,769,472 fold
  units of the 27 stages and 16,127 partition strips were folded away from home), by default (QBATCHPAR slots, removed on 2026-10-28)
  x0.958 / 0.963 / 0.975; 1,775 cDNAs (the batch loop by default) x0.950 / 0.953 / 0.891 at 8 / 16 / 32 threads; blastn
  of the 31,168 cDNAs at 16 threads in the batch loop x0.946, megablast x0.999 (rounds 0.977-1.009), megablast of E.
  coli genome pieces against LC738868 x1.023 (0.987-1.047; 0.22 s runs whose 7 stages carry no seeds and take 1 ms
  each with or without the switch; cycles x0.994). The first build of the switch (c7a29e25, every run of a
  non-home partition unbounded, window v4d11-w1) read the same in the same window (9381cc0a / c7a29e25 x0.994-1.001
  in the batch loop). Max RSS unchanged in every cell.
- `LOSAT_X_PAIRPARMERGE=1` (in the strict set since 2026-10-25 (perf/pairpar-merge, evidence at the end of this entry):
  on unless it is `0`) blastn, at the end of the partitioned stages of `LOSAT_X_PAIRPARFOLD` (above; the strips paths of
  PAIRPAR, `LOSAT_X_PAIRPARDISC` and `LOSAT_X_PAIRPARNA` with the diagonal hash): how `merge`
  (`blastn/blast_engine/x_pair_seed_part.rs`) puts the stage's hits in the order of their keys. Without the switch the
  owner sorts the (key, partition, index) of every hit of the stage on its thread while the other participants wait
  (`sort_unstable`; the AB input below: 76,600 hits per stage, 91 stages, 0.30-0.40 s of a 4.7-9.8 s run at 8-32
  threads, the only serial part of the stage, which otherwise scales in cycles). With the switch (`hits_by_strips`), every
  stage is ordered strip by strip: a partition's keys ascend (it folds the strips in order and a strip's seeds of the partition in their scan order, the
  counting sort of the strip being stable; one hit per seed at most) and a key is (strip, position in the strip's scan
  order), so the hits of strip i are one run per partition after those of the strips before i; each strip's runs are
  gathered and its few hits sorted (`hits_of_strips`). The strips are cut into `MERGE_RANGES_PER_THREAD` (2) ranges per
  thread of the pool, which the pool orders side by side (`rayon::scope`); a range starts each partition's cursor at its
  first key of the range (binary search) and writes its hits into its own slice of the list's spare capacity, after the
  hits of the ranges before it. A stage of one strip, or off a pool, is ordered on the owner's thread, the strips in turn
  (`hits_by_strips_here`). The order is the sort's (keys are unique), so the hits, and everything made from them, are the
  same. Until 2026-10-27 (perf/pairpar-merge-one-form) a stage of fewer than 4,096 hits kept the sort, a bound set above
  an estimated crossover and never measured; it is gone because no measured stage has 2,048 to 6,080 hits, so any bound
  there splits the stages the same way and the measured inputs give it no value (per-stage counts at the end of this
  entry), and the pool's form costs the small stages tens of microseconds each (`x_reach!` per
  stage, always the switch's side while it is on). Tests: `hits_by_strips_is_the_sort_of_merge` (random stages
  of 1 to 9 partitions and 1 to 300 strips, strips without hits, an empty partition, hits already in the list; pools of
  2, 3 and 8 threads, 1 to 4 x strips ranges, and the owner's thread with its `runs` left from the stage before) and `LOSAT/tests/blastn_pairpar_merge.rs` (the reports at 2 and 4 threads, and at 4 with
  `LOSAT_X_PAIRPARHOME=0`, against `=0`, the
  one-hit fold of `-task blastn` with 17,506 hits in its stage and with 412, the Na fold of `-task blastn-short`, the
  two-hit fold of `-task dc-megablast`, with the sides of `LOSAT_X_REACH`).
  Default ON since 2026-10-25 (perf/pairpar-merge; workstation (2x Xeon Gold 6338, node 1 booked through its FIFO, 8 /
  16 / 32 physical cores with the siblings idle, `numactl --membind 1`, the many-pair driver with the arm order rotated,
  every output the expected one: the 31,168-cDNA blastn NCBI's daa7e0ac, the dc-megablast NCBI's 5c357f81); a ratio is
  the median of the per-round with / `=0` of the same binary, 1d2884b1 against origin/dev e41ccb34, v4-A15
  `results/s12/ms1`): blastn of 31,168 D. melanogaster cDNAs against 3 Mb of chromosome 2L in the batch loop
  (`LOSAT_X_QBATCHPAR=0`) x0.960 / 0.949 / 0.947 at 8 / 16 / 32 threads (5 / 8 / 8 rounds, none slower; at 32 threads
  rounds 2-6 of every arm ran with 1.6 times the CPU time, a disturbance outside the run, and rounds 0, 1 and 7 read
  x0.941-0.948), against a `-db` of the chromosome x0.953 / 0.958 at 16 / 32 threads, with eight QBATCHPAR slots
  (`LOSAT_X_QBATCHPARSLOTS=8`; removed on 2026-10-28) x0.988 / 0.981, 1,775 cDNAs (the batch loop by default) x0.940 at 16 threads (20 rounds).
  Within noise where no stage reaches the bound or none is partitioned (16 threads): the 31,168 cDNAs by default
  (one-thread QBATCHPAR slots, removed on 2026-10-28) x1.003, 5,000 cDNAs x1.000, dc-megablast of the 31,168 cDNAs x0.999, of 1,775 x0.999,
  of 236 x1.001 (100 rounds), megablast x1.004 and in the batch loop x0.985, blastn-short of 236 cDNAs x0.996. CPU time
  x0.985-1.009, max RSS unchanged (`ru_maxrss`, every cell within 1 % or the arena spread of the arms). One thread (one
  node 0 core, `perf stat`, the 1,775 cDNAs with blastn and dc-megablast): instructions within 0.0001 % of `=0` and of
  dev, cycles x0.999 / x0.995 (no stage is partitioned at one thread).
  Probe (v4-A15 s12, workstation 2x Xeon Gold 6338, node 1 booked through its FIFO, 8 / 16 / 24 / 32 physical cores with
  the siblings idle, `numactl --membind 1`, dev f24d161f, every output NCBI's): blastn of 31,168 D. melanogaster cDNAs
  against 3 Mb of chromosome 2L in the batch loop (`LOSAT_X_QBATCHPAR=0`) 9.84 / 6.03 / 5.09 / 4.70 s at 8 / 16 / 24 / 32
  threads, the seed stage 4.50 / 2.47 / 1.92 / 1.68 s of it, of which the merge 0.40 / 0.36 / 0.29 / 0.30 s (its hit
  sort 0.354 of 0.360 s at 8 threads; the appended cells 1,300 per stage); the stage's other time scales with the
  threads in cycles (cycles:u of the run 201 -> 218 G from 8 to 32 threads; its CPU time grows with the clock, 3.06 ->
  2.40 GHz). The switch's form as a probe (`A15_P3=2`, exp/a15-seedscale 09aa43f4) in the same window, 6 rounds against
  the probe without it: x0.944 / 0.948 at 16 / 32 threads in the batch loop (x0.942 at 8), `-db` x0.956 / 0.955, eight
  QBATCHPAR slots (`LOSAT_X_QBATCHPARSLOTS=8`; removed on 2026-10-28) x0.991 / 0.985.
  One form without the bound (v4-A15 s17, 2026-10-27, perf/pairpar-merge-one-form): per-stage hit counts (a probe on
  origin/dev e15cc65e, exp/a15-mergeform 234ef1fd, `LOSAT_TIMING` per stage; workstation node 0, 16 CPUs shared with
  other jobs, so its times are indicative; 21 inputs at 16 and 4 threads, five forms each, every output the same across
  the forms, v4-A15 `results/s17/s17count_c1`). blastn of the 31,168 cDNAs (batch loop, slots, eight slots, `-db`): 91
  stages, 86 of 65,536-131,071 hits, 3 of 32,768-65,535, one of 6,081 and the last one of 358; of 5,000 cDNAs: 14 of
  65,536 or more, one of 22,883, the last of 358; of 1,775 cDNAs: 5 of 65,536 or more and the last of 358 (word size 9: 5
  of 65,536 or more); blastn-short of 236 cDNAs: 8,896, 97,259 and 509,202. dc-megablast of the 31,168 cDNAs
  (`-subject` and `-db`): 27 stages, none above 2,047; of 1,775 / 236 / 7 cDNAs: 0 and 1,086 / 0 and 53 / 0. megablast
  of the 31,168 cDNAs (both loops): 25 stages, none above 511; of 1,775: 0 and 198. No stage had 2,048 to 6,080 hits.
  Ordering time at 16 threads, the pool's strip order against the sort: 185 / 481 ms (31,168 cDNAs, slots), 62 / 315
  (batch loop), 42 / 74 (5,000), 4.5 / 44.5 (blastn-short); the same strip order on the owner's thread alone took 0.7-0.8
  of the sort, so the large stages need the pool. On the small stages the pool's start (`MERGE_RANGES_PER_THREAD` jobs
  per thread) cost 0.03-0.1 ms per stage more than the sort (dc-megablast of the 31,168 cDNAs 1.7 / 0.8 ms over its 27
  stages, `-db` 2.4 / 0.8, megablast in the batch loop 1.6 / 0.5 over 25), at most 0.1 % of the run, and about 5 us per
  stage at 4 threads, the same as the sort: the crossover moves with the threads and the pool's state, so a bound on the
  hit count alone has no measured value. Speed (workstation node 1 booked through its FIFO, 16 physical cores with the
  siblings idle, `numactl --membind 1`, the many-pair driver with the arm order rotated; arms origin/dev e15cc65e, a
  byte copy of it (A/A), b6b3a2d5 and b6b3a2d5 with `=0`; every output the same across the arms, dc-megablast of the
  31,168 cDNAs NCBI's 5c357f81, blastn in the batch loop NCBI's daa7e0ac; v4-A15 `results/s17/many_s17ab1`): the
  median of the per-round b6b3a2d5 / dev (A/A copy / dev) at 16 threads: dc-megablast of 236 cDNAs 1.011 (1.012, 60
  rounds), of 7 cDNAs 1.011 (0.986, 60; spread +-4 %), of 1,775 1.005 (1.001, 40), of the 31,168 1.001 (1.003, 16) and at
  8 threads 0.996 (0.997, 8), against the `-db` 1.003 (1.001, 12); megablast of 1,775 cDNAs 1.002 (1.003, 60), of the
  31,168 1.021 (1.003, 16; spread +-6 %), in the batch loop 0.998 (0.994, 16); blastn of 1,775 cDNAs 0.996 (1.001, 24),
  of 5,000 1.005 (1.004, 12), of the 31,168 in the batch loop 0.996 (0.998, 8); blastn-short of 236 cDNAs 1.003 (1.001,
  12): every case within the A/A spread. b6b3a2d5 / its `=0`: blastn of 1,775 cDNAs x0.943, of the 31,168 in the batch
  loop x0.947, of 5,000 1.002, blastn-short 0.995; the inputs with small stages only 0.976-1.035, inside their spread.
  CPU time x0.997-1.006, max RSS within 2 % (`ru_maxrss`).
- `LOSAT_X_QBATCHPAR` (BLASTN query-batch slots; REMOVED on 2026-10-28 with `LOSAT_X_QBATCHPARMEM`, `LOSAT_X_QBATCHPARSLOTS`,
  `LOSAT_X_QBATCHPARMISS`, `LOSAT_X_DCTWOSLOTS`, `LOSAT_X_BNEIGHTSLOTS`, `LOSAT_X_BNQSLOTS` and the `[X_STATS] QBATCHPAR` line;
  Owner decision: at 16 threads the slots' max RSS rise exceeded their time saving, 31,168 cDNAs x 3 Mb: blastn 107 -> 400 MiB for
  6.0 -> 4.6 s, dc-megablast 136 -> 238 MiB for 3.2 -> 1.8 s). A run ignores these variables. Every query batch is searched by
  the batch loop of `run_in_pool` on the whole search pool, as with `LOSAT_X_QBATCHPAR=0` before; outputs do not change. The
  design and the evidence are in the git history (the entries of the slot switches as of commit 70fbab3af).
- `LOSAT_X_RPTSTREAM=1` (BLASTN, all tasks, `-subject` and `-db`, memory and speed; in the strict set since 2026-10-17 (perf/perfmb-default-bundle, evidence at the end of this entry): on unless it is `0`) the batch loop of `run_in_pool` writes the
  reports of each query batch right after it has taken the batch (in batch order), finished its lists
  (`-db`) and given its gie to the mixer, and drops its final hit lists; after the loop it writes the
  epilogs (NCBI prints each batch's result sets in its batch loop and `PrintEpilog` after it,
  blastn_app.cpp:277-318). Without the switch the loop keeps the final lists of every batch and writes all
  reports once after the loop (`blastn/blast_engine/x_rpt_stream.rs`, exactness argument in its header).
  The same report code runs per batch (`write_reports`, with the batch's queries, lists numbered from 0, its
  per-query values and its part of the `-query_loc` placements); what is not per query is the outfmt 0
  prolog (before the search, as before), the outfmt 0 epilog (fixed inputs) and outfmt 7's last line, which
  counts the queries of every batch (NCBI's `m_QueriesFormatted`; `write_output_blastn_hitlists_to_writer`
  takes the count). NCBI's errors of a later batch (reader errors, an `eEOF` batch, a batch of queries
  without letters, `ErrorAfterInvalidBatches`, a header file that cannot be read) come after the reports of
  the batches before in both paths. LOSAT's own refusals come before every report without the switch, so
  with two or more batches the switch acts only where none can come after the first batch, and a subject
  set's tabular ids, outfmt 0 titles and `-lcase_masking` masks are then made once; with one batch (the
  first batch takes every query record) nothing is checked or made ahead (no refusal can follow its
  batch). Otherwise the reports stay after the loop, and `off=` names why: `custom_tabular` (a tabular
  field list or delimiter, or outfmt 10 or 20: `BlastnOutputFormat::CustomTabular`, whose refusals of
  taxonomy names and stored data come in its one report call), `outputs` (hit records, a format
  observer, a file sink, or two or more formats, which the library API allows on one stream and the batch
  writes would interleave: the Web adapter; the command line and ABI v1 are streamed), `debug` (`-verbose`,
  `BLEMIR_DEBUG`, `LOSAT_DEBUG_BLASTN`, `LOSAT_DEBUG_COORDS[_START]`, the blastn trace,
  `LOSAT_WASI_THREADS_DEBUG`: `run.rs` `x_debug_output`), `scoring`
  (no Karlin-Altschul values, or greedy gap costs LOSAT refuses), `chunk_size` (a negative `CHUNK_SIZE`),
  `search_space` (`(db + n (2 q + 1)) (2 q + 1)` above 2^63 - 1, q the longest query: the overflow refusal
  of the effective search space), `subject_ids` (a subject whose tabular id LOSAT refuses), and with
  `-db`, which the switch reads in one pass over the database's sequences before the search, in the pool
  (`database::x_streaming_refusal`): `titles` (a definition line that does not read as the reports read
  it, holds a byte outside ' ' ..= '~' or gives no tabular id, as blastp's block search checks) and
  `ambiguity` (a sequence whose ambiguity data the traceback cannot decode; the preliminary stage reads
  only the 2-bit letters, so the traceback of the batch whose lists name the sequence decodes them first;
  `SeqDb::check_ambiguity` reads only those bytes). With `LOSAT_TIMING` the `[TIMING]` blocks are printed
  when the loop takes a batch, as before, so with the switch a block comes after the query warnings of the
  batches before it (the same lines in another order). Memory: a batch's lists, flags, masks, blocks,
  search spaces and its queries' warnings are dropped after its write; kept for the whole search, as without
  the switch, are one (empty) warnings slot per query and with `-db` the stored title of every distinct
  sequence that a batch reported (`DatabaseTitles`: at most the database's sequences and at most the queries
  times the hit list size). With
  `LOSAT_X_STATS` and the switch, `[X_STATS] RPTSTREAM batches= queries= off= one_batch= check_s=`
  (`batches` written by the loop, `off` why the switch did not act or `none`, `one_batch` 1 where nothing was
  checked, `check_s` the seconds of the checks, the `-db` pass included).
  Measured on the workstation (2x Xeon Gold 6338, powersave, node 1 CPUs 32-47 with the siblings idle,
  `numactl --membind`; `-outfmt 6 -out FILE`; both arms with `LOSAT_X_PAIRPARDISC`, `QBATCHPAR` (removed on 2026-10-28),
  `DIAGLV`, `PAIRPARFOLD`, `QDUSTPAR` and `DISCLUTPAR`, one binary, the arms alternating with the order
  rotated per round and a warm-up run before each 16-thread cell; every output has its expected sha256),
  on/off as the median of the per-round ratios. D. melanogaster cDNAs against 3 Mb of chromosome 2L with
  `-subject`: 31,168 cDNAs with blastn (91 batches, 1,020,023 rows): 1 thread 71.3 -> 71.7 s (1.006, one
  run each), max RSS 502 -> 157 MiB (NCBI BLAST+ 2.17.0 `-db` 126 MiB, in another window); 4 threads 20.54
  -> 20.00 s (0.974, 2 rounds), 638 -> 295 MiB; 16 threads 7.37 -> 6.13 s (0.831, 3 rounds: the serial write
  of all reports after the loop now runs while the slots (removed on 2026-10-28) search), 1,076 -> 775 MiB. The same queries with
  dc-megablast (37 batches), 16 threads: 3.48 -> 3.47 s (0.992), 503 -> 489 MiB; 1,775 cDNAs with blastn (6
  batches), 16 threads: 0.568 -> 0.543 s (0.963), 329 / 330 MiB; with dc-megablast (3 batches), 16
  threads: 0.392 -> 0.406 s (1.033, 3 rounds) and in 8 more rounds of another window (1.8 GHz) 0.534 / 0.534
  s (0.998, 0.958..1.096), 103-104 MiB either way. What is left of the 1-thread memory is mostly the query
  file (about 94 MiB of records, held for the whole run). Shrinking a batch's lists to their length (the
  prototype `LOSAT_X_HELDSHRINK` of `tmp/split3`) adds nothing with the switch on: 157.8 / 157.5 MiB at 1
  thread, 777.6 / 782.3 MiB at 16 threads (medians of 3, spread 767-793). With `-db` (megablast; the slots
  (removed on 2026-10-28) did not run there, so the writes stayed serial): 100 16S rRNAs against all 16S (28,036 sequences,
  `BATCH_SIZE=15000`, 10 batches, 50,762 rows): the checks took 0.021-0.022 s at 16 threads (0.3 % of the
  6.4 s search) and 0.048 s at 1 thread (0.07 % of 65.7 s); wall 6.37 -> 6.46 s (1.013, 1.005..1.022, 3
  rounds), 65.6 -> 65.7 s at 1 thread; max RSS 178 -> 168 and 53 -> 33 MiB. 5,000 D. melanogaster cDNAs
  against the human cDNAs (465,769 sequences, 1.28 Gb, the largest nucleotide database on the workstation; 6
  batches, 2,201 rows), 16 threads: the checks 0.106-0.121 s (0.8 % of the 13.9 s search), wall 13.83 ->
  13.91 s (1.007, 0.990..1.031), max RSS 179 -> 173 MiB. Tests: the unit tests of `x_rpt_stream.rs` (a
  batch's part of the query input, every `off=` reason of `check` and the one-batch
  case) and of `database.rs` (`streaming_needs_every_title_and_every_ambiguity_record_to_read`: a missing
  header file and a refused title byte, and a raised ambiguity count in the old and the new format, whose
  error equals the traceback's), `LOSAT/tests/blastn_report_stream.rs` (off / on equal in exit
  status, report or `-out` file, and standard error, with `BATCH_SIZE` 4,000: the three tasks x outfmt 0,
  6, 7 x 1 and 4 threads; `-out` files;
  `-query_loc`, `-subject_loc -lcase_masking` with lower case subject runs, `-max_target_seqs -max_hsps
  -subject_besthit`; a reader error and a batch of queries without letters after earlier batches; `-db`,
  intact and with a refused title or a raised ambiguity count of a sequence that only the last batches hit
  (`titles`, `ambiguity`: the refusal, no row); and the `subject_ids`, `search_space`, `scoring`,
  `chunk_size` and `debug` cases) and `LOSAT/tests/blastn_report_stream_api.rs` (`run_local` with hit
  records and an observer, one and three formats: the formats, hits and observer ranges of one batch).
  Default ON since 2026-10-17 (perf/perfmb-default-bundle; workstation, perf-mb `logs/rptstream_report.md`,
  alternating rounds): blastn of 31,168 cDNAs against 3 Mb max RSS 502 -> 157 MiB at one thread (x1.006 wall), 638 ->
  295 MiB at 4 (x0.974), 1,076 -> 775 MiB at 16 (x0.831); dc-megablast of them x0.992 at 16 threads, blastn of 1,775
  cDNAs x0.963; dc-megablast of 1,775 at 16 threads x1.033 in three rounds and x0.998 in eight rounds of a second
  window (runs under 0.6 s).
- `LOSAT_X_QSTREAM=1` (BLASTN, all tasks, `-subject` and `-db`, memory; in the strict set since 2026-10-18
  (perf/bn-memory-rules, with a rule; evidence at the end of this entry, "Default ON"): on unless it is `0`;
  it acts only with `LOSAT_X_RPTSTREAM`, and not where its rule says `off=`)
  reads the query file twice instead of holding all its records for the whole search. Where the
  records were read until now (after `Query is Empty!`, before the search's checks and the outfmt 0
  prolog), a first reading with the same reader, by the loop of `read_queries`, keeps only each
  record's length and how the reading ended (its error and warnings): the batch planner and size
  mixer, `LOSAT_X_RPTSTREAM`'s checks and `-query_loc`'s query
  input (placements, skipped records, searched lengths, made from the lengths as `cut_queries`
  makes them) read only these, so every refusal, warning and error that comes before the first
  output comes from the same values at the same place. Then a second handle on the same open file
  (`try_clone`, at the position the first reading started from) reads the records forward as the
  batch loop takes its batches, as NCBI reads its query batches one at a time
  (blastn_app.cpp:277-279, blast_input.cpp:134-171). A batch's records are dropped after its reports
  (`LOSAT_X_RPTSTREAM`). Exactness (`blastn/blast_engine/x_query_stream.rs` header): the reader is a
  function of the bytes and its configuration; the second reading starts at the same byte with the
  same calls per record (`End()`, then `GetNextSequence`), so it reads the same records with the
  same local ids, titles, letters (lower case too) and warnings, and never reaches the error or
  `eEOF` that ended the first; each record's length is checked against the first reading's (a record
  whose length changed, a file that ends early or a reader error ends the search with LOSAT's error,
  and every later reading of the search returns it; a change that keeps every record's length is
  read as the new content, as NCBI, reading per batch, would read it); the batches are cut
  from the same lengths by the same code. Off (`off=`; the records are read before the search, as
  without the switch): `reports` (`LOSAT_X_RPTSTREAM` off: the reports after the loop read every
  record), `standard_input` (`-query -`), `not_a_file` (a named pipe, a device such as
  `/dev/stdin`: `metadata().is_file()`), `io` (no position or second handle, e.g. `try_clone` on
  WASI), and the rule (`x_query_stream::Rule::off`): `megablast` (`-task megablast`),
  `batch` (a file shorter than twice the largest query batch: `BATCH_SIZE`, or the `CBatchSizeMixer` maximum,
  999,000 letters for blastn and 4,999,000 for dc-megablast; the file's bytes against the batch's letters);
  `run_local` (the Web adapter, ABI v1) gets its
  records from its caller. Where
  `LOSAT_X_RPTSTREAM` is on but keeps the reports after the loop for a search (its `off=` reasons),
  the records of every batch are kept for them (the memory of the records read before the search).
  Since 2026-11-01 (line SF S3, family bn-query-read) the first reading is read in parts wherever the switch acts on a
  file of 4 MiB or more after the reader's start (`LOSAT_X_QSCANPARMIN`, tuning, sets the minimum), at two threads or
  more on a `unix` build with the `parallel` feature (formerly a switch of its own, `LOSAT_X_QSCANPAR`, in the strict
  set from 2026-10-26, `blastn/blast_engine/x_query_scan_par.rs`). The
  bytes from the reader's start to the end of the file are cut at targets at most 256 KiB apart (at least four parts per
  thread); a part runs from the first record start at or after its target (a `>` that begins a line, not followed by
  `?`) to the next part's, and a new reader with the same configuration reads it with the loop of `scan` on the calling
  thread and scoped helpers (the search pool does not exist yet; the parts read at offsets, `read_at`, so the search's
  reader stays where it was). A line that starts with `>` and not `>?` always starts a record, so the parts' lengths
  are the first reading's. The first reading falls back to one thread, which finds the error with its messages, on a
  `\r` or a `>?` anywhere in the file, a record start not within 256 KiB after a target, a byte of a part not read, a
  part whose reading ends otherwise than at the end of its bytes (an `eEOF`, a reader error or LOSAT's rejection of a
  Seq-id line), and for the BLAST archive,
  other targets than `unix` (WebAssembly, Windows), one thread and a file below the minimum. The minimum stays
  (2026-11-01): forced below it (`LOSAT_X_QSCANPARMIN=0`), blastn of 2.1 and 4.0 MB at two threads on the
  workstation's node 1 (an idle socket, cores warmed) ran 0.07-0.32 s longer in 11 of 32 runs, against 1 of 32 with
  the one-thread reading (medians x0.999-1.056), while at 16 threads the parts read x0.969-0.998, and on a loaded node
  0 the two-thread runs had the same user cycles and 0.04 % more user instructions in every arm; the cause is not
  known. Earlier evidence of the parts
  against one thread (workstation node 1, 16 threads, outputs identical and NCBI's): the first reading of 31,168 cDNAs
  (96 MB) 0.123 -> 0.009 s on the local i9, dc-megablast of them x0.932 wall at 16 threads (20 of 20 rounds) and x0.964
  at 8, blastn x0.972 (`-subject`) and x0.983 (`-db`), a file below 4 MiB within noise, max RSS unchanged. In the stats
  line below `scan_par=1` says that the parts' lengths were taken, `scan_par=0` the one-thread reading.
  With `LOSAT_X_STATS` and the switch, `[X_STATS] QSTREAM records= read= batches=
  kept= scan_s= scan_par= read_s= peak_window= off=` (kept: the searched records kept for a writer after their batch, 0 where
  `LOSAT_X_RPTSTREAM` writes each batch's reports; scan_s / read_s: the first
  and the second reading; peak_window: the most records read ahead of the loop at once). Measured on
  the workstation (2x Xeon Gold 6338, powersave, node 1, 16 cores with the siblings idle,
  `numactl --membind`; D. melanogaster cDNAs against 3 Mb of chromosome 2L with `-subject`,
  `-outfmt 6 -out FILE`; every arm with `LOSAT_X_PAIRPARDISC`, `QBATCHPAR` (removed on 2026-10-28), `DIAGLV`, `PAIRPARFOLD`,
  `QDUSTPAR`, `DISCLUTPAR` and `RPTSTREAM`, on vs off in the same window, one binary, the arms
  alternating with the order rotated per round and a warm-up before each 16-thread cell; every
  output has its expected sha256; on/off = median of the per-round ratios; window 4 of
  `perf-mb-20261014/logs/qstream_report.md`). Where it wins, max RSS: 31,168 cDNAs (96 MB) with
  blastn, 1 thread 158 -> 62 MiB (0.39; NCBI BLAST+ 2.17.0 `-db` 126 MiB in another window), 4
  threads 289 -> 209 (0.72), 16 threads 837 -> 755 (0.91); the same with dc-megablast, 1 thread 180
  -> 90 (0.50; NCBI 150), 16 threads 488 -> 413 (0.85); 13,193 cDNAs (40 MB), 1 thread 99 -> 61
  (0.61); wall on these 0.991-1.025 (one round each at 1 thread: 1.011, 1.025 and 1.014, and
  0.979 and 1.004 for the 96 MB file at 1 thread in another window; 0.991 at 4 threads, 2 rounds) and
  0.996-0.998 (16 threads, 3 rounds). Where it loses: at 16 threads on query files of 5-40 MB the slots (removed on 2026-10-28) read most of the file
  ahead anyway (`peak_window` 71-93 % of the records), so max RSS stays within -1.5..+0.4 %, and the
  first reading's extra parse comes before the search: 1,775 cDNAs (5.4 MB) 0.554 -> 0.565 s (1.018,
  8 rounds; 1.021 in another window; the same code without the switch's effect 1.001 / 0.995),
  5,000 cDNAs (15 MB) 1.026 and 40 MB 1.035 (noisy, flagged runs), 3,492 cDNAs (11 MB) 1.007. No
  rule found then that separates the two (8 threads and files of 40-96 MB were not measured), so the
  switch stayed opt-in (Owner's default-ON rule of 2026-10-17 with its correction of 23:0x).
  Default ON since 2026-10-18 (perf/bn-memory-rules), with the rule above; workstation, 16 physical cores of one
  node with their siblings idle, the dev default of 18b2d038 (and the candidate), a warm-up and 3-8 rounds per cell
  with the arms rotated, max RSS of both arms in the same window, cycles:u by `perf stat` (equal between the arms);
  switch / without, wall per round, and max RSS difference; v4-D `results/d10_w1_summary.md`, `d10_w2_summary.md`,
  `d10_w4_summary.md`, `d10_w5_summary.md`. Off by the rule: megablast, where the serial first reading (`scan_s`
  0.117 s for 96 MB) is several % of a fast search: 31,168 cDNAs 1.052 at 4 threads (-89.8 MiB) and 1.037 at 16
  (-116.6), 1,775 cDNAs 1.039 at 8 and 1.063 at 16 (0 / +0.9 MiB); its max RSS stays below NCBI's without the switch
  (31,168 cDNAs at 4 threads 241 MiB, NCBI 358.5); blastn at 16 threads where the slots (removed on 2026-10-28) read most of the file: 40 MB
  1.012 (-12.0 MiB), 15 MB 1.014 (-1.6). Acting: blastn of 31,168 cDNAs at 16 threads 0.998 (-77.1 MiB), 40 MB at 8
  threads 0.994 (-21.4) and at one thread 1.025 (3 rounds; -37.2), 15 MB at 4 threads 1.008 (-10.9), 1,775 cDNAs
  (5.4 MB) at 16, 4 and 1 threads 1.006, 1.000-1.006 and 1.013 (-4.7, -3.3..-3.5, -5.6); dc-megablast of 31,168 cDNAs
  at 16 threads 1.002 (-65.0). With the rule (a6258e6c; `d10_w5_summary.md`, node 1): dc-megablast of 1,775 cDNAs
  (5.4 MB, three batches of at most 5 Mb) at 16 threads says `off=batch` (the candidate without the batch part had
  acted there at 1.020 of the time for -1.8 MiB, `d10_w5o_summary.md`), megablast `off=megablast`; acting, blastn of
  1,775 cDNAs at 16 threads 1.012 (-6.4 MiB) and dc-megablast of 31,168 cDNAs 1.010 (-68.3), where the same window's
  placebo pair (the switch not acting in both arms) read 1.015; one thread on the 96 MB file, max RSS blastn 61.7 MiB
  (80.0 s) and dc-megablast 88.3 MiB (30.9 s), below NCBI BLAST+ 2.17.0 `-db` (126 / 150 MiB; 158 / 180 without the
  switch). Tests: the unit tests of `x_query_stream.rs` (the first reading equals
  `read_queries` on 40 seeded texts with LF / CR LF / mixed line ends, lines longer than the line
  reader's window, 1.1 MiB records, digits, comments and an error at the end; the second reading
  equals the first in batches of every size, with and without
  `-query_loc`; the query input from the lengths equals `cut_queries`'s; a file
  changed between the readings; `done` releases the loop's batch; the inputs that are read once) and
  `LOSAT/tests/blastn_query_stream.rs`
  (off / on equal in exit status, report or `-out` file and standard error, `BATCH_SIZE`
  4,000: the three tasks x outfmt 0, 6, 7 x 1 and 4 threads x `LOSAT_X_RPTSTREAM` off and on;
  `-out`; `-query_loc`, `-subject_loc
  -lcase_masking`, `-max_target_seqs -max_hsps -subject_besthit`; a reader error late in the file
  and in the first record, a batch of queries without letters, a batch of skipped records, LOSAT's
  refusals; standard input and `/dev/stdin`; `-db`; the switch off unless set; `kept=0` with
  `LOSAT_X_RPTSTREAM` acting).
- TBNPAR (graduated 2026-10-22; no variable, no `LOSAT_X_TBNPARSHADOW`; always on) tblastn:
  preliminary stage per (frame, chunk) and redo per match in parallel, one pool for the search. The
  unit-wise stage is the only preliminary stage of a search (at one thread, without a pool, it
  runs the units in order); the unit tests that check the stage against NCBI traces run this unit
  path through a `#[cfg(test)]` probe (`preliminary_hsps_traced`). The query batches go to `x_search_batches` at any thread count (since
  `TBNSSIDE` graduated); the one-batch-at-a-time loop searches a batch that `x_search_batches`
  did not take (a database whose sequences could not be read).
- TBNBATCH (graduated 2026-10-22; no variable; conditional: its rule stays) tblastn query batches
  searched side by side, appended in input order, with two or more threads, four or more query
  batches and at least one batch per thread (and the rule of `TBNBATCHRULE`). The batches share
  the kept subject translations (`TBNKEEPFEW` keeps them for such a search).
- TBNOIDCHUNK (in the strict set from 2026-10-13, graduated 2026-10-22; no variable; conditional:
  2+ threads and several subjects) tblastn: the subjects are cut into chunks and each chunk is
  one task that the pool gives to whichever thread is free (originally one entry into the search
  pool took four chunks per thread instead of one subject per thread; `TBNOIDNEXT` now hands
  the chunks out without entries);
  inside it the subjects of the chunk are a parallel range, so a thread that finds no chunk left
  takes part of one still being searched. A chunk is cut as NCBI's database reader cuts
  them for its search threads (`seqdbimpl.cpp:684-708` `x_FillSeqBuffer`, `seqdbatlas.hpp:496-501`
  `GetSliceSize`): consecutive subjects while their bases stay within
  `min(2^30, F) / (4 * threads) + 1`, at least one subject. NCBI's F is the largest file of the
  database; for `-subject` LOSAT takes ceil(bases of all subjects / 4), the size of the `.nsq`
  that `makeblastdb` writes for the same subjects (LOSAT's tblastn has no `-db` yet). So a subject
  longer than that is a chunk of its own and the subjects around it go to the other threads.
  Until 2026-10-14 a chunk was a hundredth of the subjects (NCBI's one-thread iterator,
  `blast_engine.c:1409`); on the fly genome (1,870 records, the seven chromosome arms first) that
  put the arms in one chunk on one thread (20 zebrafish proteins at 8 threads: 23.96 s against
  6.10 s without the switch). The subjects are searched by the same call and their lists reach
  the collector in subject order, so only the scheduling changes. It was off unless set (to any
  value) until it joined the set. Without it 28,036 subjects at 8 threads are 3,505 entries, each
  shorter than the wake-up of the workers; with it 13 (a hundredth of the subjects per chunk),
  17 (chunks bounded by residues, one per thread: about 217 subjects of 317 kb at most each) or 5
  (four such chunks per thread). An entry of one chunk per thread, the first form of the residue
  rule, lasted as long as its slowest chunk with nothing left for a thread that finished early:
  on the workstation (NUMA node 1) the zebrafish genome took 102.6 s at 8 threads against 77.8 s
  with the hundredth, and the transcriptomes and 16S 5-22 % longer, with the same CPU time. Timed with
  chunks of a hundredth on those subjects (16S rRNA genes, 40.6 Mb together) at 8 threads with other work on the
  machine: one query batch of 55 proteins 1.43-1.64 s without it and 0.76-0.85 s with it (NCBI
  BLAST+ 2.17.0 `-db` 0.68-0.78 s), one protein 0.87-1.11 s and 0.22-0.27 s; the same output and
  the same peak memory. With the chunks bounded by residues, four rounds on a loaded machine (load
  4-25; the first form, one chunk per thread): medians 1.35 s (55 proteins) and 0.61 s (one
  protein) against 1.25 s and 0.60 s with the hundredth and 2.81 s and 1.50 s without the switch;
  on the fly genome 7.0 s against 26.2 s with the hundredth and 6.8 s without the switch (NCBI
  `-db` 6.6-7.4 s; on the workstation 16.4 s against 62.1 s, NCBI 16.1 s). Nothing changes at one thread, nor where the batches are searched
  side by side (each batch on one thread). With `LOSAT_WASI_THREADS_DEBUG=1` the
  `[losat-thread-activity]` line of the subject stage ends with `pool_entries=1` (the search is
  one entry since `TBNOIDNEXT`); `LOSAT/tests/tblastn_oid_chunk.rs` reads it.
- TBNOIDNEXT (in the strict set from 2026-10-15, graduated 2026-10-22; no variable; conditional:
  2+ threads and several subjects) tblastn: no wait between entries. The threads of the pool
  take the next chunk when they are free, as NCBI's search threads ask the database for the next
  chunk (`seqdbimpl.cpp:333-344`), instead of waiting at the end of each entry for its last chunk.
  The thread that owns the search takes the lists back chunk by chunk in OID order and enters them
  into the collector as before, and searches a chunk itself while the next one in order is not
  finished (`tblastn/x_oid_flow.rs`). At most `4 * threads` chunks are taken or finished beyond the
  oldest one not taken back; a thread waits before a new chunk when that window is full, so the
  memory of the waiting lists stays bounded. The chunks and the search of each subject are those of
  `TBNOIDCHUNK`, so the report is the same. The entries of chunks that `=0` gave (an entry of four
  chunks per thread, or of one subject per thread) are gone. With
  `LOSAT_WASI_THREADS_DEBUG=1` the search is one entry (`pool_entries=1`) and the line adds
  `flow_chunks=` and `flow_peak_waiting=` (the most finished chunks that waited at once);
  `LOSAT/tests/tblastn_oid_chunk.rs` reads them.
- TBNOIDUNITS (in the strict set from 2026-10-15, graduated 2026-10-22; no variable; conditional:
  2+ threads and several subjects, a subject beyond a chunk's budget) tblastn: a subject longer
  than a chunk's budget (a chunk of its own, such as a chromosome) runs its (frame, chunk) units on
  the pool, as the single subject of a search does, instead of one after another on
  the thread that took it. Threads that have no chunk to search (at the end of an entry, or when no
  chunk is left) take its units; its merge in NCBI order, linking, E-values and reaping run on its
  thread after its last unit, so the report is the same. NCBI's thread searches a subject whole, so
  at 16 and 32 threads the longest chromosome bounds both programs; this goes beyond NCBI. Other
  subjects keep their units on their thread. The threads of the flow start
  together (each waits until all have entered the flow) so that a thread waiting for a unit can only
  take units, never a loop of the flow. With `LOSAT_WASI_THREADS_DEBUG=1` each
  subject whose units ran on the pool prints a `prelim_units` stage line with
  `parallel_selected=true`; `LOSAT/tests/tblastn_oid_chunk.rs` counts them.
- `LOSAT_X_TBXNEXT=1` (TBLASTX, speed; in the strict set since 2026-10-17, on by default; acts only
  where the parallel subject loop runs: two or more threads and more than one subject; native builds, and the
  threaded Wasm module (`wasm-threads`, LOSAT Web) since 2026-10-20) each thread of the pool takes the next subject from one
  counter when it is free (`tblastx/x_subject_next.rs`, one broadcast job per thread), with one worker
  state per thread made when it takes its first subject. Without it the loop gives rayon the subject indices
  as one indexed range (`for_each_init`), which rayon cuts into pieces before the search, so the long
  subjects a genome lists first (the chromosome arms) fall into one piece and one thread searches them in
  order. `process_subject` is unchanged and `report::final_hit_order` orders the HSP lists by query and
  subject index, so the reports are the same. It needs the windowed frames (TBXFRAMES, always on since 2026-10-26): without them it held the frames of
  several long subjects at once (3 fruit fly transcripts x the fruit fly genome database at 8 threads: 892
  MiB, NCBI 359 MiB). With both at 8 threads (node 0, CPUs 8-15, medians of 3 alternating pairs, reports
  identical): that search 11.48 -> 2.58 s (NCBI `-db` 2.73 s), max RSS 257 -> 316 MiB; with the genome FASTA
  as `-subject` 12.83 -> 3.14 s (410 MiB); 10 transcripts x the zebrafish cDNA database 8.51 -> 8.01 s.
  The threaded Wasm module (the LOSAT Web reactor under Node 26 on the i9-14900HX, both modules in one process
  alternating search by search, warm searches, every report identical): 3 query pieces of 1.5 kb against the
  Sakai and MG1655 genomes cut into 165 pieces of log-normal lengths (500 b to 600 kb) in length order (10.1
  Mb, the reference loop of each subject on one thread) at 8 threads 411 -> 160 ms (x0.391, 4/4 rounds faster)
  and at 4 threads 792 -> 327 ms (x0.569, 4/4), the CPU time divided by threads x wall 0.40-0.48 -> 1.02-1.04;
  against the 28 distinct genomes of `LOSAT/tests/fasta` in length order (25.2 Mb, the 3 E. coli first) at 4
  threads 2.41 -> 0.97 s (x0.444, 4/4); 25 genomes of 0.24-0.66 Mb x0.915 (8 threads, 4/4) and x0.960 (4
  threads, 8/8); the 28 genomes in reverse order at 4 threads x0.998 (4/8) and in length order at 8 threads,
  where PAIRPAR's units split each subject over the pool, x1.009 (3/8; 8 rounds at load 41-48, pair
  ratios 0.79-1.19).
- `LOSAT_X_TBXPREFETCH=1536` (TBLASTX, speed; a value switch, in the strict set since 2026-10-18
  from the pull request that added it, on by default; x86-64 only) the unit scan of
  PAIRPAR (`scan_unit` in `tblastx/x_pair_par.rs`) prefetches, for each hit of a
  `scansub` batch, the diagonal cell of the hit 16 entries ahead, for the subjects where one unit's
  diagonal table (4 bytes per cell: the query group's span plus the subject chunk) has at least the
  switch's number of KiB (2 to 4,194,304; built-in 1,536; `=0` off; `=2` every table of 2 KiB or
  more). A prefetch has no architectural effect, so the cells, the hits and the reports are those
  without it (`LOSAT/tests/tblastx_prefetch.rs` compares `=2` with `=0` and the default, with and
  without the probes); the other subjects run the scan without the hint, compiled as its own
  function item, whose machine code is the one before the switch (`objdump` comparison of the
  release binary with origin/dev, function by function). The `seed` time and `seed_ns_per_hit` of
  `LOSAT_TIMING`'s `x_unit` line include the hint where it acts. Measured (plan v4 P0-6; medians of
  8 rounds alternating, one warm-up, every report identical):
  - Workstation (Xeon 6338, L2 1.25 MiB per core; 8 booked cores of one node, other sessions'
    TBLASTX A/B runs on the node's other cores): Sakai x MG1655 (12.9 MiB tables) at 8 threads
    136.11 -> 116.97 s (x0.859), at one thread 776.4 -> 650.4 s (x0.838; 3 rounds, 750-881 /
    496-660 s under that load); the first 0.6 Mb of Sakai against the first 1.0 / 2.0 / 3.2 Mb of
    MG1655 (tables 2.1 / 3.4 / 4.8 MiB) at 8 threads x0.976 / 0.975 / 0.971 (2.37 -> 2.32, 4.36 ->
    4.26, 7.20 -> 6.99 s), at one thread x1.000 / 0.958 / 0.962 (12.32, 24.83 -> 23.79, 38.71 ->
    37.24 s). With the hint on every table: NZ_CP006932 against itself (1.67 MiB) x0.978 at one
    thread and x0.998 at 8, AP027131 x AP027133 (1.65 MiB) x0.987 at one thread and x0.990 at 8,
    AvCLPV x PsCLPV (1.18 MiB) at one thread x1.007: the built-in 1,536 KiB lies between the tables
    that gained and the one that lost.
  - i9-14900HX (L2 2 MiB per P-core; WSL2, threads not pinned, rounds started below load 4): with
    the hint on every table AP027131 x AP027133 x0.970 at one thread (5.85 -> 5.67 s) and x0.983 at
    8; Sakai's first 0.6 Mb x MG1655's first 2.0 Mb at 8 threads x0.981 (2.77 -> 2.71 s); Sakai x
    MG1655 at 8 threads (4 rounds, other sessions raised the load to 3-17 during the runs) 142.0 ->
    143.2 s, within the spread of the rounds (129-200 / 137-151 s; CPU time 968 -> 921 s).
  - Again after `LOSAT_X_TBXSCANBF` joined the default set (origin/dev 9e05b829 with this change,
    built-in 1,536 KiB, workstation, `=0` against the default): AP027131 x AP027133 (1.65 MiB, now
    above the built-in) x0.968 at one thread (11.30 -> 10.94 s) and x0.978 at 8; Sakai x MG1655 at 8
    threads x0.85 in the two rounds of 4 that ran without heavy load on the node (138.5 / 139.0 ->
    119.5 / 116.4 s).
- TBNSTATE (in the strict set from 2026-10-15; graduated 2026-10-22: no variable, always on, the
  fresh state per subject and the clearing per unit are gone) tblastn, the unit path: the
  scan state of a subject's (frame, chunk) units (the diagonal table of the two-hit word finder, the
  offset pairs, the initial hits, the gapped alignment memory) is kept by the thread from one
  subject to the next, as NCBI's search thread keeps its word-finder and gapped alignment structures
  for every subject of its loop (`blast_engine.c:1002-1003,1412-1413`), and each unit ends with
  NCBI's step of the diagonal offset (`Blast_ExtendWordExit`, `blast_extend.c:166-171`: the offset
  moves on by the unit's length plus the window; the table is cleared only when the offset reaches
  INT4_MAX / 4) instead of the next unit clearing the whole table. After the step every entry lies
  more than a window behind any hit of the next unit, so its first hit on a diagonal only records
  it, as on a cleared table: the same initial hits by construction. Before it each subject
  made its own state and each unit clears a table of next_pow2(query residues + window) entries of
  8 bytes (1 MB for 84,000 query residues); on the human transcriptome (465,769 subjects, 118 fly
  proteins, 8 threads) glibc's per-thread heaps were created and deleted about 250,000 times for
  these states (munmap 461,484, mmap 249,858 against 36 and 71 before #161). On the workstation (8
  threads on 8 cores of one NUMA node, other work on the machine, alternating rounds against a
  placebo variable, the same reports): that search 161.9 -> 133.2 s (CPU 1,261 -> 1,031 s, minor
  faults 3.8-22 M -> 0.87 M, munmap 415,942 -> 55), the zebrafish transcriptome 15.2 -> 12.0 s, the
  16S rRNA genes (28,036 subjects, 55 proteins) 1.18 -> 1.01 s, and at one thread 8.15 -> 6.82 s.
  The units of a long subject that run on the pool (`TBNOIDUNITS`, the single subject of a search)
  take the offset step too, on a state per task. The thread that runs the search keeps one set (the
  table and the largest offset-pair and hit buffers) until it ends or its next tblastn search takes
  it over. It acts at every thread count. Test: `tblastn/search_init.rs` (the table kept across
  units, across the clear and taken over by another query set, against a table made afresh for
  every unit).
- TBNXLATE (in the strict set from 2026-10-15; graduated 2026-10-22: no variable, no
  `LOSAT_X_TBNXLATESHADOW`, always on) tblastn: the six translated frames of a
  subject are read from tables, as NCBI makes them from the compressed subject
  (`blast_util.c:1071-1090`: one table of 64 residues per strand, `s_BlastGetTranslationTable`, and
  one table read per codon, `BLAST_TranslateCompressedSequence`). The subject's letters become
  ncbi2na in one pass; a gap or an ambiguous letter draws from the generator seeded with the
  subject's length, in subject order, as `s_Ncbi4naToNcbi2na` does
  (`blast_objmgr_tools.cpp:427-475`). The tables are filled from `GeneticCode::get`, so every
  residue is the one the letter-wise translation (`tblastx::translation::generate_frames`) gives;
  the unit tests in `x_subject_translate.rs` compare the two on every length 0-40, every genetic
  code and every IUPAC letter. The letter-wise translation of tblastn subject plans, the shadow and
  the `[X_STATS] TBNXLATE` / `[X_SHADOW] TBNXLATE` lines are gone, and so is
  `LOSAT/tests/tblastn_subject_translate.rs` (it compared the report with `=0` and the shadow).
- TBNSCANBF (graduated 2026-10-26, in the set since 2026-10-17, by the Owner's decision of that day on
  the timings below; no variable; always on) tblastn (speed; `tblastn/x_scan_bf.rs`): the
  lookup scan of the unit path (`search_seed.rs` `x_scan_unit`) takes the subject in
  blocks of 256 words. Each word's position and table index are written to a list whose length
  advances by the word's presence-vector bit, without a branch; the hits of the listed words are
  then copied in position order. When the offset array has no room for a word's hits, the scan stops
  at that word and the next call resumes there, as NCBI's `s_BlastAaScanSubject`
  (`blast_aascan.c:75-127`) does, so every call returns NCBI's pairs and scan range (test
  `x_scan_bf_equals_ncbi_scan_on_random_subjects` against `s_blast_aa_scan_subject_one_range` on
  random subjects, with offset arrays down to the longest chain). A backbone cell's hits are copied
  as three pairs, so the offset array has three entries more. Before, `LOSAT_X_TBNSCANBF=0` scanned with the reference
  loop; that call went from `x_scan_unit` (the loop stays in `tblastx/blast_aascan.rs` for the unit tests
  and the other programs). Tblastx and blastp keep their scans, and blastx scans
  this way always (BXSCANBF, graduated 2026-10-26). With `LOSAT_X_STATS=1` a line `[X_STATS] TBNSCANBF
  units= calls=` gives the scan units of the unit path this loop scanned and their scan calls
  (`tests/tblastn_scan_bf.rs`; the `_ref` counts of NCBI's loop went with the variable). perf stat on the workstation (one thread, a
  411-residue protein x the human genome `-db`, the same output): branch misses 1.96 G -> 0.46 G,
  cycles 136.0 G -> 96.7 G, instructions 285 G -> 348 G. Timed on the workstation (perf-si-20261012
  session 20, node 0, NCBI BLAST+ 2.17.0 `-db`, LOSAT without and with the switch in turn, one hash
  per case; time with the switch / without): human genome x one protein 0.73 at one thread (45.6 ->
  33.2 s, NCBI 44.9 s) and 0.77 at 8 threads (6.79 -> 5.21 s, NCBI 7.35 s); zebrafish genome x one
  protein at one thread 0.79 (31.3 -> 24.6 s, NCBI 31.3 s). With many query proteins the bit is set
  at most words and the branch is predicted: 0.98-0.99 (fly genome x one batch at 1 / 8 threads, E.
  coli MG1655 x 16 batches at 8, 16S x one batch at 8). For tblastn's table (word size 3, 5-bit
  letters) the candidate loop has its shifts as constants (`x_candidates::<3, 5>`; other tables take
  the general loop): a further 0.97 of the time on the one-protein cases (cycles 96.9 G -> 93.4 G),
  the same CPU time on the many-protein ones. With the table's mask (`0x7FFF`) it is a counted loop
  (`x_candidates_w3c5`, #273: cycles -6.0 % on the human and zebrafish genomes x one protein). An
  AVX2 form of the loop (`vpgatherdd` of the presence words) cut instructions by a third but was 2-7
  % slower (IPC 3.6 -> 2.3) and is not used.
- `LOSAT_X_TBNFLOWHELP=1` tblastn (speed; in the strict set from the pull request that added it, 2026-10-23;
  `=0` turns it off; a number from 2 to 1,000,000 sets the wait below in microseconds;
  acts where the subject flow of `TBNOIDNEXT` runs: 2+ threads and several subjects; `tblastn/x_oid_flow.rs`
  `wait`): a thread that waits in the flow after the start runs the pool's jobs meanwhile (`rayon::yield_now`)
  instead of sleeping on the flow's condition: the thread that owns the search when the next chunk in OID order
  is not finished and every chunk it could take is taken, and a thread of the pool when the window of
  `4 * threads` chunks is full. When the pool has no pending job the thread waits for a change of the flow (a
  chunk finished or taken back, the close, a failure; a change it missed while it looked at the pool ends the wait
  at once) for at most 100 us, doubled after each such wait up to 1.6 ms until it runs a job or takes a chunk, and
  looks again (where the OS timer is coarse, as on Windows, the waits are longer). The jobs it can run are the (frame, chunk) units of the subjects being searched
  (`TBNOIDUNITS`): the threads start together, so no loop of the flow is still queued in the pool once a thread
  waits, and a unit does not wait for the flow. A unit's result goes back to its subject by unit number and the
  merge stays on the subject's thread, so the report is the same. Before, on a genome database the owner that took
  a short subject first slept while the units of a chromosome another thread searched waited in the pool.
  With `LOSAT_X_STATS=1` each search by the flow (each query batch) prints `[X_STATS] TBNFLOWHELP help=<us>
  owner_idle_ms= owner_help_ms= owner_jobs= workers_idle_ms= workers_help_ms= workers_jobs=` (the time the owner and
  the other threads waited without and with a job, `help=0` with the switch off). Measured 2026-10-23 on the
  workstation (physical cores of node 1 with their HT siblings idle, `--membind=1`, cores warmed 3 s before each
  run, pairs with the order alternated; node 0 busy, load 22-77), the pull request's binary against its base, one
  protein against a genome (`-db -evalue 1e-5 -max_target_seqs 5`), every output the same bytes across arms: the
  fly genome at 8 threads 2.44 -> 2.30 s (x0.939, 8 pairs; the base's owner slept about 2 s in 4 of 8 runs, 2.52-2.60
  s, the branch in none, 2.27-2.33 s), at 16 threads 1.38 -> 1.33 s (x0.964), at 32 threads x0.984; the human genome
  at 16 threads 2.16 -> 2.09 s (x0.963), at 32 threads x0.990; the zebrafish genome at 8 / 16 / 32 threads x0.987 /
  0.998 / 1.002 (8 pairs, `perf stat` walls); 2,500 and 5,000 16S sequences with fly 2L between them (`-subject`) at 8
  / 16 threads x1.000 / 0.989; max RSS unchanged. With `LOSAT_X_STATS=1` on the fly genome at 16 threads the owner
  waited 1.17-1.20 s with `=0` and 0.06 s with the switch, during which it ran 2-3 parts of chromosomes (0.99-1.09
  s). Where the flow does not run (one thread, one subject such as an E. coli genome) the code is not reached:
  instructions the same.
- TBNFRAMEREUSE (graduated 2026-10-27, in the set since 2026-10-22 (#340), added 2026-10-19 with `TBNKEEPFEW`; no
  variable; always on; `tblastn/x_subject_translate.rs` `x_frame_buffer`, `x_reuse_frame`) tblastn (speed): a subject plan that is dropped gives its six frames to the next plan its thread
  translates from a packed subject, which writes every residue of them without clearing them first.
  A thread keeps at most six frames of at most 256 KiB (1.5 MiB). Before, every plan allocated and
  zero-filled six new frames: with `TBNKEEPFEW` a search of many subjects makes a plan per subject
  and query batch, and the allocations and zero fills were the 0.2-0.4 % that `TBNKEEPFEW` alone
  lost at one thread and side by side. The unit test `x_frames_of_a_dropped_plan_make_the_next_plan` checks that
  the next plan is made from the dropped plan's frames and is the reference plan. The residues do not change
  (`TBNXLATESHADOW` compared the
  plans with the reference translation: 266,092 plans of 5,000 fly cDNAs and of the zebrafish
  cDNAs at 8 threads, identical). With `TBNKEEPFEW`, against the default before (both `=0`), one
  query batch (no plan kept either way) of 27 fly proteins against the zebrafish cDNAs (`-db`)
  took x0.997 of the time at 1 thread (0.993-1.008, 6 pairs) and x1.000-1.004 at 8 (0.991-1.014,
  two runs of 6 pairs); the other measured cases are in `TBNKEEPFEW`.
- `LOSAT_X_TBNSCANBYTE=1` tblastn (speed; in the strict set from the pull request that added it,
  2026-10-21; `=0` turns it off; by default the AVX2 forms where the CPU has AVX2, the plain loops
  for sparse tables on x86-64 without AVX2 and in WebAssembly (from 2026-10-22), nothing on aarch64;
  `tblastn/x_scan_byte.rs`, exactness argument in its header): for tblastn's table (word size 3,
  5-bit letters) the lookup scan of the unit path (`x_scan_unit`) makes the table indexes of a block
  of 256 words from the letters directly (`l[s] << 10 | l[s+1] << 5 | l[s+2]`, 16 at a time with
  AVX2, no chain from word to word as in the rolled index) and lists the candidates from a presence
  vector with one byte per index (32 KB, built once per table and thread). The copy then takes each
  candidate once: while three hits for every candidate left fit in the offset array (one test, made
  again only after a cell of more hits), a cell of up to three hits is written without a test (one
  32-byte store with AVX2); a cell of more hits takes NCBI's test and is copied from the overflow
  array four pairs at a time; a candidate met while the block could fill the array takes NCBI's test,
  so the scan stops at the same word and every call returns NCBI's pairs and scan range (test
  `x_scan_byte_equals_ncbi_scan_on_random_subjects` against `s_blast_aa_scan_subject_one_range`,
  both forms). The offset array has four entries more. `LOSAT_X_TBNSCANBYTE=2` takes plain loops of
  the same form on any CPU for every table (one core on the workstation: one protein x the human
  genome x0.79, the LOSAT Web serial reactor one protein x E. coli x0.83, but x1.01-1.06 where most
  words are candidates); `=3` takes them on any CPU for sparse tables only, those where at most three
  tenths of the backbone cells have hits (the default without AVX2 on x86-64 and in WebAssembly). The
  first k proteins of AP027131 x MG1655 (one core, plain loops against the scan before): the LOSAT
  Web serial reactor x0.84-0.95 and x86-64 cycles x0.84-0.90 at a fill of 0.09-0.24 (k = 1-8), x1.00
  and x0.95 at 0.34 (k = 16), x1.02-1.03 and x0.99-1.00 at 0.42-0.46. Tables of other shapes, tables
  that are not sparse without AVX2, aarch64 and `LOSAT_X_TBNSCANBYTE=0` keep the scan before
  (`TBNSCANBF`).
  With `LOSAT_X_STATS=1` a line `[X_STATS] TBNSCANBYTE units= calls= avx2=` gives the units and scan
  calls taken this way (`tests/tblastn_scan_byte.rs` reads it). Timed on the workstation (one booked
  core, against `TBNSCANBF`, same output): one protein x the human genome `-db` x0.754 (20.75 ->
  15.65 s; cycles x0.720), fly proteins x the fly genome `-db` x0.988-0.990, 16 batches x E. coli
  `-db` x0.993-0.995, AP027131 x MG1655 x0.978-0.980, 55 proteins x 16S x0.984. At 8 threads (8 booked cores, warmed 3 s,
  4 rounds, medians of the per-round ratios): x0.828 (3.43 -> 2.84 s), x0.981, x0.997 (x0.977-1.021),
  x0.985, x0.982 in the same order.
  Probe: line F P3 (V2) and P3b (arm 7, `exp/f2-tbn-scan`).
- TBNHITBF (graduated 2026-10-26, in the set since the pull request that added it, 2026-10-19; no
  variable; always on) tblastn (speed) the two-hit loop of the unit path (`XTwoHit::pairs`, every offset array
  the scan fills; `tblastn/x_tbn_hitbf.rs`, exactness argument in its header) takes the two outcomes most
  hits end in, too far (`diff >= window`: `last_hit = s + diag_offset`) and overlap (`diff < wordsize`:
  nothing), without a branch: `diff = (s + diag_offset) - last_hit` (the reference's value, both
  wrapping in a release build), a candidate is `(diff - wordsize) as u32 < window - wordsize` (one unsigned
  compare), and the cell's next `last_hit` is selected (`s + diag_offset`, or the value it read). A
  flag-1 cell or a candidate stops the loop at that hit, which runs the body of `pairs` on its own, in
  order. Acts when `0 <= wordsize < window` (otherwise `pairs` runs). On the four ledger cases at one
  thread 83-87 % of the hits are too far, 12-15 % overlap, 0.5-0.7 % meet a flag-1 cell and 1.2-1.6 %
  are candidates, HsG_p1 included (v4-A10 probe 1). Measured (workstation, node 0: one physical core
  for 1 thread, 8 physical cores for 8 threads, HT siblings idle, the booked cores warmed 3 s before every
  run, arm order reversed every other round, one warm-up round + 4 rounds, medians of the per-round
  ratios against the dev binary before the switch; every report identical, max RSS equal; a gate ran on
  other node-0 cores meanwhile): AP027131 proteins x MG1655 (`-subject`) 1 thread x0.831 (5.47 -> 4.55 s)
  and 8 threads x0.852; Bac_d16 x0.832 and x0.817; FlyG_b1 x0.803 (20.93 -> 16.88 s) and x0.837; HsG_p1
  x0.890 (23.32 -> 20.74 s) and x0.909. Cycles x0.80-0.88, instructions x0.96-0.99, mispredicted branches
  x0.32-0.66. With `=0` (before it graduated) the reports were the same and the instruction counts equal
  those of the binary before the switch (x1.000), time x0.98-1.01. The 2026-10-17 conditional-move form of these tests
  (perf-si-20261012 decision 180, not kept) added 7 % instructions on HsG_p1 and removed 2.5 % of the
  mispredicted branches; this one keeps the instruction count and removes 68 % there.
- `LOSAT_X_TBNUNGAP=1` tblastn (speed; in the strict set from the pull request that added it, 2026-10-24;
  `=0` turns it off) the two-hit ungapped extension of the word finder (`XTwoHit::pairs`, every candidate of the
  unit path) takes BPUNGAP's walks (`blastp/x_ungapped_fast.rs`: one signed byte per residue pair from
  the 32 x 32 BLOSUM62 table, the bounds from the slice lengths as the reference computes them, so a chunk that
  ends inside a frame stops the walks where it stops the reference's) instead of the matrix-generic
  `extend_two_hit`; an extension with an X-drop outside 1..=127 runs the reference (`tblastn/search_init.rs`,
  `x_tbn_ungap`). Conditional: BLOSUM62 (tblastn's default) takes the walks, another `-matrix` (BLOSUM45) keeps
  the reference (`x_reach!` where the word finder is made). Found by the gate-gap audit (2026-10-19,
  `blastx_shared.md`: tblastn had no counterpart of BPUNGAP); ungapped extension was 13 % of the cycles of a
  short `-subject` search (50 AP027133 proteins x MG1655, local perf). Measured on the workstation (node 0 CPUs
  20-23, HT siblings idle, cores warmed 3 s, one warm-up round + 3 rounds, arm order rotating, medians of the
  per-round ratios against the dev binary before the switch, a806eaa5; reports and max RSS identical; other
  lines' gates on node 0 meanwhile, load 31-54): AP027131 proteins x MG1655 (`-subject`) 1 thread cycles x0.968
  (4.56 -> 4.43 s), 4 threads x0.971; HsG_p1 (`-db`) 1 thread x0.974 (14.30 -> 13.90 s); FlyG_b1 (`-db`)
  1 thread x0.964 (17.05 -> 16.58 s), 4 threads x0.964. Instructions x0.952-0.977. With `=0`: cycles
  x0.994-0.999 and instructions x1.004-1.009 of the dev binary.
  WebAssembly (serial WASI command, one Node process per case, CPU 20, 3 rounds, against `=0`): x0.981
  (AP027131 x MG1655, 5.56 -> 5.42 s) and x0.996 (FlyG_b1, 20.46 -> 20.39 s), linear memory equal.
- `LOSAT_X_TBNREDOLPT=1` tblastn (speed; in the strict set from the pull request that added it, 2026-10-25;
  `=0` turns it off; `tblastn/stage_d_pipeline.rs` `x_redo_lpt`, `x_redo_lpt_order`): the composition-based redo
  (`-comp_based_stats 2`, the default) of a chunk of matches on the search pool runs in one `rayon::scope` of
  `min(threads, matches to redo)` tasks instead of a `par_iter` over the chunk: each task holds one of the
  per-thread redo states and takes the chunk's matches one at a time from one counter, those with the most HSPs in
  their preliminary list first (equal counts in list order), and puts each result in the match's slot; the
  results are then applied in list order as before (early termination tested again, heap insertion). The
  `par_iter` had no `with_max_len(1)`, so rayon cut each chunk (32 matches per thread) into runs of adjacent
  matches, and the list is in query order: on the human genome the 24 chromosome matches of one protein
  (30,000-95,000 HSPs each) sat in one run and one thread redid them one after another while the other
  threads waited (B11f probe, HsG_b4 at 8 threads: 618 s of the redo's 1,087 s were waits). A match's result is
  a function of its own preliminary list and a redo state carries no value from one match to the next, so the
  reports are the same (unit tests `the_redo_lpt_order_takes_the_wanted_matches_largest_first`,
  `the_redo_lpt_scope_returns_the_par_iter_results_by_index`, `a_panicking_redo_match_propagates_and_the_pool_stays_usable`;
  `tests/tblastn_redo_lpt.rs`). NCBI runs the same matches in one `omp for schedule(static)` (`blast_kappa.c`
  3448-3451), whose fixed shares have the same problem. Conditional: two or more threads and two or more matches,
  on each chunk with two or more matches to redo (`x_reach!`); one thread (also each query batch searched side by
  side with others, on one thread) keeps the loop on the calling thread. Until 2026-10-27 a pool of more than 16
  threads kept the `par_iter` too (the last part of this entry). Measured 2026-10-25 on the workstation (2x Xeon
  Gold 6338, cores warmed 3 s before each run, pairs alternating in order, `perf stat` walls; the first commit of
  the pull request that added it against its base 313fd4e7; every output of a case the same bytes in both arms),
  proteins against the human genome (`-db hsap_dna -evalue 1e-5 -max_target_seqs 5`): the 122 proteins of four
  query batches (HsG_b4) at 8 threads 1,879.9 -> 1,218.7 s (x0.648, one pair, node 0 cores 24-31 with their
  siblings idle, load 9-50, output sha 7f761059 as before); the first batch alone (35 proteins, 396 matches) at 8
  threads 586.6 -> 349.3 s (x0.630, 3 pairs, same cores), at 16 threads 397.2 -> 284.0 s (x0.715, 3 pairs, node 1).
  At 32 threads (node 1) the batch took 271.5 -> 280.2 s (x1.035; pairs x1.038, x1.035, x0.985), hence a limit of
  16 threads (`X_REDO_LPT_MAX_THREADS`): a probe of the items (one line per match) showed the chunk's time is its
  longest match's (66,169 HSPs, 97-121 s) in both arms at 16 and 32 threads; at 16 the `par_iter` gave one thread
  a run of adjacent heavy matches (redo 226 s against 112 s), at 32 its runs are short and it started the longest
  match within 1-8 s, while running the 32 heaviest matches at once raised the redo's work by 5 %. Matches with few
  HSPs (node 1, 8 pairs at 8 / 16 threads): 16 proteins x E. coli x1.007 / x1.000, fly proteins x 5,000 16S x0.998
  / x1.014, the same against 2,500 + 5,000 16S around fly chromosome arm 2L (`-subject`) x0.981 / x0.980, against
  the fly genome x0.987 / x0.994, one protein against the zebrafish genome x1.002 / x1.001. Instructions equal
  those of the base in every case and arm (x1.0000), `=0` included, and at one thread (the same path). Max RSS
  within 2 % of the base (HsG_b4 2,813 -> 2,775 MiB).
  The limit came off in a later pull request (2026-10-27, `perf/tbn-redo-lpt-all`). While the redo still ran loops
  quadratic in a match's HSPs, the heaviest matches set its time at 32 threads in either order: with
  `REDOCONTAINIDX` (the 35 proteins at 32 threads 173.6 s) the switch allowed up to 64 threads still ran x0.991
  (3 pairs). After `TBNREAPIDX` (35.8 s) the B11e probe of v4-B10 session 14 (one line per match and per
  thread, node 1, one run each) found the 32-thread tail in the `par_iter` again: its redo took 14.49 s against
  5.13 s of work per thread (299.5 thread-s idle), one thread redoing six chromosome matches of query context 0 one
  after another, while largest first took 8.05 s, the longest match (context 31 x chromosome 1, 66,169 HSPs), and
  the run 35.86 -> 29.63 s (x0.826) for 1.7 % more user time (763 -> 776 s). Measured for that pull request on the
  workstation (as above; base 3a83c151, its commit 718c2b1f; node 1 whole, 3 pairs unless said; every output of a
  case the same bytes in every arm, the 35 proteins' 52de36f0 as on dev): the 35 proteins at 24 threads (CPUs
  32-55) 53.3 -> 33.2 s (x0.622, 3/3; user 697 -> 705 s), at 32 threads 36.1 -> 29.7 s (x0.823, 3/3; user 773 ->
  780 s; `=0` 36.2 s, one run), at 16 threads 41.6 -> 41.7 s (one pair, the same path, instructions x1.0000), and
  at 64 threads (the 32 cores with their siblings) 29.1 -> 29.9 s (x1.029; pairs x1.031, x0.995, x1.029): with
  every logical CPU busy, the heaviest matches started together share cores with one another. Max RSS 2,707 ->
  3,007 MiB at 24 threads (x1.110), 3,258 -> 3,449 MiB at 32 (x1.070), the same at 16 and 64: the heaviest
  matches are in flight together. Matches with few HSPs at 32 threads (8 pairs): 16 proteins x E. coli x1.003,
  fly proteins x 5,000 16S x1.028 (no hits, 0.09 s), the same around fly 2L (`-subject`) x0.996, against the fly
  genome x0.984, one protein against the zebrafish genome x1.010; instructions x1.0000 (x0.9999). At one thread
  (node 0, one core) 16 proteins x E. coli 81.411 G and fly proteins x the fly genome 127.161 G instructions in
  both arms.
- `LOSAT_X_REDOCONTAINIDX=1` blastp, blastx and tblastn (speed; in the strict set from the pull request that added
  it, 2026-10-25; `=0` turns it off; `core/composition_adjustment/x_align_index.rs` `x_contains`, `x_add_distinct`,
  `x_finish`, used by `redo_alignment.rs` `blast_redo_one_match_with_workspace_queries_and_matrix_observed`): the
  composition-based redo of one match keeps, for each query, the alignments it has redone so far (NCBI's list
  `alignments[query_index]` of `Blast_RedoOneMatch`, newest first) and asks two questions of it for every HSP:
  whether the HSP is contained in an alignment that scores enough more (`s_IsContained`, `redo_alignment.c`
  965-996) and, for every new alignment, whether one with a shared end scores at least as well, else which ones it
  replaces (`s_WithDistinctEnds`, 385-450). NCBI walks the whole list each time, so a match with tens of thousands
  of HSPs costs the square of that: on the human genome the B11g probe (35 proteins of HsG_b4, 8 threads, 396
  matches, 1.35 M HSPs) counted 63 G list nodes visited, 15 ns each (966 thread-s), 99.5 % of them in the 44
  matches with lists of 10,000 or more. With the switch a query's list moves, once it holds 32 alignments
  (`X_INDEX_FROM`), into an `XAlignIndex` of its own: the alignments in slots in the order they were kept plus a
  `BTreeSet` of (frame, match start, slot), and the two questions look only at the alignments whose match start
  lies where an alignment that answers them can start (bounded by the largest and smallest `match_end -
  match_start` added so far; the frames of the same sign, or the same frame, as NCBI tests; `s_IsContained`'s
  candidates nearest start first): the probe's shadow index visited 133,594 candidates in 1.9 thread-s and agreed
  on every call and every final list. The answers do not depend on the order the alignments are met in, the
  removals are the same alignments, and at the end of the match the slots that hold an alignment are linked newest
  first, which is the order of NCBI's list (the new alignment goes to the head, the others keep their order), so
  every value and the report stay the same (unit tests `the_index_answers_as_ncbis_list` (600 random runs of 1,200
  calls on three queries against `is_contained` / `with_distinct_ends`), `coordinates_at_the_ends_of_i32_give_the_same_answers`,
  `a_list_moves_into_the_index_once_it_reaches_the_threshold`; `tests/redo_contain_idx.rs`: one input per program,
  the switch on and off at 1 and 2 threads). Shorter lists stay NCBI's lists, the original functions answer, and the
  vector of indexes is not allocated, so matches with few HSPs (most blastp and blastx matches) cost what they
  did: a first form that wrapped every query's list in a 128-byte set (built, converted and dropped per redo call)
  added 0.27 % to the instructions of blastp WSSV x PajaWSV, about 50,000 redo calls (callgrind: the redo function
  itself, `BTreeSet` drops, `memcpy`, `malloc`); the lists-plus-indexes form is within 0.02 % of the base there,
  with the switch on and off. Conditional: the new side is a match where some query's list reaches 32 alignments
  (`x_reach!`, once per match), the reference side one where none does.
  Measured 2026-10-25 on the workstation (2x Xeon Gold 6338, cores warmed 3 s before each run, pairs alternating
  in order, `perf stat` walls; base 51234c2a; every output of a case the same bytes in every arm; tblastn options
  `-evalue 1e-5 -max_target_seqs 5 -outfmt 6`). Round 1, this pull request's first form ca0688e1 (the same index;
  every query's list wrapped in a 128-byte set): the first batch of HsG_b4 (35 human proteins x `-db hsap_dna`,
  396 matches) at 8 threads 330.8 -> 204.0 s (x0.619, 3 pairs, node 0 cores 8-15 with their siblings idle, load
  9-58; user 2,358 -> 1,477 s, instructions x0.934), at 16 threads 278.4 -> 185.4 s (x0.670, 3 pairs, node 1), at
  32 threads 269.3 -> 173.6 s (x0.650, 3 pairs, node 1); with `TBNREDOLPT` allowed up to 64 threads, 32 threads
  ran x0.991 against this form (3 pairs), so its limit of 16 stayed (it came off after `TBNREAPIDX`: see
  `TBNREDOLPT`). Matches with few HSPs (node 1, 8 pairs at 8 /
  16 threads): 16 proteins x E. coli x0.991 / x1.001, fly proteins x 5,000 16S x1.008 / x1.001, the same around fly
  2L (`-subject`) x0.998 / x1.006, against the fly genome x1.007 / x1.009, one protein against the zebrafish genome
  x1.001 / x0.992, instructions x1.0000 in every case. blastx of 50 kb of fly 2L x fly proteins (`-db`) at 8
  threads x1.000 (8 pairs). Round 2, the head 87dbbb02 (node 0 cores 16-23, load 30-89 while another gate ran):
  the 35 proteins at 8 threads 359.7 -> 209.4 s (x0.587, 3 pairs); blastp WSSV x PajaWSV at 1 thread x0.997
  (8 pairs), instructions x1.0001 (`=0` x1.0000; the first form +0.27 %); the blastx case x0.990 (8 pairs),
  instructions x1.0000; at 1 thread 16 proteins x E. coli and fly proteins x fly genome 81.420 G / 127.235 G
  instructions in all three arms (base, default, `=0`). Max RSS within 3 % of the base everywhere.
- `LOSAT_X_TBNKEEPFEW=1` tblastn (memory; in the strict set from the pull request that added it,
  2026-10-19, with `TBNFRAMEREUSE` (graduated 2026-10-27); `=0` turns it off; `tblastn/search_seed.rs`
  `XSubjectSide::kept_within_threads`, `tblastn/stage_d_pipeline.rs` `x_searched_side_by_side`): a
  search of four or more query batches that searches its batches one after another
  keeps the subjects' six-frame translations (the plans) only when they hold no more than its
  threads could have in flight, the thread count times the plan of the longest subject. One genome
  or a few chromosomes are kept; many subjects are translated again for each batch and dropped
  after the subject, as NCBI does (`blast_engine.c:627-631`). A search that searches batches side
  by side (at most 24 Mb of subjects, or fewer than two subjects per thread) or on pools
  of their own (`TBNOWNPOOLS`) keeps them as before: batches searched at the same time share them,
  and the batches after the first ones read them (5,000 fly cDNAs, 7 batches at 4 threads:
  4 side by side, then 3 one after another, took 1.004 of the time without them, 8 pairs). The
  layout of the batches does not change. Measured with `TBNFRAMEREUSE` (then a switch), the default against
  `TBNKEEPFEW=0 TBNFRAMEREUSE=0` (the default before, same binary; workstation, 2026-10-20, node 1,
  physical cores, booked cores warmed 3 s, arm order reversed every round, medians of 4-8 pairs;
  the fly proteins of `dmel_b4.faa`, five batches, `-evalue 1e-5 -max_target_seqs 5`): 57,775
  zebrafish cDNAs (122 Mb, `-db`) x1.001 at 1 thread (0.986-1.034), x0.989 at 8, `-subject` x0.990
  at 8; 28,036 16S records (40.6 Mb, `-db`) x0.985 at 8; 7 batches (`dmel_b7.faa`) against 5,000
  fly cDNAs (14.3 Mb) x1.000 at 2, x0.996 at 4 (both side by side, kept) and x0.998 at 8 (one
  after another); every report the same. Peak anonymous memory (RssAnon, one run each, node 0):
  the zebrafish cDNAs `-db` 361 -> 52 MiB at 1 thread and 390 -> 83 at 8 (NCBI BLAST+ 2.17.0 16
  and 28; the rest is the stored database), `-subject` at 8 417 -> 108 MiB (NCBI `-subject` holds
  every subject, about 800), the 16S set 134 -> 19 and 152 -> 38 MiB (NCBI 9 and 11), the 5,000
  cDNAs at 8 threads 70 -> 37 MiB (NCBI 16). Before (dev 23db0a4f, `TBNSSIDE=0` (a variable that graduated 2026-10-22) against the
  default, 38 pairs): x0.980-0.994 at 4, 8 and 16 threads on the two large sets, no pair slower.
- `LOSAT_X_TBNPAIRCAP=1` tblastn (memory; in the strict set from the pull request that added it, 2026-10-26;
  `=0` turns it off; `tblastn/search_seed.rs` `x_scan_plan`, `X_PAIR_ROOM`): the offset-pair array that a unit's scan
  calls fill (one per running unit; the units a thread runs one after another share theirs, `XUnitBuffers`) holds
  32,768 + the longest chain of the lookup table, never more
  than the 1024 x the longest chain it held before (a table whose longest chain is below 33 keeps that). NCBI's array
  holds 4,096 + the longest chain (`lookup_wrap.c:264-266` `GetOffsetArraySize`; `blast_engine.c:992,1040-1041`
  allocates one per thread, `aa_ungapped.c:504-505` passes its size to scansub). The size only decides where a scan
  call stops on a full array; the next call resumes there, so the pairs reach the two-hit loop in the same order at
  any size (unit tests `x_scan_bf_equals_ncbi_scan_on_random_subjects`, which also compares the pairs of all sizes in
  scan order, and the two scan tests of `x_scan_byte.rs`, at 4,096 and 32,768 + the longest chain among others;
  `tests/tblastn_scan_byte.rs` `the_offset_array_size_changes_only_the_scan_calls`: 115 proteins x MelaMJNV, the same
  report with and without the switch at one and four threads, 646 scan calls against 339). The
  1024 x of an early diagnostic was 864 KiB for the fly proteins of `dmel_b4.faa` (longest chain 108), against NCBI's
  33 KiB, and the thread count multiplied it. NCBI's own size costs time here: a call that stops starts the rest of its
  block of words again (`TBNSCANBYTE`), and at one thread `dmel_b4.faa` x zebrafish cDNAs (`-db`) took 1.0041 of the
  instructions with 4,096 + the longest chain, 1.0016 with 8,192, 1.0004 with 16,384 and 1.0001 with 32,768 (the
  outputs the same; with 4,096 the time at 8 / 16 threads was x1.0032 / x1.0076, 8 pairs). Always: every tblastn scan
  plan whose longest chain is 33 or more.
  Measured 2026-10-26 on the workstation (2x Xeon Gold 6338; base dev ce281c9a; every output of a case the same bytes in
  every arm). Peak memory (max of `/usr/bin/time` %M and a 2 ms VmHWM sample, node 0, untimed, medians of 3 rounds in
  rotating order, NCBI BLAST+ 2.17.0 in the same rounds): `dmel_b4.faa` x zebrafish cDNAs (`-db`, `-evalue 1e-5
  -max_target_seqs 5`) at 1 / 16 / 32 / 40 / 48 threads 28 / 89 / 151 / 176 / 207 -> 26 / 69 / 104 / 122 / 139 MiB
  (NCBI 118 / 137 / 149 / 156 / 162); the 595 proteins of AP027131 x AP027133 (`-db_gencode 4`) `-db` at 16 / 32 / 48
  63 / 97 / 114 -> 60 / 88 / 107 MiB (NCBI 115 / 117 / 116), `-subject` 62 / 90 / 120 -> 55 / 74 / 93 MiB (NCBI
  `-subject` runs one thread: 81). Time (node 1, cores warmed 3 s before each run, arms alternating in
  order, `perf stat`): the `dmel_b4.faa` case at 1 thread x1.0000 (8 pairs; instructions x1.00006, `=0` x1.00002); at
  8 / 16 / 32 threads x0.9994 / x1.0025 / x0.9979 against the base (16 rounds), where `=0` in the same binary, the code
  path of the base, read x1.0001 / x1.0055 / x1.0039 (the layout of the new binary), so that the switch against `=0`
  in one binary is x0.9986 / x0.9978 / x0.9900 (11 / 11 / 12 of 16 faster; instructions x1.00002); AP027131 `-subject`
  x0.9963 at 1 thread (8 of 8 faster) and x0.9986 at 16.
- `LOSAT_X_TBNDIAGPACK=1` tblastn (memory and speed; in the strict set from the pull request that added it,
  2026-10-29; `=0` turns it off; `tblastn/search_init.rs` `XDiagCells`, `x_packs`, `tblastn/x_tbn_hitbf.rs`
  `XDiagCell`): the diagonal table of the two-hit word finder (one per running unit; the units a thread runs one
  after another share it, `XUnitBuffers`) holds NCBI's `DiagStruct` in four bytes (`blast_extend.h:51-55`:
  `last_hit` in bits 0-30, `flag` in bit 31) instead of eight-byte `(last_hit, flag)` pairs: 128 KiB instead of 256
  KiB for a query batch of 20,000 letters. The values fit: `last_hit` starts at 0 and is set to `s + offset` or
  `s_last_off - (wordsize - 1) + offset`, the offset runs from the window to below INT4_MAX / 4 + the window + the
  longest unit (at most 5,000,000 residues), so with a window above the word size and at most 2^28 every value is
  below 2^30 and not negative, inside NCBI's signed 31 bits too; another window keeps the pairs (LOSAT's TBLASTN
  takes window 40 with word size 3 and 60 with 2, both packed). A flag-0 cell read as `i32` is then
  its `last_hit` and a flag-1 cell is negative, so the branch-free loop (`x_hits_bf`) reads one word per hit. The
  tables give the same init HSPs and the same values in every cell (unit tests
  `x_hits_bf_on_packed_cells_equals_pair_cells`, and `x_pairs_bf_gives_the_hits_and_table_of_pairs` on both
  layouts against `pairs` on pairs). Always: every tblastn search with the default window.
  Measured 2026-10-29 on the workstation (2x Xeon Gold 6338; base dev 593fff35; every output of a case the same
  bytes in every arm). Peak memory (max of `/usr/bin/time` %M and a 2 ms VmHWM sample, node 0, untimed, at most 16
  CPUs for every run, medians of 4 rounds in rotating order, NCBI BLAST+ 2.17.0 in the same rounds): the 595 proteins
  of AP027131 x AP027133 (`-db_gencode 4`) `-subject` at 16 / 32 / 40 / 48 threads 56 / 75 / 83 / 94 -> 55 / 72 / 80
  / 88 MiB (NCBI `-subject` runs one thread: 81), `-db` at 48 104 -> 101 (NCBI 117); `dmel_b4.faa` x zebrafish cDNAs
  (`-db`) at 16 / 32 / 40 / 48 69 / 103 / 120 / 138 -> 67 / 101 / 116 / 132 (NCBI 137 / 150 / 156 / 162); with `=0`
  in the same binary 56 / 75 / 83 / 93, 105, 68 / 104 / 120 / 137 (the base's). At one thread (`perf stat` of a
  probe build with the same cells, 3 runs each, node 0, untimed): AP027131 `-subject` instructions x0.99995, cycles
  x0.989; the `dmel_b4.faa` case instructions x1.00004, cycles x0.946 (64.8 -> 61.2 s): half the table's bytes on
  the random reads of the two-hit test. Time (node 1, cores warmed 3 s before each run, arms alternating in order,
  `perf stat`; the PR's binary against the base, and `=0` in the PR's binary): `dmel_b4.faa` at 1 / 8 / 16 / 32
  threads x0.953 / x0.955 / x0.956 / x0.958 (8 of 8, 12 of 12 each faster; instructions x0.9987; `=0` x0.9992 /
  x1.0005 / x1.0002 / x0.9998); AP027131 `-subject` at 1 / 16 / 32 x0.991 / x0.995 / x0.990 (12, 7, 8 of 12 faster;
  `=0` x1.0001 / x1.0007 / x1.0000).
- `LOSAT_X_TBNSCANRESUME=1` tblastn (memory; in the strict set from the pull request that added it, 2026-10-29;
  `=0` turns it off; `tblastn/x_scan_byte.rs` `XScanHeld`, `x_scan_one_range_held`, `tblastn/search_seed.rs`
  `x_scan_unit`, `OFFSET_ARRAY_SIZE`): a unit scanned on presence bytes (`LOSAT_X_TBNSCANBYTE`) keeps the block of
  words where a scan call stops on a full offset array (the block's indexes and candidates and the candidate that did
  not fit), and the next call of the same range copies on from that candidate instead of making the block's
  candidates again; the offset array of such a unit then holds NCBI's 4,096 + the longest chain
  (`lookup_wrap.c:264-266` `GetOffsetArraySize`) instead of `LOSAT_X_TBNPAIRCAP`'s 32,768 + it (never more than the
  plan's size): 33 KiB instead of 257 KiB per running unit for the fly proteins of `dmel_b4.faa` (longest chain 108).
  A stop then costs one call of the scan and of the two-hit loop, not a block of up to 256 words made again (which
  took 1.0041 of the instructions at this size, `LOSAT_X_TBNPAIRCAP`). The pairs and the `s_range` of every call are
  NCBI's at any size (unit tests `x_scan_byte_equals_ncbi_scan_on_random_subjects` and
  `the_last_overflow_chain_is_copied_exactly`, both with and without the kept block; the first counts that the call
  after every stop finds its block kept for it; `tests/tblastn_scan_byte.rs` `the_offset_array_size_changes_only_the_scan_calls`: the same
  report at one and four threads with both switches, `=0` and neither). The units of the branch-free scan
  (`x_scan_bf.rs`: aarch64, and tables that are not sparse without AVX2) keep `LOSAT_X_TBNPAIRCAP`'s size. Where:
  the units that take TBNSCANBYTE's scan (platform).
  Measured 2026-10-29 on the workstation (2x Xeon Gold 6338; base dev d6208db9; every output of a case the same bytes
  in every arm). Peak memory (max of `/usr/bin/time` %M and a 2 ms VmHWM sample, node 0, untimed, at most 16 CPUs for
  every run, medians of 8 rounds in rotating order, NCBI BLAST+ 2.17.0 in the same rounds): the 595 proteins of
  AP027131 x AP027133 (`-db_gencode 4`) `-subject` at 1 / 16 / 32 / 40 / 48 threads 19 / 55 / 72 / 80 / 88 -> 18 / 53
  / 69 / 76 / 84 MiB (NCBI `-subject` runs one thread: 81), `-db` at 48 101 -> 94 (NCBI 115); `dmel_b4.faa` x
  zebrafish cDNAs (`-db`, `-evalue 1e-5 -max_target_seqs 5`) at 16 / 32 / 40 / 48 66 / 101 / 116 / 133 -> 62 / 90 /
  102 / 117 (NCBI 137 / 149 / 156 / 162); lower in 79 of the 80 pairs; `=0` in the same binary x0.990-1.008 of the
  base. The `dmel_b4.faa` case takes 3.08 M scan calls instead of 1.40 M (1.39 M units): at one thread x1.0019 of the
  instructions (about 380 per extra stop, the calls of the scan and of the two-hit loop), AP027131 `-subject` x1.0003;
  blastx (`BXSCANBYTE`, AP027131 x AP027133 proteins and 800 fly cDNAs x yeast proteins `-subject`) x1.00000. Time
  (node 1, cores warmed 3 s before each run, arms alternating in order, `perf stat`; the PR's binary against the base,
  and `=0` in the PR's binary): `dmel_b4.faa` at 1 / 8 / 16 / 32 threads x0.9962 / x0.9974 / x0.9958 / x0.9991 (7 of
  8, 11, 11 and 7 of 12 faster), `=0` x0.9919 / x0.9915 / x0.9920 / x0.9910; AP027131 `-subject` at 1 / 16 / 32
  x0.9986 / x1.0014 / x0.9985 (10, 4, 6 of 12 faster), `=0` x0.9974 / x0.9972 / x0.9999. Against `=0` in the same
  binary the switch costs x1.004-1.008 on `dmel_b4.faa` (the extra stops), which the binary's layout hides against
  the base.
- `LOSAT_X_TBNLUTLEAN=1` tblastn (memory; in the strict set from the pull request that added it, 2026-10-30;
  `=0` turns it off; `tblastx/lookup/x_lut_direct.rs` `x_build_lean`, `tblastx/lookup/backbone.rs`
  `x_lookup_table_for_profile`): a query batch's lookup table is built with the passes of the in-place builder of
  BLASTX, BLASTP and TBLASTX (`x_build_direct`: count, layout, fill, no chains) with every temporary freed by the
  build: the exact backbone as a counting sort of the query's words, the neighbour lists of its words made once into
  one vector (no memo kept by the thread), and the per-cell counts reused as the fill cursors. The table's backbone
  and overflow are mapped (on 64-bit glibc a capacity above its mmap ceiling, `utils/x_mapped.rs`
  `x_unthresholded_capacity`; only the pages written are resident; elsewhere the capacity `Vec` gives), so that a
  batch's table goes back when its search ends.
  The chain builder kept the exact and thin backbones' cells (446 KiB each), their chains' slabs (about 1.4 MiB) and
  the counts next to the new table, 3.4 MiB at once for a batch of 18,000 letters; the lean build peaks at about
  1.5 MiB. TBLASTN builds each batch's table on the thread of the batch's first subject, and a thread's glibc arena
  keeps the high water of what it allocated: with the batches searched side by side (`-subject`) the difference
  stayed in each building thread's arena under what its search allocated later. With the batches searched one after
  the other (`-db`) glibc trimmed the chain builder's large free top with the table at the batch's end; the lean
  build's smaller one is not trimmed, and a table left in the arena kept the peak 1-3 % above the chain builder's
  until the table was mapped.
  The table is the chain builder's (unit test `lean_table_equals_reference_for_thresholds_masks_matrices_and_word_lengths`:
  BLOSUM62 word lengths 2 and 3, BLOSUM45 at 2 and 3, PAM30, thresholds 0-1000, masks, the overflow array and a query
  set without words; `tests/tblastn_lut_lean.rs`: the same report as `=0` with BLOSUM62 word size 3 and BLOSUM45 word
  size 2, at one and four threads, every table lean by default and with chains with `=0`, read from
  `[X_STATS] TBNLUTLEAN lean= chain=` with `LOSAT_X_STATS`). Always: every TBLASTN lookup table (LOSAT's TBLASTN takes
  word size 3, and 2 with BLOSUM45) outside diagnostics mode, whose prints read the chains; a table of 4 or more
  letters a word (backbone indexes over 16 bits, which the CLI refuses) would take the chain builder, as the in-place
  builders decline it.
  Measured 2026-10-30 on the workstation (2x Xeon Gold 6338; base dev a2ee6ff3; every output of a case the same
  bytes in every arm). Peak memory (max of `/usr/bin/time` %M and a 2 ms VmHWM sample, node 0, untimed, at most 16
  CPUs for every run, medians in rotating order, NCBI BLAST+ 2.17.0 in the same rounds): the 595 proteins of AP027131
  x AP027133 (`-db_gencode 4`) `-subject` at 1 / 16 / 32 / 40 / 48 threads 17 / 53 / 68 / 76 / 83 -> 17 / 45 / 60 /
  67 / 75 MiB (8 rounds, lower in every pair; NCBI `-subject` runs one thread: 81), `-db` at 48 95 -> 93 (NCBI 115);
  `dmel_b4.faa` x zebrafish cDNAs (`-db`, `-evalue 1e-5 -max_target_seqs 5`) at 16 / 32 / 40 / 48 62 / 90 / 102 / 116
  -> 61 / 90 / 101 / 116 (x0.985 / x0.994 / x0.997 / x0.994, 16 rounds, lower in 8 / 12 / 9 / 11 of the pairs; NCBI
  128 / 139 / 145 / 151); `=0` in the same binary x0.994-1.012 of the base. At one thread x0.9957 of the instructions
  (AP027131 `-subject`) and x0.9998 (`dmel_b4.faa`); BLASTX, BLASTP and TBLASTX x1.00000 (their code is the same).
  Time (node 1, cores warmed 3 s before each run, arms alternating in order; the PR's binary against the base, and
  `=0` in the PR's binary): AP027131 `-subject` at 1 / 16 / 32 threads x0.9949 / x0.9986 / x0.9840 (11, 7, 10 of 12
  faster), against `=0` x0.9937 / x0.9957 / x0.9844; `dmel_b4.faa` at 1 / 8 / 16 / 32 x0.9986 / x0.9963 / x0.9965 /
  x0.9946 (6 of 8, 11, 9, 11 of 12 faster), against `=0` x0.9989 / x0.9987 / x0.9984 / x0.9919.
- `LOSAT_X_NEWTONDROP=1` blastp, tblastn (memory; in the strict set from the pull request that added it, 2026-11-01;
  `=0` turns it off; `core/composition_adjustment/x_newton_lanes.rs` `prefetch`): the Newton lanes (NEWTONLANES)
  keep the results they solve ahead in a store of the thread (6,736 bytes a result: its input, matched by bits, and
  its final iterate), and a result not taken stayed there until the store held 512 (3.3 MiB). A prefetch comes after
  the windows (TBLASTN: the coming windows of a match) or the matches (BLASTP: the next 16 matches of a query, or a
  job's matches in the windowed redo) of the thread's earlier prefetches were walked, so their results not taken are
  not asked for again; with the switch a prefetch drops them first. TBLASTN leaves such results for the HSPs its redo
  finds contained (their problems are solved ahead with the others): 1,667 of the 50,620 results solved ahead for
  `dmel_b4.faa` x zebrafish cDNAs (`-db`, 10.7 MiB of stores when four or more threads share them; one thread keeps
  512 at most), and 44 of 40,738 for AP027131 x AP027133; BLASTP left none (WSSV x PajaWSV, 3,000 SwissProt proteins
  x Species6 `-db`). A dropped result asked for after all would be solved by its
  call with the same bits (`take` returns none), so the report cannot change; the results taken are the same
  (`[X_STATS] NEWTONLANES prefetched= hits= misses=` equal with and without the switch; `[X_STATS] NEWTONDROP
  dropped=` counts the results dropped). Unit test `x_newton_drop_empties_the_store_at_the_next_prefetch` (both sides
  in one process: the same results taken, the earlier ones not taken dropped or kept);
  `tests/tblastn_newton_drop.rs`: the same report as `=0` at one and four threads, the same counts of results solved
  ahead and taken, some results dropped by default and none with `=0`. Always: every prefetch (where NEWTONLANES
  runs).
  Measured 2026-11-01 on the workstation (2x Xeon Gold 6338; base dev b94e9844; every output of a case the same
  bytes in every arm). Peak memory (max of `/usr/bin/time` %M and a 2 ms VmHWM sample, node 0, untimed, at most 16
  CPUs for every run, medians in rotating order, NCBI BLAST+ 2.17.0 in the same rounds): `dmel_b4.faa` x zebrafish
  cDNAs (`-db`, `-evalue 1e-5 -max_target_seqs 5`) at 16 / 32 / 40 / 48 threads 61 / 90 / 102 / 116 -> 54 / 85 / 96 /
  111 MiB (x0.890 / x0.943 / x0.940 / x0.953, 16 rounds, lower in every pair; NCBI 128 / 140 / 145 / 152; the 10 ms
  RssAnon series lower from 40 % of the run on); the 595 proteins of AP027131 x AP027133 (`-db_gencode 4`) `-subject`
  at 1 / 16 / 32 / 40 / 48 x1.000 / x0.990 / x0.992 / x0.998 / x0.997 (8 rounds; 75 MiB at 48, NCBI 81 at one
  thread), `-db` at 48 x0.996 (16 rounds, lower in 9; 93 MiB, NCBI 113); BLASTP WSSV x PajaWSV `-subject` at 16
  x0.995; `=0` in the same binary x1.000-1.012 of the base. At one thread x0.9998 of the instructions (AP027131
  `-subject`) and x0.9996 (`dmel_b4.faa`); BLASTP x1.00002, BLASTX and TBLASTX x1.00000.
  Time (node 1, cores warmed 3 s before each run, arms alternating in order; the PR's binary against the base, and
  `=0` in the PR's binary): AP027131 `-subject` at 1 / 16 / 32 threads x0.9992 / x1.0000 / x0.9956 (8, 6, 7 of 12
  faster), against `=0` x0.9992 / x1.0000 / x1.0043; `dmel_b4.faa` at 1 / 8 / 16 / 32 x0.9998 / x0.9989 / x1.0006 /
  x1.0011 (5 of 8, 7, 5, 5 of 12 faster), against `=0` x0.9984 / x0.9983 / x0.9986 / x0.9989.
- `LOSAT_X_NEWTONBATCH=1` blastp, tblastn (memory; in the strict set from the pull request that added it, 2026-11-01;
  `=0` turns it off; `core/composition_adjustment/x_newton_lanes.rs` `prefetch` / `take`): with NEWTONDROP the results
  a thread solved ahead are the ones of its last prefetch only, so a prefetch keeps them in the input and output
  vectors they were solved in (one batch a thread; `take` searches it from the first result not taken, and the batch
  goes when its last result is taken or at the next prefetch) instead of copying each into a box of the store (6,736
  bytes a result). The copies held the caller's inputs, the outputs and the boxes at once at the end of every prefetch
  (about 13.2 KiB a problem; TBLASTN's redo solves the problems of the coming windows of a match until 32 or more are
  gathered), and this transient set the high water of most threads' heaps in the redo. The results taken are the
  same (`[X_STATS] NEWTONLANES prefetched= hits= misses=` and `[X_STATS] NEWTONDROP dropped=` equal with and without
  the switch, unless one prefetch solves more than the 512 results the store holds: the store keeps the last 512 of
  them, the batch all, and more are taken; `[X_STATS] NEWTONBATCH kept=` counts the prefetches kept as their batch);
  with `LOSAT_X_NEWTONDROP=0` the results are stored (the batch holds one prefetch). The search of the batch and its comparison of an input are
  functions of their own: inlined into `composition_matrix_adjust` they changed its code on every path (one-thread
  instructions x1.001 with the switch off), and the comparison inlined into the search kept the index of its loop
  over `q` on the stack (twice the instructions of the store's search). Unit tests
  `x_newton_batch_gives_each_result_once_and_goes_when_all_are_taken` (any order, a repeated input, the store empty,
  the batch gone with its last result) and `x_newton_drop_empties_the_store_at_the_next_prefetch` (both forms in one
  process); `tests/tblastn_newton_drop.rs`: the same report as `=0` at one and four threads, the same counts of
  results solved ahead, taken and dropped, batches kept by default and none with `=0` or `LOSAT_X_NEWTONDROP=0`.
  Always: every prefetch (where NEWTONLANES runs, with NEWTONDROP).
  Measured 2026-11-01 on the workstation (2x Xeon Gold 6338; base dev 7a835d1f; every output of a case the same
  bytes in every arm). Peak memory (max of `/usr/bin/time` %M and a 2 ms VmHWM sample, node 0, untimed, at most 16
  CPUs for every run, medians in rotating order, NCBI BLAST+ 2.17.0 in the same rounds): `dmel_b4.faa` x zebrafish
  cDNAs (`-db`, `-evalue 1e-5 -max_target_seqs 5`) at 16 / 32 / 40 / 48 threads 55 / 84 / 96 / 110 -> 54 / 81 / 92 /
  106 MiB (x0.981 / x0.966 / x0.959 / x0.963, 16 rounds, lower in 13, 15, 16, 16; NCBI 128 / 139 / 145 / 151); the
  595 proteins of AP027131 x AP027133 (`-db_gencode 4`) `-subject` at 1 / 16 / 32 / 40 / 48 x0.998 / x0.992 / x1.002 /
  x1.002 / x0.989 (8 rounds; 74 MiB at 48, NCBI 81 at one thread), `-db` at 48 x0.994 (93 MiB, NCBI 114); BLASTP WSSV
  x PajaWSV `-subject` at 16 x0.992; `=0` in the same binary x0.999-1.005 of the base. At one thread x0.99922 of the
  instructions (AP027131 `-subject`) and x0.99999 (`dmel_b4.faa`); BLASTP x0.99960, BLASTX and TBLASTX x1.00000.
  Time (node 1, cores warmed 3 s before each run, arms alternating in order; the PR's binary against the base, and
  `=0` in the PR's binary): AP027131 `-subject` at 1 / 16 / 32 threads x0.9954 / x0.9942 / x0.9941 (12, 6, 8 of 12
  faster), against `=0` x0.9961 / x0.9986 / x1.0015; `dmel_b4.faa` at 1 / 8 / 16 / 32 x0.9992 / x0.9984 / x1.0008 /
  x0.9992 (5 of 8, 10, 6, 6 of 12 faster), against `=0` x1.0003 / x0.9999 / x1.0001 / x1.0004.
- TBNLINKWIN (graduated 2026-10-26, in the set since 2026-10-17, by the Owner's rule of that day that a
  speed-up lands on by default, on the timings below; no variable; always on) tblastn (speed;
  `tblastn/stage_d_linking/x_link_window.rs`): the uneven-gap linking of a subject's HSP list
  (`stage_d_linking.rs` `link_preliminary_hsps`, NCBI's `s_BlastUnevenGapLinkHSPs`) tests only the
  candidates of a head set that can be admissible. NCBI tests every HSP in the head's query-offset
  window; with a short query that is nearly the whole list (a 411-residue protein x the human
  genome: 85,908 HSPs, 159 million admissibility tests, 13 links). An admissible merged set has
  increasing subject starts, each at most `longest_intron` after the end before it, so a set whose
  subject range lies more than `longest_intron` from the head set's cannot be admissible. The HSPs
  alone in their set are looked up by context, frame sign and subject start, the HSPs of linked sets
  are always kept (their ranks sorted, so a head takes those in its window by two binary searches;
  since `TBNLINKSETS`, 2026-10-26, only the sets the lookup meets, unless `LOSAT_X_TBNLINKSETS=0`),
  and the candidates are tested in NCBI's order with NCBI's tests; a window with no more HSPs than
  the subject starts looked up is tested whole, as NCBI does, so every link, number and E-value is
  NCBI's (unit test `x_chosen_candidates_give_the_links_of_the_whole_window` against the whole
  window on random lists; the whole-window loop stays in `link_preliminary_hsps` for that test only).
  With `LOSAT_X_STATS=1` a line `[X_STATS] TBNLINKWIN heads= window= chosen= sets=` gives the heads, the
  candidates of their windows, the candidates tested and the sets the lookup brought in. blastx
  keeps its own linking.
  Timed on the workstation (perf-si-20261012 session 21, node 1, NCBI BLAST+ 2.17.0 `-db` and LOSAT
  with TBNSCANBF (then a switch) without and with the switch in turn, three rounds, one hash per case;
  time with the switch / without): zebrafish genome x one protein 0.60 at one thread (24.8 -> 14.8
  s, NCBI 31.3 s); human genome x one protein 0.87 at one thread (32.3 -> 28.0 s, NCBI 44.9 s) and
  at 8 threads (5.10 -> 4.46 s); fly genome x one batch at 8 threads 0.76 (4.04 -> 3.06 s); E. coli
  MG1655 x 16 batches at 8 threads 0.95. Again in session 22 (node 0, the default with all three
  tblastn switches of 2026-10-17 against `=0` in rotating rounds, 3-8 rounds, one hash per case):
  zebrafish 0.58 (23.2 -> 13.5 s), human 0.82 at one thread and 0.84 at 8 threads, fly 0.81, MG1655
  0.98.
- `LOSAT_X_TBNLINKCTX=1` tblastn (speed; in the strict set from the pull request that added it, 2026-10-27;
  `=0` turns it off; `tblastn/stage_d_linking/x_link_ctx.rs` `x_link_by_context`): a subject's HSPs are linked one
  query context at a time, and a subject whose units ran on the pool (the single subject of a search, or one longer
  than a chunk's budget) links its contexts there side by side, those with the most HSPs first. Before, the whole list was linked in one pass on the
  thread that merged the subject's units: with many queries the long subjects ended in a serial link while the
  other threads waited (the B11e / B10-E probes, 35 human proteins x the human genome `-db`: chr1 936 k HSPs over
  35 contexts, link 1.38-1.47 s, of it 1.03 s the head loop, 0.23 s sorts; the six longest chromosomes ended the
  preliminary stage in links of 0.8-1.5 s with 2-6 threads busy, at 16 and at 32 threads). The chains of NCBI's
  pass never leave a context (a head's candidates are HSPs of its context, `link_hsps.c` 1048-1060 and 1407-1491),
  so each context's HSPs linked alone, in the sum-score order restricted to the context, give the same chains,
  numbers and e-values; the values go back to the list's order before the final stable sort by score (unit test
  `contexts_linked_apart_give_the_list_of_one_pass`: random lists of 1 to 3,000 HSPs over 1 to 35 contexts, on one
  thread and on a pool of 4, against the one pass, every field and the e-value bits). A list of one context is one
  call as before (single-query searches and the post-redo linking of one match do not change). Always: every
  linked list of two or more contexts takes this path. Measured 2026-10-27 on the workstation (2x Xeon Gold 6338,
  node 1 whole, cores warmed 3 s before each run, pairs alternating in order, `perf stat` walls; the first commit of
  the pull request that added it, 1a4ef2a6, against its base c780af86; every output of a case the same bytes in
  every arm; ratio = median of the per-pair ratios; `-evalue 1e-5 -max_target_seqs 5`): 35 human proteins x the
  human genome (`-db hsap_dna`) at 32 threads 29.78 -> 28.49 s (x0.963, 3 of 3 pairs; `=0` 29.86 s), at 16 threads
  41.35 -> 40.20 s (x0.972, 2 of 2; the probe model put -1.17 s and -1.15 s, the runs -1.29 s and -1.15 s); 794
  proteins x E. coli MG1655 at 32 threads x0.989 (8 of 8), 27 fly proteins x the fly genome x0.974 (6 of 8). The same
  where a list has one context or the subjects are short: one human protein x the human genome at 16 threads x0.998
  (5 of 8 faster), one zebrafish protein x the zebrafish genome at 32 threads x1.009 (2 of 8; per pair 0.98-1.08,
  instructions x1.0001), the 35 proteins x human cDNA (`-db hsap_cdna`, short subjects, each linked on its own
  thread) at 16 threads x1.001 (instructions x1.0019, the grouping by context). One thread (node 0, one core):
  794 x MG1655 9.58 -> 9.53 s, 27 fly proteins x the fly genome 15.46 -> 15.34 s, instructions x0.9994 and x0.9996.
  Max RSS (`time` %M): 35 x human genome x0.994 at 32 threads and x0.991 at 16, x human cDNA 143 -> 121 MiB (x0.85,
  3 of 3), other cases x0.98-1.01 at 32 threads; at one thread the fly case reads 49.0 or 54.6 MiB in either arm
  (the A/B's pair 49.0 -> 54.6 MiB, a later run of 3 pairs on the same core 54.6 -> 49.0 MiB in each pair).
- `LOSAT_X_TBNLINKSETS=1` tblastn (speed; in the strict set from the pull request that added it, 2026-10-26;
  `=0` turns it off; `tblastn/stage_d_linking/x_link_window.rs` `XLinkWindow::candidates`): the linking candidates
  of a head whose query-offset window holds more HSPs than its subject lookup (TBNLINKWIN, graduated) take a set of
  two or more only when the lookup meets one of the set's HSPs, and then every HSP of that set in the window. Before,
  every HSP of every set in the window was a candidate: repeat-family queries against a genome form thousands of
  small sets, and every head tested all of them again (the B11i probe, 35 human proteins x the human genome `-db`
  at 8 threads: 12.47 M heads chose 18.85 G candidates, 18.83 G of them HSPs of sets, and of the 9.60 G set HSPs
  tested past NCBI's skip of 1719, 0.47 M (0.005 %) belonged to a set whose subject range meets the head set's
  range +- `longest_intron`; chr1 x one query: 13,272 sets of at most 3, 11,734 candidates per head). A set of two or
  more is itself an admissible merge, so each gap inside it is at most `longest_intron`; when its subject range
  meets [lowest start - `longest_intron`, highest end + `longest_intron`] of the head set (which every admissible
  set does), one of its HSPs meets that range, and the lookup (subject starts within `longest_intron` plus the
  longest HSP of the list) finds it. The set is brought in once, by its first HSP met (the starts of a set increase
  along its chain); the HSPs of a set share its sum score, number and E-value, and all of its HSPs in the window
  stay candidates in NCBI's order, so NCBI's skip (`link_hsps.c` 1719) and tests (1407-1491) choose the same link
  (unit test `x_chosen_candidates_give_the_links_of_the_whole_window`: random lists and lists of repeat copies of two
  to five pieces over the whole subject, switch on and off, against NCBI's whole window; `tests/tblastn_link_window.rs`:
  the 115 proteins of PemoMJNVA x MelaMJNV, the switch on and off at 1 and 2 threads). The sorted list of set HSPs
  and its upkeep after each link are not kept while the switch is on. With `LOSAT_X_STATS=1` the
  `[X_STATS] TBNLINKWIN` line ends with `sets=`, the sets the lookup brought in (0 with `=0`). Always: every
  linking call whose heads take the lookup (a head's window longer than its lookup) takes this path.
  Cost and floor: the preliminary linking was 902 of the 1,405 thread-s of the preliminary stage (35 proteins at 8
  threads, session 10 of v4-B10); the floor is the lookup, about 1 us per head (298 thread-s of lookup and 631 of
  tests for 12.47 M heads, loaded). Measured 2026-10-26 on the workstation (2x Xeon Gold 6338, cores warmed 3 s
  before each run, pairs alternating in order, `perf stat` walls; base 9d477f4f; every output of a case the same
  bytes in every arm; options `-evalue 1e-5 -max_target_seqs 5 -outfmt 6`): the first batch of HsG_b4 (35 human
  proteins x `-db hsap_dna`) at 8 threads 212.1 -> 90.8 s (x0.433, 3 pairs, node 0 cores 24-31 with their siblings
  idle, load 15-42; user 1,487 -> 706 s, instructions x0.515, max RSS x0.918; `=0` 211.0 s), at 16 threads 186.3 ->
  48.3 s (x0.259) and at 32 threads 176.1 -> 42.9 s (x0.243; 3 pairs each, node 1; max RSS x1.00-1.01); at 1 thread
  16 proteins x E. coli and fly proteins x the fly genome 81.416 / 81.404 / 81.416 G and 127.222 / 127.132 / 127.221 G
  instructions (base, default, `=0`). Light cases (node 1, 8 pairs at 8 / 16 threads, instructions x0.9985-0.9999,
  max RSS within 2 %): 16 proteins x E. coli x0.992 / x0.995, fly proteins x 5,000 16S x1.006 / x1.016 (0.10-0.15 s
  runs; per-round ratios 0.92-1.09), the same around fly 2L (`-subject`) x0.999 / x0.979, against the fly genome
  x0.995 / x0.991, one protein against the zebrafish genome x0.996 / x1.002.
- `LOSAT_X_TBNREAPIDX=1` tblastn (speed; in the strict set from the pull request that added it, 2026-10-26;
  `=0` turns it off; `tblastn/kappa/x_reap_index.rs` `x_reap_keep`, used by `kappa.rs`
  `reap_contained_postredo_hsps`): after the composition-based redo of a match, NCBI removes every HSP whose two end
  points lie inside an earlier HSP still kept, of the same context and frame, that scores at least as much
  (`s_HitlistReapContained`, `blast_kappa.c` 224-273; `CONTAINED_IN_HSP`, `blast_hits_priv.h` 65-70), testing each
  HSP against every earlier one, so a match with tens of thousands of HSPs costs the square of that: the B11h probe
  (35 human proteins x the human genome `-db` at 8 threads, session 10 of v4-B10) counted 21.11 G iterations, 3.48 ns
  each (73.6 thread-s, about 10 % of the user time once `TBNLINKSETS` had cut the linking). With the switch a list of
  32 HSPs or more (`X_REAP_INDEX_FROM`) keeps the HSPs it has kept so far in a `BTreeSet` of (context, frame, subject
  start, position), and each HSP tests only those whose subject start lies in [max(start, end) - dmax, min(start,
  end)] (dmax: the largest `s_end - s_start` kept so far), with NCBI's whole test, the score included: an HSP that
  contains both end points starts at or before both and ends at or after both, so it lies in that range. The answer
  is whether some such HSP exists, which does not depend on the order they are met in, so every kept HSP, the
  survivors' order and the report stay the same (unit tests `x_reap_keep_keeps_what_ncbi_keeps` (1,152 random lists
  of 2 to 800 HSPs: dense and sparse coordinates, score order and none, tied scores, nested HSPs, ends one past the
  edges, other frames and contexts, ends before starts, against NCBI's loop `reap_keep`), `x_reap_keep_edges`,
  `the_reap_is_the_same_with_and_without_the_index`; six broken variants of the index fail them: dmax never grown,
  the range one short at either end, the score not compared, frames of one sign sharing an index, the range start
  wrapped instead of clamped at the ends of `i32`; a seventh, which also indexes the HSPs already reaped, gives the
  same answers: an HSP inside a reaped one is inside the HSP that reaped it; `tests/tblastn_reap_idx.rs`: four
  AP027131 proteins with long redo lists and one with a short list against the AP027133 genome, the switch on and
  off at 1 and 2 threads, both sides reached). The probe's shadow of this index (the same range and order, lowest
  start first) visited 6,326 candidates in 0.30 thread-s with no answer differing. Shorter lists keep NCBI's loop.
  Conditional: the new side is a match whose redo gives 32 HSPs or more (`x_reach!`, once per match), the reference
  side one with 2 to 31. Cost and floor: the reap was 73.6 thread-s (the probe) of the 706-730 s of user time of the
  35 proteins at 8 threads after `TBNLINKSETS` (here user 730 -> 662 s); the floor is the index, about 0.3 thread-s.
  Measured 2026-10-26 on the workstation (2x Xeon Gold 6338, cores warmed 3 s
  before each run, pairs alternating in order, `perf stat` walls; base f5326a59; every output of a case the same
  bytes in every arm; options `-evalue 1e-5 -max_target_seqs 5 -outfmt 6`): the first batch of HsG_b4 (35 human
  proteins x `-db hsap_dna`) at 8 threads 93.7 -> 85.3 s (x0.904, 3 pairs, node 0 cores 24-31 with their siblings
  idle, load 57-87; user 730 -> 662 s, instructions x0.955, max RSS x0.999; `=0` one run 105.5 s at the same cycles
  as the base (x1.004) and instructions x1.0000, a lower clock under load), at 16 threads 47.2 -> 41.5 s (x0.877)
  and at 32 threads 41.7 -> 36.1 s (x0.865; 3 pairs each, node 1; max RSS x0.95-0.99); at 1 thread 16 proteins x
  E. coli and fly proteins x the fly genome 81.455 / 81.455 / 81.456 G and 127.219 G instructions in all three arms
  (base, default, `=0`). Light cases (node 1, 8 pairs at 8 / 16 threads, instructions x1.0000, max RSS within 3 %):
  16 proteins x E. coli x1.001 / x1.002, fly proteins x 5,000 16S x0.995 / x1.019 (no hits, 0.10-0.15 s runs;
  per-round ratios 0.97-1.08), the same around fly 2L (`-subject`) x0.995 / x1.003, against the fly genome x1.001 /
  x0.980, one protein against the zebrafish genome x1.001 / x0.999.
- `LOSAT_X_TBNXLATESIMD=1` tblastn (speed; in the strict set since 2026-10-17, by the Owner's rule
  of that day that a speed-up lands on by default, on the timings below: on by default, `=0` turns
  it off; `tblastn/x_subject_translate.rs`): a packed subject's six frames (`TBNSTORE`, every `-db`
  sequence) are read 32 codons of a phase at a time with AVX2: the 24
  bytes after bit 2 x phase (two loads one byte apart), every three bytes split into four six-bit
  codons (a byte shuffle and two multiplications, as base64 encoders split bytes), each 64-residue
  codon table of the table translation (`TBNXLATE`) read as four `vpshufb` tables chosen by bits 4 and 5 of the codon, and
  the minus residues stored in reverse order. Without AVX2, and with `=0`, the codons come eight at
  a time from one big-endian word of eight packed bytes (`x_phase_words`; NCBI's
  `BLAST_TranslateCompressedSequence` reads three bytes per four codons, `blast_util.c:584-596`);
  that loop replaced, with no switch, the unpacking of the codes into a buffer and the codon indices
  made from it in a second pass. Both give the frames of `x_generate_frames` on the store's codes
  (unit test `x_packed_frames_read_from_the_bytes_are_the_frames_of_the_codes`, every length from 0
  to 400 bases and five long ones, three genetic codes); `LOSAT_X_TBNXLATESHADOW=1` (gone) compared 6.2 G
  residues of the human genome with the reference translation, none different. With
  `LOSAT_X_STATS=1` a line `[X_STATS] TBNXLATESIMD codons_vector= codons_scalar=` counts the codons
  of each loop (`tests/tblastn_xlate_simd.rs`). perf stat on the workstation (one thread, a
  411-residue protein x the human genome `-db`, `TBNSCANBF` and `TBNLINKWIN` on, the same output):
  cycles 75.6 G before -> 72.0 G with the word loop -> 66.0 G with the vector loop (instructions 271
  -> 245 -> 219 G); zebrafish genome x one protein 41.8 -> 39.9 -> 37.4 G. Timed in session 22
  (node 0, the default against `=0` in rotating rounds, 3-8 rounds, one hash per case; time on /
  off): human genome x one protein 0.88 at one thread (27.9 -> 24.5 s; NCBI BLAST+ 2.17.0 `-db` 47.1
  s) and 0.94 at 8 threads; zebrafish genome x one protein 0.99; fly genome x one batch and MG1655 x
  16 batches at 8 threads 0.99-1.00.
- `LOSAT_X_TBXSCANBF=1` TBLASTX (speed; in the strict set since 2026-10-17, on by default; `=0` turns
  it off; `tblastx/x_scan_bf.rs`): the lookup scans of the unit path
  (`x_pair_par.rs`) and of the reference frame loop (`blast_engine/run_impl.rs`) walk the subject
  ranges as NCBI's `s_BlastAaScanSubject` (`blast_aascan.c:48-131`) does and scan each range with
  the branch-free loop of tblastn (TBNSCANBF, graduated; `tblastn/x_scan_bf.rs` `x_scan_one_range_bf`), passing
  on the room the earlier ranges of the call left, so every call returns NCBI's pairs and scan range
  (test `x_scan_subject_bf_equals_ncbi_scan_on_random_subjects` against `s_blast_aa_scan_subject`:
  TBLASTX tables, subjects of one or many ranges, offset arrays down to the longest chain). The
  offset arrays get three entries more while the switch is on. Only a table whose presence bits
  are set in fewer than half of its backbone cells is scanned this way (`x_scan_bf_for`, decided once
  per unit group and once per reference batch): the table of a translated query of more than about
  60 kb holds every word of the 25 common letters, its bit is set at nearly every subject word, the
  branch is predicted, and the list only adds work (0.4-0.8 % slower on those tables). Timed on the
  workstation (perf-si-20261012 session SCANBF-1, node 1, one binary, without and with the variable
  in turn, one output hash per case; time with / without): the first L kb of MelaMJNV against
  PemoMJNVA at one thread 0.72-0.73 at 2 and 5 kb (presence density 0.16-0.29), 0.81 at 10 kb
  (0.37), 0.91 at 20 kb (0.44), 0.98 at 40 kb (0.49); from 80 kb (0.51; MelaMJNV x PemoMJNVA,
  AP027131 x AP027133) NCBI's loop runs in both and the time is the same (0.999-1.004; before the
  density rule 1.004-1.007); 3 fruit fly transcripts (2.0 kb) against the fruit fly genome database
  0.84 at one thread (9.74 -> 8.22 s, NCBI BLAST+ 2.17.0 `-db` 9.56 s; branch misses 531 M -> 361 M,
  cycles 30.3 G -> 25.6 G) and 0.87 at 8 (2.29 -> 2.00 s, NCBI 2.48 s); 10 fruit fly transcripts
  (18.8 kb) against the zebrafish cDNA database 0.85 at one thread (59.95 -> 50.81 s) and at 8
  (7.83 -> 6.68 s, NCBI 5.88 s). With `LOSAT_X_PAIRPARGATE` these 10 transcripts run in the
  reference frame loop, which scans each subject frame once with the batch's table: 0.99 (47.33 ->
  46.91 s). Joined the strict set on these timings (no case slower; Owner decision of 2026-10-17 that
  the branch-free scans of blastp / blastx / tblastx join the set once each shows no loss).
- `LOSAT_X_TBNUNITXLATE=1` tblastn (speed and memory; in the strict set from the pull request that added it,
  2026-10-23; `=0` turns it off; `tblastn/x_subject_translate.rs` `x_packed_lean_plan`, `x_unit_residues`,
  `tblastn/search_gapped.rs` `XUnitResiduesPool`): the plan of a packed subject (every `-db` sequence, `-subject`
  with `TBNSTORE`) whose first frame has 1 MiB or more (3,145,722 bases) and which is made for its use alone
  (a search that keeps no plan of it: fewer than 4 query batches, or, since 2026-10-28, a plan that the 256 MiB
  of kept plans of `TBNSSIDE` refuses and every batch makes per use; `XSubjectSide::keeps`) has no residues: its
  frames have their lengths and its units
  their chunks, and each (frame, chunk) unit translates the residues it reads, from the chunk's first to 4,096
  past its last (or the frame's end sentinel), into a buffer the threads share, from the packed bytes of its
  phase and strand with the tables and loops of `TBNXLATESIMD` (32 codons with AVX2, else eight from one word).
  Before, the thread that took a chromosome translated its six frames (up to 83 MB each on the human genome)
  before the chromosome's units started, glibc mapped them fresh and unmapped them after the subject
  (`x_frame_capacity`), and the frames of the chromosomes searched at the same time made the peak memory. A unit
  reads its chunk and nothing else and its residues are the frame's bytes, so the reports are the same (unit
  tests `x_one_strand_of_a_phase_is_that_strand_of_both`, `x_frame_residues_of_a_range_are_the_frames`,
  `x_lean_plans_are_the_plans_without_residues`, `x_units_of_a_frame_of_three_chunks_read_the_frames_residues`;
  `long_packed_subject_chunk_wordfinder_hsps_match_ncbi`, NCBI's initial HSPs of a subject of two chunks per
  frame; `tests/tblastn_unit_xlate.rs`). The unit takes an idle buffer of its query batch's pool
  (`search_gapped.rs` `XUnitResiduesPool`, part of the batch's query-side state with its lookup table), or a new
  one, and gives it back after it, so a unit the thread runs while it waits (`TBNFLOWHELP`) takes another; a
  buffer only grows (up to a chunk and the margin, about 5 MB; from 1 MiB on with a capacity glibc maps,
  `x_reserve_mapped`), and the batch's buffers are freed at the end of its preliminary stage. Before 2026-10-28
  (B11m) each thread kept its buffer until it ended, so the buffers of every thread that had run such a unit
  were still there at the peak, which comes after the units. Measured 2026-10-29 on the workstation (node 1,
  `--membind=1`, the pull request's binary against dev ddb8c7ad in the same run, 8 rounds with the order rotated, VmHWM
  sampled every 2 ms; every output of a case the same bytes across arms): the fly genome with `dmel_b4.faa` (4
  batches) at 16 threads max RSS 454 -> 381 MiB (-66 MiB per round, 8 of 8 lower; dev against itself -4), one batch
  of 35 human proteins against the human genome at 16 threads 2,494 -> 2,441 MiB (-53, 7 of 8; dev against itself
  -21), the zebrafish genome with `drer_b4.faa` (4 batches) at 32 threads 1,359 -> 1,319 MiB (-45, 8 of 8), the human
  genome with `hsap_b4.faa` (4 batches) at 32 threads 4,743 -> 4,751 MiB (-40 per round, 5 of 8: no change beyond the
  rounds' spread of -227 to +227); time x0.991 / 0.999 / 0.999 / 0.998 (dev against itself x0.995 / 1.001), minor
  faults +5 % on the 4-batch human search (its buffers made once per batch). The units' codons count in
  `[X_STATS] TBNXLATESIMD`. Subjects with shorter frames keep the frames of their plans (those of 256 KiB or less reused
  by `TBNFRAMEREUSE`). With `LOSAT_X_STATS=1` `[X_STATS] TBNUNITXLATE plans= units= residues= buffers= freed=`
  counts the plans without residues, the units' residues, and the buffers made and freed (every buffer made is
  freed at the end of its batch's preliminary stage). Measured 2026-10-23/24 on the workstation (physical cores
  of node 1 with their HT siblings idle, `--membind=1`, cores warmed 3 s before each run, pairs with the order
  alternated, `perf stat` walls; node 0 busy, load 13-54), the pull request's binary against its base, one protein
  against a genome (`-db -evalue 1e-5 -max_target_seqs 5`), every output of a case the same bytes across arms (the
  human genome's equal to NCBI BLAST+ 2.17.0 `-db -num_threads 1`): the human genome at 1 / 8 / 16 / 32 threads
  15.71 -> 13.33 s (x0.848), 2.82 -> 2.41 (x0.854), 2.10 -> 1.64 (x0.783), 2.06 -> 1.36 (x0.663), max RSS 1,264 /
  3,725 / 5,349 / 6,685 -> 849 / 916 / 985 / 1,094 MiB (NCBI 1,089 / 3,657 / 5,529 MiB at 1 / 8 / 16), page faults
  1.71 M -> 0.20-0.27 M, instructions +2.3 %, cycles -2.6 %; the zebrafish genome x0.881 / 0.865 / 0.798 / 0.663, max
  RSS 517 / 1,295 / 1,919 / 2,942 -> 427 / 494 / 562 / 679 MiB; the fly genome (one batch of proteins) x0.992 at one
  thread, x0.994 at 8, x0.976-1.001 at 16 (three runs), max RSS 324 / 352 / 398 -> 155 / 227 / 356 MiB at 8 / 16 / 32, and at 32
  threads x1.011 (median of 28 pairs over five runs whose medians were x0.983-1.020; 11 of 28 faster; the same
  cycles and user time, the flow's tail: 21-22 of 32 CPUs busy in both). E. coli x 16 protein batches (the plans of
  `TBNSSIDE`) and fly proteins against 5,000 16S sequences (frames below 1 MiB) take the frames as before
  (x1.000-1.006, the same instructions). With `=0` the binary ran as its base (the human genome at one thread
  x1.002, the same instructions). Widened 2026-10-28 to the plans a search of 4+ batches does not keep: whether
  a subject's plan is kept is decided from its length before the plan is made (`XSubjectSide::keeps`, the
  same bytes and limit as before, `plan_bytes_of`), and a plan the limit refuses, which every batch made per
  use with its six frames (a human chromosome of 128 Mb or more never fits; a search of 4 batches at 16 or 32
  threads keeps its side, `TBNKEEPFEW`), is made as a plan without a side is: without residues when its
  first frame has 1 MiB or more. Same reports; no new rule (`=0` makes every plan with its frames, as before
  the widening for such a search). Measured 2026-10-28 on the workstation (node 1, physical cores with their HT
  siblings idle, `--membind=1`, cores warmed 3 s, pairs with the order alternated, walls by the driver's monotonic
  clock, peak anonymous memory sampled every 50 ms from the process; load 22-138 from node 0), the pull request's
  binary (82db3a82) against its base (dev 2087618e), `-db -evalue 1e-5 -max_target_seqs 5`, every report the same
  bytes in every arm: 122 human proteins (4 batches) against the human genome at 16 threads 137.2 -> 135.7 s
  (x0.990, 2 of 2 pairs faster), peak RssAnon 7,267 -> 3,477 MiB (-3,791 MiB; -3,856 and -3,725 per pair), page
  faults 8.30 M -> 2.65 M; at 32 threads 96.8 -> 95.3 s (x0.984, 2 of 2), 9,483 -> 4,700 MiB (-4,783), faults
  8.73 M -> 3.00 M, where the same binary with `=0` ran as the base (x1.000, -140 MiB); 4 batches of zebrafish
  proteins against the zebrafish genome at 32 threads (23 of 25 chromosomes refused) 26.7 -> 25.8 s (x0.965, 2 of
  2), 3,318 -> 1,327 MiB (-1,991); one batch (35 human proteins, no side) x1.001 and -64 MiB (the unchanged path).
  **Where it does not hold:** 4 batches of fly proteins against the fly genome at 16 threads (one refused
  chromosome of about 32 Mb whose six frames were not at the peak) x0.996 but peak RssAnon 384 -> 425 MiB (+39;
  +39, +16 and +49 per pair): the per-thread unit buffers of its lean plan (a chunk of its frame, about 3.6 MB per
  thread that ran one of its units) stay for the search where the refused frames were freed after the subject.
  NCBI BLAST+ 2.17.0 on the first batch alone (35 proteins) at 16 threads holds 6,240 MiB of anonymous memory
  (capped at 900 s, unfinished; node 0), above the 4-batch search with the change (3,477) and the base's one batch
  (2,474).
- `LOSAT_X_TBNSTORE=1` tblastn `-subject`: each subject record's letters are packed as soon as the
  reader has read it (`sequence/packed_subject.rs`, the shape of a BLAST database volume: ncbi2na 4
  bases per byte with the remainder count in the last byte, the runs of letters other than `ACGT`
  as ncbi4na codes, the lower-case intervals), so the letters of the whole set are never held.
  With `-subject_loc` a record is first cut in place to its interval, into the record
  `seq_range::cut_subjects` makes, whose length seeds the resolution of ambiguous bases as NCBI's
  subject set-up does; the records keep their deflines. With two or more threads a second
  thread packs the records while the reader reads the next ones (`x_packed_while_read`, since
  2026-10-25 in `blastinput/x_packed_subjects.rs`, shared with `LOSAT_X_TBXSTORE`: batches of
  up to 1,024 records or 4 MiB of letters, handed over one at a time), so the packing does not add
  to the read; the records and their stores are those of the one-thread read, in the same order.
  A WebAssembly build packs on the reader at every thread count (since 2026-10-18): there the
  second thread was one host worker more than `-num_threads` asked for (`check_wasm_threading.py`
  saw two `spawn_attempt` lines at `-num_threads 2`). The unit path of the preliminary
  stage translates a subject's six frames straight from its packed bases
  (`x_packed_plan` in `tblastn/x_subject_translate.rs`: the codon tables of `TBNXLATE`, graduated, the codons
  read from the packed bytes, each phase filling one plus frame forward and one minus frame
  backward; see `TBNXLATESIMD`), as NCBI's
  `BLAST_GetAllTranslations` translates the compressed subject, so a thread holds the frames of
  the subject in hand and no letters, resolved bases or reverse complement; the letter-wise
  reference translation went with `TBNXLATE` (2026-10-22), and the serial preliminary stage that
  read the packed bases as letters is a `#[cfg(test)]` harness since `TBNPAR` (2026-10-22); with
  `-lcase_masking` the stage reads the letters. The redo (composition mode 2), its translations,
  the traceback of composition mode 0 (also every BLOSUM45 run) and the pairwise rows (outfmt 0)
  read windows of the subject (`SubjectNt` in `packed_subject.rs`: the HSP's codons, or the window
  around an HSP that NCBI translates), made in upper case from the store when read, exactly as the
  reader gave them but for the case (every reader ignores the case; the comparison tests in
  `packed_subject.rs`), so no stage after the preliminary one makes a subject's whole letters again
  but the traceback of a list retried after a fence, whose frames NCBI translates whole (one window
  of the whole subject for all of them). Every `-db` sequence is held packed too, with no switch.
  The packed bases are those NCBI's subject set-up makes
  (`encode_subject_ncbi2na_packed`, `s_Ncbi4naToNcbi2na` in `blast_objmgr_tools.cpp:427-475`; one
  pass over the letters). Memory: about a quarter of a byte per base held instead of one; per
  thread the frames of the subject in hand (two bytes per base) instead of four bytes per base
  (letters or resolved bases, reverse complement, frames). A set with a byte that has no ncbi4na
  code (only an API caller can make one) keeps the letters; the library entry `run_local` reads the
  caller's records. The plans of letters are made from them packed (below).
  With `-outfmt 11` the reading takes the `=0` side (the records keep their letters, one byte per
  base, and with `-subject_loc` the whole subjects), because every archive repeats the whole
  subjects; the plans still pack per plan (+1 byte per subject base: 33 against 19 MiB on 20,000
  subjects, 18.8 Mb; NCBI 282 MiB).
  In the set since 2026-10-16 (on unless `0`). With
  `LOSAT_X_STATS=1` a `[X_STATS] TBNSTORE
  subjects= stored_bytes= letters_made= packed_plans=` line gives the subjects held
  packed, the bytes of their stores, the times their whole letters (outside the windows above)
  were made, and the plans translated from the packed bases; `LOSAT/tests/tblastn_subject_store.rs`
  reads it. The records `run_local` (LOSAT Web) is given keep their letters; since 2026-10-23 the
  worker that makes a subject's plan packs its letters (`x_plan_store` in
  `tblastn/x_subject_store.rs`) and makes the plan from the packed bases, a quarter of a byte per
  base in place of the plan's ncbi2na copy of one byte per base, while the plan is made (LOSAT
  Web, one protein against the fly genome: the run x0.85-0.89 at 1, 4 and 8 threads, linear
  memory at 8 threads 510 -> 452 MiB); a subject with lower-case letters searched with
  `-lcase_masking` keeps its letters. `packed_plans` counts these plans too.
- `LOSAT_X_BPGAP=1` blastp preliminary stage: the restricted gapped alignment (NCBI
  `s_RestrictedGappedAlign`, `blast_gapalign.c:991-1300`, which `s_BlastProtGappedAlignment` runs
  for a query whose best ungapped score is below 0.68 of the cutoff) keeps the two score rows and
  the column gaps of the columns that are multiples of 10 in separate arrays, and on x86_64 with
  AVX2 computes eight cells of a row at a time (`blastp/x_restricted_gap.rs`, exactness argument in
  its header): in the rows where no gap opens in A the cells depend on each other only through the
  X-drop test, whose `best_score` is a prefix maximum; in the other rows the row gap is a decayed
  prefix maximum. 32-bit values as in the reference, the matrix rows read from 16-bit tables built
  from the reference's score function. Without AVX2 (other CPUs, aarch64, Wasm) the reference
  function runs (the scalar loops on the same arrays were slower than it). A letter
  code of 28 or more, a negative gap cost, a gap cost above 65,536 or an X-drop above 2^28 leaves
  the call to the reference function. In the strict set since 2026-10-14 (`=0` turns it off). With
  `LOSAT_X_STATS=1` a `[X_STATS] BPGAP fast= not_taken=` line on stderr at the end.
- `LOSAT_X_BPGAP16=1` (BLASTP, speed; in the strict set from the pull request that added it,
  2026-10-23; `=0` turns it off; acts only where `LOSAT_X_BPGAP` runs its vector alignment) the
  restricted gapped alignment's rows in 16-bit lanes, sixteen cells per step instead of eight
  (`restricted_rows16_avx2` in `blastp/x_restricted_gap.rs`, exactness argument in its header,
  "Sixteen-bit rows"): cells, column gaps and the row gap in 16 bits with saturating sums,
  `i16::MIN` where the 32-bit rows hold `MININT`. Taken where the X-drop and `gap_open + gap_extend`
  are at most 4,096, `gap_extend` at most 256 and the largest matrix entry at most 4,096 (every
  matrix and gap cost pair blastp supports); a row starts in 16 bits only while the best score is at
  most 32,767 - largest entry - 16 x `gap_extend` - 1 (32,739 for BLOSUM62 with 11/1), and past that
  the 32-bit rows compute the alignment again from row 0. x86_64 with AVX2, BMI1 and LZCNT only:
  aarch64 and WebAssembly keep the reference function (no AVX-512 on the default path).
  Measured 2026-10-23 on the workstation (node 1, CPU 40 for one thread and CPUs 40-47 for eight,
  their HT siblings idle, warmed 3 s before every run, arm order reversed every other round, medians
  of the per-round ratios against dev 7cfec941, the binary before the switch; every report identical
  between the arms and equal to NCBI BLAST+ 2.17.0's). OrthoFinder jobs of 7 species at one thread
  (`-subject`, `-evalue 0.001`, 4 rounds): Blast0_6 x0.959 (54.56 -> 52.35 s, 4 of 4 rounds below 1;
  cycles x0.959, instructions x0.974), Blast6_6 x0.957 (18.32 -> 17.55 s). With the prototype of the
  same rows (8 rounds at one thread, 2 at eight): Blast0_6 x0.957, Blast6_6 x0.956, 5,000 x 5,000
  UniRef50 proteins (`-db`, `-evalue 0.001`, 8 threads) x0.962, a 42-letter UniRef50 query against
  the first UniRef50 volume (`-evalue 1e-5`, one thread, few restricted alignments) x1.001. In
  Blast0_6 the restricted rows went from 20.1 to 13.3 G cycles (11.6 % of dev's cycles in `perf
  record`). With `=0` the instruction counts equal dev's (x1.00006) at x1.002-1.003 of its time. Max
  RSS within the spread of repeated runs of one binary.
- `LOSAT_X_BPGAPBAND=1` (BLASTP, memory; in the strict set from the pull request that added it (perf/bp-gap-band,
  v4-D12 s23-s25, `DECISIONS` 66-70): on unless it is `0`; acts only where `LOSAT_X_BPGAP` runs its vector rows) the
  restricted gapped alignment's arrays (`XRestrictedScratch` in `blastp/x_restricted_gap.rs`: the score rows `h`, the
  column gaps `e10`, the column letters `col`, and the 16-bit rows' `h16` and `e16`) hold the cells the band reaches
  instead of the subject letters left or right of the seed (`len2`). NCBI starts `dp_mem` at `num_extra_cells + 100`
  cells, `num_extra_cells = x_dropoff / gap_extend + 3` (`blast_gapalign.c:1040-1051`), and in the row loop grows it to
  `MAX(last_b_index + num_extra_cells + 100, 2 * dp_mem_alloc)` (`1261-1269`). The scratch sized and filled its
  arrays for `len2` and never shrank them, so each thread kept arrays for the longest subject it had aligned, about
  15.4 bytes per letter (UniRef50 q100 at 32 threads: 16.3 MiB held, 26.6 MiB of capacity; sized by the furthest cell
  each scratch's band reached, 0.7 MiB). With the switch `restricted_fast` starts the arrays at `ahead + 100` cells,
  `ahead` = `(x_dropoff + 2 * max) / gap_extend + 3 + ROW_AHEAD` (`max` the largest matrix entry: a row's band
  extension starts from a score at most `best_score + 2 * max`; `ROW_AHEAD` = 47 cells the rows read past `b_size`;
  `len2` without a gap extension cost, as NCBI's `N + 3`), and gives `col` `held - ahead` letters (`band_col`):
  `fill_col` returns false past them, so the rows stop (None) before a row's band passes the cells held; the arrays
  then double (NCBI's `2 * dp_mem_alloc`, up to `len2 + ahead + 2 + TAIL` cells, where the rows copy every letter) and
  the rows run again from row 0, as the 16-bit rows do when they hand over to the 32-bit rows (3,000 human proteins x
  `BlastDBSpecies6`: 77 stops in 1.54 M alignments). `col`'s length depends on the arrays only; with `=0`, or without a
  gap extension cost, the arrays hold `len2` cells and `col` only grows, as before. The row functions are those of
  the build before the switch: growing the arrays inside the row loop (a comparison per row, or in the letter-copy
  block) changed the register allocation of the 16-bit rows, +2.0 to +3.0 % of their instructions with the switch on
  or off; this form reads x1.00029 of the dev base's instructions there and x1.00030 in OrthoFinder Blast0_6
  (callgrind, one thread). The cells the rows read and write, and every value, are those of the arrays of `len2`
  (lib tests `x_bpgapband_reused_scratch_matches_reference`, `x_bpgapband_arrays_hold_the_band_not_the_subject`, and
  every restricted-alignment test of `check()` also with the switch on, on a scratch of its own;
  `tests/blastp_restricted_gapped.rs` `restricted_gapped_alignments_match_ncbi_with_and_without_the_switch` compares
  NCBI BLAST+ 2.17.0's rows with `=0` too).
  Evidence 2026-10-29/30 (v4-D12 s24-s25, workstation, head d78bb8a8 on dev 593fff35; arms default, `=0`, and the
  dev base built the same way, 8 rounds rotated, peak RssAnon from 2 ms samples; every report identical across
  arms): UniRef50 q100 `-evalue 0.001` at 32 threads (node 1, membind 1) 141.2 [139.3-143.5] MiB against the dev
  base 161.2 [159.0-162.7] and `=0` 161.6 [156.8-162.4] (17.9-22.0 MiB lower in 8 of 8 rounds; wall x1.0027
  [0.9955-1.0138] of the dev base, `=0` x1.0000 [0.9969-1.0042], in runs not warmed beside the sampler); q100
  against a third of UniRef50 at 16 threads (two runs side by side, 12 slots) 67.8 [59.0-86.5] against 69.6
  [63.2-77.8] and `=0` 66.2 [63.3-72.2] (1.0-5.6 MiB lower than the dev base in the same slot, 4 of 4; this
  case's peak spreads by +-16 MiB between slots, the arrays are 0.1 MiB of it); 3,000 human x `BlastDBSpecies6`
  at 16 threads 37.4 [37.1-37.5] against 38.2 [38.1-38.6] and `=0` 38.3 (8 of 8 lower); one thread (node 0)
  human 12.9 against 13.0, OrthoFinder Blast0_6 19.2 against 19.2 MiB. Time: one thread (`perf stat` `:u`, node
  0, CPUs 18-19) instructions x1.0003 of the dev base in both cases (callgrind: x1.00029 and x1.00030, the
  16-bit rows +0.02 %), cycles human x0.9932 (7 of 8 rounds below 1; `=0` x0.9978), Blast0_6 x0.9405 with `=0`
  x0.9423, the binary's placement (both arms of the branch 139.1 G cycles, the dev base 146.5-148.3 G in three
  runs), not a gain of the switch; 16 threads (node 1 CPUs 32-47 warmed 3 s, 8 rounds) Swiss-Prot q100
  instructions x1.0000, cycles x1.0006 (`=0` x1.0012), wall x1.0008 (`=0` x1.0000); human instructions x1.0003,
  cycles x1.0012 (`=0` x1.0031), wall x1.0015 (`=0` x0.9985).
- PAIRPAR, TBLASTX half (graduated 2026-10-22; no variable, no `LOSAT_X_PAIRPARSHADOW`; class
  conditional, the name stays in `GRADUATED` for the reach table; it was in the strict set; the
  blastn half is above) TBLASTX seed stage (scan, two-hit test, ungapped extension) of each subject
  as units of (subject frame, chunk, query context), 36 per one-query pair, on the search pool (in
  order on the calling thread with one thread); each unit scans with a lookup table holding only its
  context's chains (order kept) on a fresh unmasked diagonal table, and the unit lists go through the
  reference's sort and merge in the reference order (`tblastx/x_pair_par.rs`, exactness argument in
  its header). The rule stays, per batch, in this order (`decide`): the reference frame loop runs
  ("fallback: gate") with the trace, debug, diagnostics and chunk switches
  (`LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS`, `LOSAT_TBLASTX_PARALLEL_CHUNKS` at 2+ threads, threaded Wasm
  with a subject of more than one chunk); then ("fallback: guard") when `slack` (diagonal table size
  minus query length) plus twice the longest subject chunk reaches 2^29 (only with `-window_size`
  near 2^29) or, without `LOSAT_X_TBXBAND` (`=0`), no query
  context and no subject chunk is at most `slack` long; then ("fallback: cost") by
  `LOSAT_X_PAIRPARGATE` (below). Otherwise the units run. With `LOSAT_X_PAIRPARGATE=0` every batch the gates and the guard let through takes the units.
  `LOSAT_X_STATS=1` counts `fallback_guard` / `fallback_gate` / `fallback_cost` and `units` in
  `[X_STATS] PAIRPAR`, which is how the tests tell the frame loop from the units; `LOSAT_TIMING`
  prints one `x_pairpar` line. The shadow (`LOSAT_X_PAIRPARSHADOW=1`: the unit path and the reference
  loop for every subject, the subject's combined HSP list compared field by field, one `[X_SHADOW]
  PAIRPAR` summary line) went with the variable; the unit tests of `x_pair_par.rs` still compare the
  per-unit hits with `masked_reference`, and unit path against frame loop is checked only through the
  frozen NCBI bytes of the fast fixtures (at one thread the frame loop for batches under 5,000 letters,
  `LOSAT_X_PAIRPARN1`, the units above; the units at 4 threads). The tests that need the frame loop reach it at one thread
  with `LOSAT_X_PAIRPARGATE=1` (the cost rule; batches under 5,000 letters, or any below 2^18 cells with
  `LOSAT_X_PAIRPARN1=0`), with several subjects without a dominant one and `LOSAT_X_TBXFEWGROUPS=0` (or a
  batch under 5,000 letters), or with
  `LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS=1`, and assert `units=0`.
- `LOSAT_X_TBXBAND=1` (TBLASTX, speed; in the strict set from the pull request that added it,
  2026-10-19; `=0` turns it off; acts where the unit path runs) the unit path no longer applies
  its first guard (every query context, or every subject chunk, at most `slack` long). That guard
  came from a looser form of item 3 of the exactness argument in `tblastx/x_pair_par.rs`; the
  argument holds without it: the last right extension of a raw diagonal, started at query offset
  q0, ends at the end of q0's context, at most `query_length - q0` letters right of its seed, while
  a hit of a diagonal that shares its masked cell lies at least `table size - q0` letters right of
  it, so that hit is never masked. The batches the guard sent to the reference frame loop, which
  searches a subject on one thread, are those whose longest context and longest subject chunk are
  both longer than `slack`: for one query of N bases, N in (3/7 2^k, 2^(k-1)] for the table of 2^k
  cells (1.76-2.0 kb, 7.0-8.2 kb, ..., 56-65 kb, 112-131 kb, 225-262 kb, 0.45-0.52 Mb, 0.90-1.05 Mb,
  1.8-2.1 Mb, 3.6-4.2 Mb), and multi-query batches whose letters lie just below a power of two,
  each against a subject frame longer than `slack`. Those whose longest subject chunk has at most
  as many letters as the reference table has cells (subjects of up to about six times the query's
  bases) now go through `LOSAT_X_PAIRPARGATE` like any other batch (batches under 5,000 letters keep
  the reference loop at one thread, `LOSAT_X_PAIRPARN1`). The others keep the reference loop without `LOSAT_X_TBXUNITMASK` (with it,
  masked unit tables; next entry): a memory rule (`x_band_fits`), since a
  unit's diagonal table has a cell for every raw diagonal of its query group and subject chunk, one
  table per thread, against the reference's one table of 2^k cells per subject in flight.
  `LOSAT_X_STATS=1` adds `band=` (such batches the units took) to the `[X_STATS] PAIRPAR` line.
  `tests/tblastx_band.rs` runs a
  band batch through the binary. Measured (workstation node 1, 16 or 8 booked cores or one, each
  run after a 3 s warm-up of its cores, alternating order, medians of 3 runs unless noted, every
  report equal across arms), against `=0`: Sakai 1 Mb x MG1655 1 Mb at 1 thread 27.27 -> 18.59 s
  (x0.682), 8 threads 24.85 -> 3.28 s (x0.132), 16 threads 24.69 -> 2.42 s (x0.098), max RSS 184 ->
  129 / 148 / 167 MiB; Sakai 1-3.9 Mb x MG1655 at 1 thread 668 -> 325 s (x0.49, one run each),
  8 threads 616 -> 63.3 s (x0.10, 1 and 2 runs), 16 threads 624 -> 52.2 s (x0.084, 2 runs), max
  RSS 1,035 -> 700 MiB at 1 thread and 1,389 -> 1,348 MiB at 8 and 16. Outside the band, Sakai
  0.88 Mb x MG1655 1 Mb against the binary before the switch x0.965-0.977 (the same path; build
  differences), max RSS equal. Without the rule, 3 fruit fly transcripts (2,005 bases) against the
  fruit fly genome database ran 1.91 -> 1.41 s at 8 threads and 1.93 -> 0.87 s at 16 but its max
  RSS grew from 318 to 600-650 and 873 MiB (NCBI BLAST+ 2.17.0: 368 MiB at 8 threads), and a 1.9 kb
  query against MG1655 0.30 -> 0.10 s and 0.31 -> 0.09 s with 21 -> 68-82 and 111 MiB; without
  `LOSAT_X_TBXUNITMASK` the rule leaves both on the reference loop (as before); with it (the default)
  they take the units on masked tables (next entry).
- `LOSAT_X_TBXUNITMASK=256` (TBLASTX, speed and memory; a value switch, in the strict set from the
  pull request that added it, 2026-10-21; `=0` turns it off, `1` or unset is the built-in 256, `=N`
  N KiB; acts where the unit path runs) a unit of the unit path takes a diagonal table of P
  cells, P the smallest power of two of at least its query group's extent plus the window, and finds
  a hit's cell as (q - s - lo) & (P - 1), instead of a table with a cell for every raw diagonal of its
  group and subject chunk: where the unmasked one has at least four times P cells and P is at most
  the switch's KiB at eight threads or more, at most 4,096 cells at two to seven, and, in a batch
  whose first guard fails
  (`LOSAT_X_TBXBAND`), wherever the subject chunk is longer than the reference table has cells, so
  that `LOSAT_X_TBXBAND`'s memory rule (`x_band_fits`) no longer keeps such batches on the reference
  frame loop. Exact by item 3's masked-aliasing argument applied to the group (`tblastx/x_pair_par.rs`
  module doc): two raw diagonals that share a masked cell lie at least P minus the group's extent,
  at least the window, apart along the subject, so the first hit of the later one leaves the cell as
  it leaves a fresh cell. A batch whose groups `LOSAT_X_TBXQSLICE`'s rule cuts keeps unmasked tables
  (the scan in slices brings a cell's diagonals out of subject order). A masked table costs one
  instruction per hit and saves the clearing of the unmasked table's every cell per unit, hence the
  size rule. `LOSAT_X_STATS=1` prints `[X_STATS] TBXUNITMASK kib= units_masked=` (only after a TBLASTX batch; the
  `[X_SHADOW] TBXUNITMASK` line that `LOSAT_X_PAIRPARSHADOW=1` printed went with that variable).
  `tests/tblastx_band.rs` runs a band batch against a genome through the binary. Measured
  (workstation node 1, 16, 8, 4 or 2 booked cores or one, each run after a 3 s warm-up of its cores,
  order rotated, medians of 8 runs, every report equal across arms and to NCBI BLAST+ 2.17.0 where
  it ran), against `=0`: band batches against longer subjects, which ran the reference frame loop:
  Sakai 1-500,000 x MG1655 at 16 threads 57.27 -> 5.34 s (x0.093; max RSS 217 -> 226 MiB, NCBI 315 MiB
  in 3,433 s), 3 fruit fly transcripts (2,005 bases) against the fruit fly genome database x0.351
  at 16 threads and x0.549 at 8 (311 MiB, NCBI 360), Sakai 100,001-101,900 x MG1655 x0.273 / x0.297 /
  x0.531 at 16 / 8 / 2 (21 -> 24 MiB, NCBI 91); small groups against long subjects, which took
  unmasked tables: Sakai 100,001-101,500 against the fruit fly genome database x0.853 / x0.865 at 16
  / 8 threads with max RSS 784 / 570 -> 309 MiB (NCBI 356 MiB: above NCBI before), against MG1655
  x0.800 (110 -> 24 MiB), Sakai 1-100,000 x MG1655 x0.987 at 16 threads (143 -> 61 MiB, NCBI 130),
  x1.000 at 8, Sakai 1-50,000 x MG1655 x1.000 (125 -> 40 MiB); at four and two threads Sakai
  100,001-103,000 and 1-10,000 x MG1655 (tables of 2 and 4 Ki cells) x0.941 / x0.923 and x0.971 /
  x0.957, Sakai 100,001-101,500 x0.964 / x0.895 (50 / 41 -> 24-25 MiB), while 1-50,000 and 1-100,000
  (32 and 64 Ki cells, left unmasked there) lose 0.2-0.9 % masked; Sakai 1-30,000 x MG1655, a band
  batch, x0.256 / x0.447. The groups the rule
  leaves unmasked (Sakai 1-300,000 and 1 Mb x MG1655) ran as before (x1.000-1.001; against the binary
  before the switch x0.992-1.004, instructions +0.06-0.25 %). With `=512` the groups of 100,000
  letters (Sakai 1-300,000 x MG1655) ran x1.008-1.012 (165 -> 143 MiB; NCBI 224); at one thread
  masking only adds the instruction (x1.006 on Sakai 1-100,000). Five to seven threads were not
  measured. Prefetching the masked tables (`LOSAT_X_TBXPREFETCH=2`) was slower (x1.008-1.028).
- `LOSAT_X_PAIRPARGATE=1` (TBLASTX, speed; in the strict set since 2026-10-18, on by default, with
  `LOSAT_X_TBXLUT`; added off unless set in #255; since PAIRPAR graduated on 2026-10-22 it is the
  only reader of the cost rule) a query batch that the
  gates and the guard of PAIRPAR let through takes the unit path only where it is expected to be faster than the
  reference frame loop (`units_gain` in `tblastx/x_pair_par.rs`): the reference's diagonal table has
  2^18 cells or more (query length plus window above 131,072 letters, a translated query batch of
  more than about 65 kb), or the search has two or more threads and one subject holds more than one
  and a half thread's share of the bases (`2 * longest * threads > 3 * total`), or two or more threads
  and a batch of 5,000 letters or more (`LOSAT_X_TBXFEWGROUPS`, the line of `LOSAT_X_PAIRPARN1`), or the lookup
  threshold is below TBLASTX's default of 13 (not measured; a lower threshold puts more words in
  every position, and the full table grows faster than the group tables: with `-threshold +inf`,
  which NCBI converts to INT_MIN, every word, the option-sweep inputs took 14.2 GB in the reference
  loop, 3.3 GB in the units and 5.9 GB in NCBI). Otherwise the batch
  runs the reference loop with the full lookup table, as for the guard (`LOSAT_X_STATS=1` counts it
  as `fallback_cost` on the `[X_STATS] PAIRPAR` line). A unit scans its subject frame
  with one group's table, so each frame is scanned once per group (up to 6 times), and every unit
  makes its own diagonal table and lists: on short queries this costs more than the small tables
  save. Both paths give the reference's list (the exactness argument of `tblastx/x_pair_par.rs`), so
  the report does not change; `LOSAT/tests/tblastx_unit_gate.rs` compares reports and reads the
  counters. Measured on the workstation (alternating pairs, every report equal across arms;
  `LOSAT_X_TBXNEXT` and TBXFRAMES on in both arms): 10 fruit fly transcripts (two query
  batches of 10.8 and 7.9 kb) against the zebrafish cDNA database (57,775 sequences) 61.3 -> 46.6 s
  at one thread (NCBI BLAST+ 2.17.0 `-db` 45.1 s), 15.75 -> 11.87 s at 4 (NCBI 11.70), 8.01 -> 6.00 s
  at 8 (NCBI 6.05), 4.07 -> 3.08 s at 16 (NCBI 3.22), max RSS unchanged (45 MiB); the first 40 kb of
  MelaMJNV against 3,000 zebrafish cDNAs 4.14 -> 3.51 s at one thread. Unchanged (the units kept):
  queries of 80 kb and more (MelaMJNV x PemoMJNVA, AP027131 x AP027133, where the units take
  0.74-0.90 of the reference's time and half its memory), and a 5 kb query against one 294 kb
  genome at 4 and 8 threads. Before the gate, the units' time over the reference's at one thread was
  1.48-1.72 for 5 kb queries, 1.19-1.41 for 20 kb, 1.03-1.18 for 40 kb and 0.90-0.99 for 80 kb (one
  294 kb subject / 3,000 cDNAs). Timed again for the strict set on origin/dev 7a5a7f87 + `LOSAT_X_TBXLUT`
  (workstation, 16 cores of node 1 with their HT siblings idle, a warm-up then 2-3 rounds with the
  order reversed in every second round; every report equal across arms, and equal to NCBI where NCBI
  ran), the gate and `LOSAT_X_TBXLUT` both on against both off: the 10 transcripts against the
  zebrafish cDNA database 60.48 -> 46.65 s at one thread (NCBI 45.16), 7.94 -> 6.09 s at 8 (NCBI
  5.90), 4.43 -> 3.24 s at 16 (NCBI 3.27); 40 kb x 3,000 cDNAs 4.04 -> 3.53 s at one thread; 5 kb x
  3,000 cDNAs 0.188 -> 0.122 s at 8; unchanged (x0.99-1.00): 3 transcripts against the fruit fly
  genome database at 8 threads (the guard), 80 kb x 3,000 cDNAs, MelaMJNV x PemoMJNVA at 8 threads,
  AP027131 x AP027133 (code 4) at one thread; max RSS equal or lower in every case.
- `LOSAT_X_PAIRPARN1=1` (TBLASTX, speed; in the strict set from the pull request that added it,
  2026-10-25; `=0` restores the one-thread line of 2^18 cells) at one thread `LOSAT_X_PAIRPARGATE`
  keeps the unit path for a query batch of 5,000 letters or more (`query_length`, the translated frames with their
  sentinels: a nucleotide query of about 2.5 kb; `X_GATE_N1_LETTERS` in `tblastx/x_pair_par.rs`)
  instead of only from a reference diagonal table of 2^18 cells; several threads keep the rule
  above. The 2^18 line was measured before the units' branch-free hit loop and masked tables; timed
  again, the units are the faster path at one thread for every query of 3 kb or more measured, and
  the reference loop for queries of 2.1 kb or less, where each unit's scan of its subject frame
  (once per group) outweighs the hits. Both paths give the reference's list, so the report does not
  change (`LOSAT/tests/tblastx_unit_gate.rs`). Measured on the workstation (dev 313fd4e7, one core of
  node 1 or of node 0 with its HT sibling idle, warmed 3 s before every run, a warm-up then 4-6
  rounds in alternating order, every report equal across arms; the units forced with
  `LOSAT_X_PAIRPARGATE=0` against the default, which is what the switch changes above the line):
  the first 5 / 20 / 40 kb of MelaMJNV against MG1655 (4.6 Mb) x0.881 / x0.895 / x0.820 (0.444 ->
  0.391 s, 1.260 -> 1.127 s, 2.398 -> 1.966 s); the same genome cut into 1,000 / 100 / 20 / 5 / 2 kb
  subjects x0.86-0.93; 3,000 zebrafish cDNAs x0.943 / x0.931 / x0.900 / x0.895 / x0.910 / x0.828
  (3 / 4 / 5 / 10 / 20 / 40 kb; 3.485 -> 2.886 s at 40 kb); 10 fruit fly transcripts (18,760 nt, two
  query batches) against the zebrafish cDNA database (57,775 sequences, `-db`) 45.86 -> 41.56 s
  (x0.906; the search that set the 2^18 line, 61.3 -> 46.6 s, in #255), against 3,000 cDNAs x0.906,
  against MG1655 x0.856. Below the line: 2 / 2.1 kb x0.989-1.028, 1 kb x1.15-1.28, 300 bases
  x1.94-1.96. WebAssembly (serial WASI command under Node 26, one process, arms alternating, process
  CPU): 5 / 20 kb x MG1655 x0.753 / x0.690, 5 kb x 3,000 cDNAs x0.782, 5 kb x MG1655 in 2 kb pieces
  x0.792, the 10 transcripts x 3,000 cDNAs x0.744; 1 kb x1.158 and 300 bases x2.03 (below the line).
  Max RSS +2-8 MiB (14-27 -> 17-33 MiB; NCBI BLAST+ 2.17.0 91 MiB on the 5 kb x MG1655 search, which
  it runs in 0.629 s against LOSAT's 0.444 before and 0.391 s after). At eight threads with the
  subjects spread over the pool (the rule above sends them to the reference loop) the units used
  less CPU but lost wall time on many short subjects (3,000 cDNAs x1.02-1.06), so several threads
  are left as they were.
- `LOSAT_X_TBXLUT=1` (TBLASTX, memory; in the strict set since 2026-10-18, on by default, with
  `LOSAT_X_PAIRPARGATE`; added off unless set in #262) the full lookup table of a query batch, the table of the
  reference frame loop (the batches the unit path's gate, guard or
  `LOSAT_X_PAIRPARGATE` sends back, and the batches of `LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS`), filled in place from per-cell counts by the builder of BXLUT /
  BPLUT (`x_tbx_build_ncbi_lookup` in `tblastx/lookup/backbone.rs`, `x_build_direct` in
  `x_lut_direct.rs`; same appends in the same order, so the same table field by field) instead of
  being copied out of one chain per cell. The chain builder holds every chain, and
  the arena blocks the chains left when they doubled (LUTARENA), until the table is complete (NCBI
  frees each chain as it copies it); with `-threshold +inf` (NCBI's INT_MIN: every word) all chains
  grow in step and the dead blocks are never reused. Diagnostics mode keeps the chain builder.
  `LOSAT/tests/tblastx_ref_lut.rs` compares reports with and without the switch (the reach file shows the
  in-place build); a unit test of `x_lut_direct.rs` compares the table with the chain builder's field by field.
  Measured: 3 kb of the option-sweep query against the sweep subjects at `-threshold +inf` with the
  reference loop (then run with `LOSAT_X_PAIRPAR=0`, a variable that no longer exists), one thread: max RSS 1,069 -> 475 MB, 8.20 -> 7.00 s (NCBI
  589 MB, 9.19 s); the whole sweep input 14.2 -> 3.31 GB (NCBI 5.86 GB). On the batches the gate
  sends back (the timings of `LOSAT_X_PAIRPARGATE` above) it takes x0.95-1.00 of the chain
  builder's time; max RSS 23 -> 19 MiB on 40 kb x 3,000 cDNAs, 15 -> 17 MiB on 5 kb x 3,000 cDNAs
  at 8 threads (the default before both: 21 and 19 MiB).
- TBXUNGAP (graduated 2026-10-26, in the set since 2026-10-18; no variable; always on, every target, the
  wasm32 builds (LOSAT Web) since 2026-10-19; the reference walks stay as the fallback below) the
  two-hit ungapped extension (`extend_hit_two_hit`, NCBI `s_BlastAaExtendTwoHit`: the word loop and the
  `s_BlastAaExtendLeft` / `s_BlastAaExtendRight` walks) reads its scores from one 256 x 256 byte table
  (`q << 8 | s`, the scores of `get_score` for every two bytes, `BLAST_SCORE_MIN` stored as -128) and runs
  each walk in a function of its own (`tblastx/extension/x_two_hit_fast.rs`, exactness argument in its
  header): 18 instructions a step on the left and 19 on the right, against 22 and about 30 (the reference
  calls `get_score`, which tests both codes and multiplies for the row, and its right walk reloads its
  bounds from the stack). Same counts, positions, statements and exit tests; an X-drop outside 1..=127 or a
  word that does not fit runs the reference, and so does the debug output. All callers (the unit path and
  the reference frame loops) go through it. Measured against the dev binary before it (workstation, node 1
  without other timed runs, 8 alternating rounds, medians of the per-round ratios, reports and max RSS
  identical): MelaMJNV x PemoMJNVA 1 thread x0.978, AP027131 x AP027133 (code 4) 1 thread x0.971 and 8
  threads x0.984, Sakai x MG1655[0:1 Mb] (code 11) 8 threads x0.988. perf stat at 1 thread: AP027131 x
  AP027133 instructions 72.2 -> 69.0 G, cycles 34.9 -> 34.1 G (about 100 instructions and 25 cycles fewer
  per extension). What is left of an extension is the latency of its first query loads and the
  mispredicted exits of the walks (`LOSAT_X_TBXEXTBATCH` takes the first).
  wasm32 (plan v4 line W, item W7): until 2026-10-19 the wasm32 builds kept the 16-residue v128 walks of
  `two_hit.rs`, which take each lane out of the vector and call `get_score` for it, so every step was the
  reference's scalar step. With the table, Node 26 (V8 TurboFan), warm searches of the LOSAT Web threaded
  and serial reactors, all arms in one process alternating search by search (6 rounds), CPU time per stage
  from `perf` (cpu-clock) in each search's window, divided by the stages the change does not touch
  (linking, SEG, reports, statistics, sorting) because this machine's CPU time moves by up to 2x with core
  placement. Against dev fd75dd9f (medians of 6 rounds, rounds faster in brackets): AP027131 x AP027133 (code 4)
  8 threads extension x0.772 (5/6), whole search x0.916 (5/6); serial reactor x0.807 (5/6), x0.913 (5/6);
  MelaMJNV x PemoMJNVA 8 threads x0.889 (6/6), x0.977 (4/6); serial reactor x0.846 (6/6), x0.957 (6/6);
  every report identical. (The method's noise: the same unchanged code read x0.92-1.02.)
- `LOSAT_X_TBXFEWGROUPS=1` (TBLASTX, speed and memory; in the strict set from the pull request that added it,
  2026-10-31; `=0` turns it off) the unit path takes as few context groups as the pool can use: each group
  scans every subject frame and has its own tables, and the groups serve only parallelism. A batch wants enough
  groups for two units per thread on the longest subject's share of the bases (`x_groups_wanted`: one at up
  to 3 threads, three at 8 threads on one subject, one with subjects enough), and takes them where it extends over at most
  2^17 letters (`X_FEW_GROUPS_LETTERS`); a longer batch keeps its six groups (its diagonal tables would outgrow the
  cache: 100-300 kb queries x1.03-1.05 in one or three groups). With fewer groups several threads take the
  one-thread line of the cost rule (`units_gain`, `LOSAT_X_PAIRPARN1`: 5,000 letters), so a search of many
  subjects runs the units there too. A group boundary only decides where the concatenated query is cut, so reports do
  not change. A group build's blocks of 64 KiB or more take the capacity glibc maps
  (`x_group_capacity`): batches are built on the pool's threads, and those blocks held the arenas' high water.
  Measured 2026-10-31 on the workstation against dev 305acb26 (a node 0 core with its HT sibling idle at 1 thread,
  node 1 CPUs 32-39 at 8, warmed 3 s before every run, arm order rotated, one warm-up round + 8 rounds, medians of the
  per-round ratios, time / max RSS; every report identical between the arms; with this switch come the unit path for
  other matrices and their sentinel walks). 40 x 40 phage genomes of 43 kb (`Autographivirales.40`, BLOSUM45,
  ViPTreeGen's options): 1 thread x0.581 / x0.746 (80.2 -> 46.7 s, NCBI 81.4 s), 8 threads x0.505 / x0.819 (13.4 -> 6.8 s,
  NCBI 18.0 s). The same with BLOSUM62: x0.802 / x0.970, x0.566 / x0.673. `Begomovirus.40`: x0.648 / x0.906, x0.721 / x1.201.
  MelaMJNV x 2 viruses: x0.525 / x0.581, x0.213 / x0.527. Sakai 50 kb x MG1655: x0.855 / x0.925, x0.998 / x0.898.
  10 fly transcripts x zebrafish cDNA database: x0.738 / x0.957, x0.681 / x0.957. 3 fly transcripts x fly genome database:
  x0.984 / x1.000, x0.894 / x0.989. 50 cDNAs x 1 Mb: BLOSUM62 x0.804 / x0.814, x0.966 / x0.626; the 8 other matrices x0.576-0.824
  / x0.751-1.228, x0.145-0.325 / x0.732-1.527. AP027131 x AP027133 (code 4) keeps six groups: BLOSUM62 x1.008 / x0.987,
  x1.012 / x1.005 with x0.9991 of the instructions and x1.006 of the cycles, the same with `=0` (code layout); other matrices
  x0.510-0.849 / x0.421-0.749, x0.092-0.518 / x0.445-1.008.
- `LOSAT_X_TBXUNGAPSENT=1` (TBLASTX, speed; in the strict set from the pull request that added it,
  2026-10-25; `=0` turns it off) the two-hit ungapped extension (TBXUNGAP's table walks, above) with
  walks that the sentinels stop instead of a bound (`two_hit_sent` in
  `tblastx/extension/x_two_hit_fast.rs`, exactness argument in its header): a query context's buffer
  and a subject frame's are `[0, residues.., 0]` (NCBI's NULLB at each end, blast_util.c:438-453), and
  the table scores code 0 against anything as -128, so a step on a sentinel stops either walk where the
  reference's bound stops it (the step does not raise the best and makes `best - score >= 128 >
  x_drop`), with the same best, length and stop position. The word loop and the walks (two steps a
  loop) are inlined into the callers' extension code (the unit path's `x_hit`, the batched `resolve`,
  the reference frame loops). Acts in native builds where the batch's query contexts carry both
  sentinels (checked once per batch) and the unit's subject chunk is its whole frame (frames of at
  most 5,000,000 residues, `MAX_DBSEQ_LEN`, so subjects up to 15 Mb; a longer frame is searched in
  chunks whose ends are residues) of at most `i32::MAX / 64` bytes (no score wraps at +63, the largest
  score of a table). A matrix other than BLOSUM62 takes the same walks on its byte table
  (`two_hit_sent_in`, `TbxScore::sent_two_hit`, read at `q << 5 | s`; since DI-16, 2026-10-31) where code 0
  scores -128 against every letter in it (`MatrixWalkTable::x_sentinels_stop`). An X-drop
  outside 1..=127 or a word that does not fit runs the reference, and so does the debug output. In
  the unit path the switch acts on the units of the branch-free hit loops (`LOSAT_X_TBXHITBF`): each
  unit takes one of two instances of the loop's function (`scan_unit_bf`, `scan_unit_batched_bf`,
  `scan_unit_sliced_bf`, `const SENT`, chosen once per unit), and the instances without the walks
  are dev's machine code (normalized `objdump`: `scan_unit_bf`, the masked unit scans, `flush`,
  `resolve_early`; the instances with batched extensions store one field more per unit). The units
  of the other hit loop (`LOSAT_X_TBXHITBF=0`, the probes of `LOSAT_TIMING` and `LOSAT_X_STATS`) run
  without it, so timings taken with those variables do not include it; the reference frame loops test
  at every call. wasm32 (LOSAT Web) keeps TBXUNGAP's walks (plan v4 line W).
  Measured 2026-10-25 on the workstation against dev f24d161f, the binary before the switch (node 1
  without other timed runs, CPU 40 at 1 thread and CPUs 40-47 at 8, their HT siblings idle, warmed
  3 s before every run, arm order reversed every other round, one warm-up round + 8 rounds, medians
  of the per-round ratios; every report identical between the arms): AP027131 x AP027133 (code 4)
  1 thread x0.950 (7.56 -> 7.19 s, 8 of 8 rounds below 1; cycles x0.950, instructions x0.924) and
  8 threads x0.968 (1.524 -> 1.474 s, 8 of 8; cycles x0.957); MelaMJNV x PemoMJNVA 1 thread x0.963
  (1.030 -> 0.992 s, 8 of 8; instructions x0.936); Sakai x MG1655[0:1 Mb] (code 11, the batched
  scan in slices) 8 threads x0.979 (11.96 -> 11.71 s, 8 of 8; cycles x0.967, instructions x0.919).
  With `=0` the instructions equal dev's (x1.0000 in all four) and so does the time (x0.997-1.002).
  Max RSS equal but for Sakai at 8 threads, 743 MiB with the switch against 884 MiB for dev and
  `=0` in every round (the switch allocates nothing). On AP027131 x AP027133 at 1 thread the
  extension (31.9 million calls) took 235 cycles a call in dev (`walk_left`, `two_hit` and
  `walk_right`, 30.6 % of the cycles in `perf record`) and 38 fewer with the switch (from the cycle
  totals), about 197; what is left is the latency of the first query loads and the mispredicted
  exits of the walks (as in TBXUNGAP).
- `LOSAT_X_TBXHITBF=1` (TBLASTX, speed; in the strict set from the pull request that added it,
  2026-10-19; `=0` turns it off) the hit loop of the unit scan of the unit path (`scan_unit`, the
  batched scan of `LOSAT_X_TBXEXTBATCH`, and, from the next pull request (2026-10-19), the scan in slices
  of `LOSAT_X_TBXQSLICE`, `scan_unit_sliced`; `tblastx/x_tbx_hitbf.rs`, exactness argument in its header) takes
  the two outcomes most hits end in, too far (`diff >= window`: `last_hit = s + window`) and overlap
  (`diff < wordsize`: nothing), from one subtraction, `d = s - last_hit` (too far: `d >= 0`; a
  candidate: `d + window - wordsize` in `[0, window - wordsize)`, one unsigned compare), and writes the
  cell's next value chosen without a branch (`s + window`, or the value it read);
  a flag-1 cell or a hit past both tests goes to the per-hit body as before (one branch, rarely taken).
  The loop walks the batch by pointer and the prefetch of `LOSAT_X_TBXPREFETCH` reads its entry ahead
  without a bound test. Acts on a unit when `1 <= wordsize < window` and its subject chunk + window is
  below 2^29 (the batch guard `last_hit_fits` already keeps the unit path inside it); a unit outside runs
  the loop before the switch. With `LOSAT_TIMING` or `LOSAT_X_STATS` (their probes count every outcome)
  the units run the loop before the switch, so timings taken with those variables do not include it.
  On AP027131 x AP027133 at one thread 80 % of the hits are too far and 16 % overlap, and the too-far
  test was mispredicted 0.14 times per hit (v4-R6 REPORT). Measured (workstation, node 0, the booked
  cores warmed 3 s before every run, arm order reversed every other round, 6 rounds, medians of the
  per-round ratios against the dev binary before the switch, every report identical, max RSS equal):
  AP027131 x AP027133 (code 4) 1 thread x0.792 (10.59 -> 8.39 s; cycles x0.792, instructions x0.867,
  mispredicted branches x0.537) and 8 threads x0.851; AvCLPV x PsCLPV 1 thread x0.873 and 8 threads
  x0.937; MelaMJNV x PemoMJNVA 1 thread x0.786 and 8 threads x0.874; Sakai x MG1655[0:1 Mb] (code 11,
  the batched scan) 8 threads x0.857 (22.15 -> 19.14 s), measured before `LOSAT_X_TBXQSLICE`, which now
  cuts that query's groups at 8 threads. With the loop in the scan in slices too (workstation node 1, 8
  booked cores warmed 3 s before every run, arm order reversed every other round, 6 rounds, medians of
  the per-round ratios, every report identical, max RSS equal), against dev with the loop in the two
  other scans only: Sakai x MG1655[0:1 Mb] 8 threads x0.784 (15.76 -> 12.34 s; cycles x0.742,
  instructions x0.890, mispredicted branches x0.412); AP027131 x AP027133 8 threads (not cut) x1.001;
  `=0` on Sakai x0.994 (instructions x1.000). With `=0` the reports are the same and the
  loop before the switch runs, but its build differs in inlining from the binary before the switch:
  instructions x1.000 / x0.989 / x0.982 on the three pairs above at one thread, time x0.98-0.99.
  A first form of the test from the v4-R6 probe (`2s - 2 last_hit` against two constants) kept both
  constants on the stack and was 3-4 % slower at one thread.
  On wasm32 (LOSAT Web; the switch acts there too) the loop runs, from 2026-10-20, as
  `x_hits_bf_split`: a function of its own, `x_hits_bf_run`, takes the too-far and overlap hits up
  to the first hit for the per-hit body, which its caller then runs: the same cell reads and writes
  and calls of the per-hit body, in the same order. In the one loop V8's optimized code (TurboFan)
  kept the table, the bounds and the loop pointer in stack slots (the inlined per-hit body needs
  every register): 28 instructions per too-far or overlap hit, 6 of them stack reads and writes;
  split, V8 inlines `x_hits_bf_run` into its caller and unrolls its loop 4 times, 18 instructions
  per hit without stack slots (no value of that loop lives across a call). Native builds keep the
  one loop (same machine code). Measured on Wasm (Node 26, the threaded and serial reactors, every
  arm in one process, warm searches alternating, the `perf` CPU of a stage divided by that of the
  stages no arm changes (linking, SEG, report, statistics, sort) in the same search, per-round
  ratios, every report identical; v4-W DECISIONS 33 and 36): the switch against `=0` on dev
  7e73ec9c, scan stage x0.62-0.75 and whole search x0.77-0.91 (AP027131 x AP027133, MelaMJNV x
  PemoMJNVA, AvCLPV x PsCLPV at 1 and 8 threads, NZ_CP006932 against itself at 8 threads); the split
  loop against the one loop, the loop's function x0.87-0.93 (AP027131 x AP027133 at 1 and 8 threads,
  MelaMJNV x PemoMJNVA at 8 threads) and whole searches x0.97-0.99.
- `LOSAT_X_TBXEXTBATCH=256` (TBLASTX, speed; a value switch, in the strict set from the pull request
  that added it, 2026-10-18; x86-64 only; `=0` turns it off) the unit scan of the unit path
  (`tblastx/x_tbx_extbatch.rs`, exactness argument in its header) holds back the two-hit extensions of
  the subjects where one unit's query-group span (the letters of the query frames the unit's hits fall
  in) has at least the switch's number of Ki letters (2 to 4,194,304; built-in 256; `=2` every unit of
  2 Ki letters or more):
  a hit that passes the two-hit test records its extension in a batch of 16, prefetches the query and
  subject lines the extension reads first and marks its cell (flag 1, the slot in `last_hit`); the
  batch is resolved by `extend_hit_two_hit` when it is full and at the end of the unit, its init HSPs
  saved in slot order, and a later hit on a marked cell resolves that one extension first. Each cell
  goes through the states of `scan_unit` in its order and the unit's list is the same element for
  element. The extension of a hit starts at a random position of the query group, so on a large query
  its first loads wait on memory; the batch lets them arrive while the scan goes on. Acts where the
  window leaves room for the marks (`window >= wordsize + 16`, the default 40). With `LOSAT_TIMING` or
  `LOSAT_X_STATS` (their probes) the units run `scan_unit`, so timings taken with those variables do
  not include the switch; the other subjects run `scan_unit`, whose machine code is the one before the
  switch (`=0` and the binary before it run the same instruction counts). For one query the rule holds
  from about 0.8 Mb; a batch of more than 6 contexts merges contexts into a group (`X_MAX_GROUPS`), so
  smaller queries can reach it there (not measured).
  Measured (workstation, 8 booked cores of a node without other timed runs, alternating rounds after a
  warm-up, medians of the per-round ratios, every report identical, max RSS equal), against the binary
  before the switch: Sakai x MG1655 (code 11, frames of 1.83 M letters) 8 threads x0.940 (112.9 ->
  105.9 s, 4 rounds), Sakai x MG1655[0:1 Mb] x0.926 (21.48 -> 19.89 s; x0.935 against `=0`), Sakai's
  first 0.8 Mb against it (frames of 267 k letters) x0.987; below the rule AP027131 x AP027133 and
  MelaMJNV x PemoMJNVA run the same instructions as before (76.99 G and 9.906 G at one thread). How
  the built-in was chosen (the path on every unit against `=0`, 8 threads): Sakai's first 0.4 / 0.8 /
  1.6 / 3.2 Mb and all of it against MG1655[0:1 Mb] (frames of 133 k / 267 k / 533 k / 1.07 M / 1.83 M
  letters) x1.005 / 0.985 / 0.977 / 0.956 / 0.918; MelaMJNV x PemoMJNVA (96 k letters) 1 thread
  x1.015, AP027131 x AP027133 (67 k, code 4) 1 thread x1.004 and 8 threads x0.999: 256 Ki lies between
  the span that lost and the first that gained. An 8 x i32 AVX2 lane kernel for the batch (the
  walks of 8 extensions side by side) was built and measured, and lost to this scalar resolution on
  every case (branch `exp/a5-lanes-kernel`; v4-A DESIGN_A5 §4.2).
- `LOSAT_X_TBXQSLICE=1` (TBLASTX, speed; a value switch, in the strict set from the pull request that
  added it, 2026-10-19; x86-64 only; `=0` turns it off) the unit scan of the unit path
  (`tblastx/x_tbx_qslice.rs`, rule and exactness argument in its header) cuts a long query group into
  slices of consecutive query offsets, each with its own lookup table (the group's chains filtered to
  the slice), and scans the unit's subject chunk once per slice, in slice order, on the unit's one
  diagonal table, which is not cleared between slices. At one subject position a slice's hits touch
  cells in a window of the slice's width, which moves one cell per position: with slices of 128 Ki
  letters it is 512 KiB and stays in the core's L2, where the whole group's table (4 bytes a cell: the
  group's span plus the chunk; 7 MB for one frame of Sakai) is read from the L3 or, with 8 or 16
  threads, from memory. Each diagonal sees its hits in the same order, the table ends the same cell
  for cell and the unit's records are the same after the chunk's stable sort, so the reports are
  those without it. Unset (or `=1`), a group is cut into slices of at most 131,072 letters where its
  span is longer and its span times the search's threads reaches 6 Mi letters (24 MiB of tables at
  once; for one query, a frame of 1.5 Mi letters at 4 threads, 768 Ki at 8, 384 Ki at 16); below that
  the whole group's scan (with `LOSAT_X_TBXPREFETCH` where the table is large) is as fast or faster.
  `=N` (2 to 2^30) cuts every group of more than N letters into slices of at most N, whatever the
  threads (tuning and tests); at most 64 slices a group. A batch with no cut group runs `scan_unit`,
  whose machine code is the one before the switch (the binary before it, the default and `=0` run the
  same instruction counts: AP027131 x AP027133 76.99 G, MelaMJNV x PemoMJNVA 9.906 G, AvCLPV x
  PsCLPV 43.61 G at one thread); in a batch with a cut group every unit takes the sliced scan (a
  group that is not cut is one slice). In the slices `LOSAT_X_TBXPREFETCH` and `LOSAT_X_TBXEXTBATCH`
  apply their rules to the widest slice instead of the group: with 128 Ki slices neither acts
  (forced on, they lost 1.4-1.8 % at 8 threads). Each slice table has its own backbone (about 460 KB).
  With `LOSAT_TIMING` or `LOSAT_X_STATS` (their probes) nothing is cut.
  Measured (workstation, Xeon 6338, L2 1.25 MiB per core, L3 48 MiB per node; booked cores of a node
  without other timed runs, one warm-up round, then rounds alternating the order, medians; every
  report identical): Sakai x MG1655 (code 11) at 8 threads 106.1 -> 85.6 s against the binary before
  (x0.807; x0.818 against `=0`), at 16 threads 86.5 -> 62.8 s (x0.730), cycles 2.0 -> 1.5 T at 8;
  Sakai x MG1655[0:1 Mb] at 16 threads x0.690, at 8 x0.801; Sakai's first 3.2 Mb (frames of 1.07 M
  letters) at 8 threads x0.961; at one thread, below the rule, the same code. How the rule was chosen
  (slices on every group against `=0`): 64 / 128 / 256 / 512 Ki letters on Sakai x MG1655[0:1 Mb] at
  16 threads x0.696 / 0.690 / 0.716 / 0.693 (max RSS +149 / +26 / +8 / +1 MB), at 8 x0.797 / 0.801 /
  0.836 / 0.797; with 128 Ki, Sakai's first 1.6 Mb at 8 threads (4.3 M letters x threads) x1.005, its
  first 0.8 Mb x0.997, Sakai x MG1655[0:1 Mb] at one thread x1.015 (2 threads x1.006; 4 threads x0.724,
  a node with other runs). Max RSS of Sakai x MG1655 at 8 threads: +51 MB on one node, +262 MB on the
  other (2.0 -> 2.26 GB), equal with `MALLOC_ARENA_MAX=1` (1956 vs 1948 MB): glibc's per-thread arenas
  keep what the units free, and the peak lies in the stage after the scan.
- PAIRLUT (graduated 2026-10-22; no variable, no `LOSAT_X_PAIRLUTSHADOW`; class conditional, the name
  stays in `GRADUATED` for the reach table; it was in the strict set) TBLASTX: the lookup tables of
  the unit path, one per context group, built directly from the query instead of being filtered out
  of the batch's full table, which is then not built: for each group the exact-match index of the
  group's lookup locations, then the count, layout and fill passes of BXLUT's builder
  (`tblastx/lookup/x_pair_lut.rs`, exactness argument in its header), the groups in order on the
  calling thread (tuning only: `LOSAT_X_PAIRLUT_POOL=1` builds them on the search pool, a few ms
  faster per batch and more memory). The reference code holds a stub (the full table's header fields
  and `longest_chain`, no entry). The rule stays: it acts for a batch that PAIRPAR's `decide` sends to
  the units, decided before the lookup build by the function `prepare` decides with; the run stops if
  `prepare` disagrees or if the reference frame loop would scan the stub. The fallback stays: the
  reference builder and `build_groups` run when the unit path does not run (gate, guard, cost rule),
  in diagnostics mode, and when a lookup location does not lie inside one group (`prepare` then builds
  the groups from the full table; no input is known to reach that case, so the branch is kept as the
  fallback). `LOSAT_X_STATS=1` prints `[X_STATS] PAIRLUT batches_direct= groups= fallback=`;
  `LOSAT_TIMING` prints one `x_pairlut` line (the direct build; it is part of `build_lookup`, and
  `sublookup` of the `x_pairpar` line is 0). The shadow (`LOSAT_X_PAIRLUTSHADOW=1`: the full table and
  `build_groups` as well, every group table, `q_min` / `q_max`, the stub's header and the contexts
  compared, the first difference stopping the run) went with the variable; the unit tests of
  `x_pair_par.rs` (`direct_groups_equal_build_groups_*`) still compare the direct groups with
  `build_groups`.
- BXPAR, BXPOOL, BXBATCH, BXCHUNK (graduated 2026-10-22; no variables; class conditional, the names stay in
  `GRADUATED` for the reach table; they were in the strict set) BLASTX scheduling. Their rules stay:
  - BXPAR: query-chunk, redo and per-context SEG parallelism on the search pool of 2+ threads (each stage with its
    own condition: several chunks, several matches, a large batch); 1 thread runs the serial reference stages;
  - BXPOOL: one pool for all query batches of the run at 2+ threads (1 thread: the per-batch path);
  - BXBATCH: with that pool, several small query batches searched together, output in input order (2+ threads and
    several query batches);
  - BXCHUNK: 64 subjects per worker in a wave of the pool (`search::X_BX_CHUNK`; 1 thread: one subject at a time).
  `LOSAT_X_BXWAVE` (tuning of the waves without BXCHUNK) went with them: nothing reached it.
- `LOSAT_X_BXFLOW=1` (BLASTX, speed; in the strict set since 2026-10-17, by the Owner's rule of that
  day that a speed-up lands on by default (DEFAULT_ON_POLICY_20261017 D4): on by default, `=0` turns it
  off; acts with several query batches at two or
  more threads, natively; the Web reactors keep the waves): the queued query batches are searched as a flow
  instead of waves of `threads x 4` (`blastx/native.rs` `x_flow`). Each batch becomes a rayon task once
  `threads` batches follow it; the last `threads` batches of the input are searched with the pool
  inside, side by side (the waves search their tail one batch after another). The reading thread
  writes warnings and reports in the reader's order as the results arrive, keeps at most
  `threads x 4` batches read and unwritten, and runs pending tasks itself while it waits. Reports,
  warnings, their interleaving and the stop at a failing batch do not change
  (`tests/blastx_flow.rs`). With `LOSAT_TIMING=1`: `[TIMING] blastx_flow: batches= alone= main_jobs=
  main_job_time .. blocked ..`. Measured against the yeast proteome (`-db`, `-evalue 1e-5`, ws node 1,
  three rounds): 800 fly 2L cDNAs at 16 threads 0.71-0.84 of the waves' time (and 0.87-0.95 of NCBI
  BLAST+ 2.17.0 `-db`), at 8 threads 0.84-0.87; 50 cDNAs (9 batches) at 16 threads 0.56-0.58, at 8
  threads 0.76; 5 cDNAs (one batch) 1.00. Timed again for the default set (2026-10-17; dev 7a5a7f87 with
  `BXSCANBF` and `BXDIAG` on, the build before the short-subject fix; node 1 CPUs 56-63, 8 threads, one arm at a
  time, the order rotated, 3 rounds; every output the bytes of NCBI BLAST+ 2.17.0 `-db`): 800 cDNAs 7.93-8.18 s
  against 9.59-9.75 s with `LOSAT_X_BXFLOW=0` (0.82-0.84; NCBI 13.09-14.24 s), 50 cDNAs 0.430-0.435 s against
  0.586-0.592 s (0.73-0.74; NCBI 1.09-1.10 s), 5 cDNAs (one batch) 0.086-0.089 s against 0.088-0.092 s (1.00); max
  RSS 69.6-71.8 MB against 66.9-67.2 MB on 800 cDNAs (NCBI 124-128 MB), 43-44 MB against 42-43 MB on 50.
- `LOSAT_X_BXSUBJFLOW=1` or `=k` (BLASTX, speed; in the strict set since 2026-10-17, in the pull request that added
  it: on unless `=0`; `k` from 2 to 65,536 subjects per task, otherwise 64; acts on the subject
  loop of a query batch searched as one chunk, `-db` and `-subject` alike, at two or more threads when there are more
  than `k` subjects, no diagnostic stage traces the subjects and `-window_size` plus twice the longest subject is at
  most 2^28; the chunks of a longer query are searched side by side on the pool, one chunk per thread
  (BXPAR's rule), unless `LOSAT_X_BXCHUNKFLOW` below keeps them in order, each with the flow): the subjects
  are searched as a flow of tasks of `k` consecutive subjects instead of waves of `threads x 64`
  (`blastx/x_subject_flow.rs`, exactness argument in its header). The other threads of the pool search tasks up to `4 x threads` ahead (`utils::xahead`, each helper
  on its own diagonal table) and the loop's thread takes the tasks in order, or searches one itself when no helper has
  started it, and hands the subjects to the sink in subject order with the waves' calls; a wave waits for its slowest
  worker and then for the hand-over on one thread. Above the 2^28 bound (a diagonal table of 1 GiB or more) the waves
  run: there an offset could pass 2^30, where the bit-field wrap of the entries makes a table's history matter.
  Shadow: `LOSAT_X_BXSUBJFLOWSHADOW=1` searches every task again on the loop's own table and stops at a difference
  of the subject results or of the first refused subject; at exit `[X_SHADOW] BXSUBJFLOW tasks_compared=N (subject
  results identical)` (`tests/blastx_subject_flow.rs`). Measured 2026-10-17 on the workstation (16 / 8 physical cores
  of one node with their HT siblings idle; DB-6's BX-t10: the first 10 zebrafish transcripts of the DB-5 q100 set,
  10.9 kb, against Swiss-Prot `-db`, one query batch; a warm-up then 3 rounds, the order reversed in round 2; every
  output the same bytes as NCBI BLAST+ 2.17.0 `-db -num_threads 1`): at 16 threads 5.23-5.30 s with `=0` (the waves),
  2.30-2.33 s with `k = 64` (2.36-2.45 s with 16, 2.44-2.51 s with 4), NCBI `-db` 2.53-2.63 s; at 8 threads
  7.09-7.14 s, 3.887-3.896 s and 4.116-4.122 s; max RSS unchanged (226 MB, NCBI 289-297 MB). The waves' own lines
  (`LOSAT_TIMING`) at 16 threads: 563 waves, 4.27 s, of which the waves' longest tasks 4.00 s against a mean of
  2.02 s (with 1,024 subjects per task still 3.34 against 2.03 s: the work clusters in oid ranges).
- `LOSAT_X_BXCHUNKFLOW=1` (BLASTX, speed and memory; in the strict set since 2026-10-23, in the pull request that
  added it: on unless `=0`; acts on a query batch split into chunks (one query of at least
  19,410 letters, `blastx/split.rs`), `-db` and `-subject` alike, where `LOSAT_X_BXSUBJFLOW` would run (above) and
  there are at least 256 subjects per thread and 2,048 subjects in all): the chunks are searched one after another
  in the serial chunk loop of `search::search_core`, each chunk's subjects as the subject flow on the whole pool,
  instead of side by side, one chunk per thread (BXPAR's rule). A query of five chunks used five threads of
  sixteen. Exactness: the switch only chooses between two schedules of the same chunk loop; the chunk waves hand
  the sink every chunk, subject and chunk end in the order of the serial loop, and the flow is checked by
  `LOSAT_X_BXSUBJFLOWSHADOW` (which now also compares every task of every chunk of a split query;
  `tests/blastx_subject_flow.rs`). Memory: one chunk's lookup table and the pool's diagonal tables at a time,
  instead of one of each per chunk of a wave. The rule (`x_subject_flow::chunks_in_order`, in subjects so that
  `LOSAT_X_BXSUBJFLOW=k` does not move it): with fewer subjects per thread the flow leaves threads idle in each
  chunk, and with fewer subjects in all each chunk's setup on the loop's thread (its lookup table, about 5 ms,
  and the helpers' start) outweighs it at two to four threads. Measured 2026-10-23 on the workstation (16 / 8 /
  4 / 2 / 1 physical cores of node 1 with their HT siblings idle, `--membind=1`, cores warmed 3 s before each
  run, 3 pairs with the order alternated; every output the same bytes across arms, and for bx50k the same bytes
  as NCBI BLAST+ 2.17.0 `-db -num_threads 1`; a probe binary that chose the schedule by a variable): bx50k
  (fly 2L 800,001-850,000, five chunks, against the 30,802 fly proteins `-db -evalue 1e-5`) at 16 threads
  2.44 s with the chunk waves, 0.95 s in order (NCBI `-db` 1.63 s), at 8 threads 2.45 and 1.63 s (NCBI
  2.76 s), max RSS 71 and 34 MiB (NCBI 91-96 MiB); 1 thread unchanged (11.28 / 11.31 s); instructions the same
  (93.0 G); 2L slices of 20-300 kb (2-30 chunks) x0.18-0.91 at 8 and 16 threads, max RSS 31-43 MiB against
  45-182 MiB. With few subjects the waves win and the rule keeps them: stride subsets of the fly proteins with
  512 / 1,024 / 2,048 subjects in order lost up to x1.80 at 16 threads (30 chunks x 512), x1.18 at 8 (512), and
  x1.04 / x1.01 at 2 threads (30 chunks x 512 / 1,024); `AP027131.fasta` against the 543 proteins of
  `AP027133.faa` (`-subject`, 68 chunks) x1.20 at 8, x1.49 at 16 threads. Under the rule no measured case is
  slower; the subsets with 2,048 subjects at 2-8 threads x0.68-0.99, with 4,096 at 16 threads x0.51-0.97.
  This pull request's binary against its base, 3 pairs on the same cores (whole node 1 for 32 threads): bx50k at 32 /
  16 / 8 threads 2.45 -> 0.76 / 0.95 / 1.62 s (x0.310 / x0.388 / x0.656), max RSS 71 -> 34-36 MiB, at one thread
  x1.004 with the same instructions; with one chunk per batch (the 10 transcripts of BX-t10 against Swiss-Prot at 8
  and 16 threads, 1,000 short reads against Swiss-Prot at 16) x0.996-1.000; the `-subject` case above x0.988-1.002
  (kept on the waves).
- `LOSAT_X_BXCHUNKCROSS=1` (BLASTX, speed at 2+ threads; in the strict set from the pull request that added it: on
  unless `=0` or `LOSAT_X_DEFAULT=0`; no condition of its own: it acts where the serial chunk loop of
  `search::search_core` would search two or more chunks of a split query one after another, each with the subject
  flow, which by default is `LOSAT_X_BXCHUNKFLOW`'s in-order case): the chunks run as one flow of tasks across the
  chunks (`x_subject_flow::run_chunks`). Task `c x tasks + t` is task `t` (64 subjects) of chunk `c`; the helpers
  search up to 16 tasks per thread ahead of the loop's thread across a chunk's end, so they no longer stop at each
  chunk's end while that thread hands the chunk's last tasks to the sink and sets up the next chunk (before: one
  `utils::xahead` loop per chunk with 4 tasks per thread ahead, then the next chunk's Karlin blocks and lookup table
  on that thread alone). The tasks offered are those of the loop thread's chunk and the next one. A chunk's setup
  (the loop's statements before its subject loop: `score_block`, the lookup table, the loop thread's diagonal
  table, `subject_parameters`; a function of the chunk alone) is made once, on one thread, by the first thread
  that needs it or ahead by a helper within 2 x threads tasks of the end of the loop thread's chunk, so at most two
  chunks' setups are alive at a time. The first chunk's is made before the helpers start, its lookup table in one
  part per thread (at most 8) on the idle pool (the partitioned build of `LOSAT_X_BXBUILDPAR`); a later one is not
  built in parts, because a helper waiting at the end of such a build for a part another thread took could run a
  helper loop not yet started and wait there for the setup it is building. Exactness: a task's value depends on its chunk and subjects alone (the module
  comment's argument for the diagonal tables holds across chunks), the loop's thread makes the sink calls in the
  serial loop's order (per chunk `chunk_start`, its subjects, `chunk_end`; a chunk without a valid context gets
  `chunk_start` only), and a setup that fails returns its error after that chunk's `chunk_start`, where the loop met
  it. `LOSAT_X_BXSUBJFLOWSHADOW` compares every task of every chunk (`tests/blastx_subject_flow.rs`: 30,000 letters
  x 3,422 subjects at 2 / 3 / 8 threads, with `LOSAT_X_BXBUILDPAR=0`, and the switch `=0`, all equal to one
  thread). Why 16 tasks per thread: in the ordered flow a task much slower than the others (a subject with many
  hits) holds the loop's thread at its place while the helpers search the rest of the window; with 4 per thread
  they ran out of tasks and spun in `sched_yield` (fly 2L 50 kb at 16 threads: 1.6 s of system time in a 0.69 s
  run); with 16 they seldom do (0.08 s), and the tasks held at once stay few (only subjects with hits are kept).
  The one-chunk subject flow (`x_subject_flow::run`) keeps 4 per thread. Measured 2026-10-27 on the workstation
  (node 1, the first n physical cores with their HT siblings idle, `--membind=1`, cores warmed 3 s per case,
  rounds in rotated arm order, medians of the per-round ratios; every output the same bytes as dev's and, for the
  50 kb slice, as the progress board's d26), dev dfe8a461 against this change (d7eb5b30, the gate's binary): fly
  2L 800,001-850,000 (five chunks) x the 30,802 fly proteins `-db -evalue 1e-5` at 16 threads 0.684 -> 0.570 s
  (x0.833, 8 of 8 rounds), system time 1.58 -> 0.08 s, max RSS 16 -> 25 MiB (NCBI `-db` 91-96 MiB); 32 threads
  0.555 -> 0.401 s (x0.723); 8 threads x0.882; 4 / 2 / 1 threads x0.947 / x1.003 (rounds 0.997-1.010) / x1.000;
  S(16) (one thread / 16 threads) 11.43 -> 13.71. 2L slices of 20 kb (two chunks) and 150 kb (15 chunks) at 16
  threads x0.894 and x0.799 (system time of the 150 kb slice 4.93 -> 0.13 s). The cross-chunk flow with 4 tasks
  per thread ahead (a prototype) gave only x0.973 at 16 threads (system time 1.59 s); 32 or 64 per thread gave the
  same times as 16. 800 cDNAs x the yeast proteins at 16 threads `-subject` x1.000 (0.999-1.001), `-db` x1.000
  (1.000-1.003), `-db` at 32 threads x0.997 (0.995-1.006); max RSS at 16 threads +7 MiB (70 -> 79, 65 -> 72), at
  32 threads 90 -> 91 (NCBI 143-154 MiB; the v4-MF floor of the 800-cDNA `-db` row at 16 threads is about 85).
- `LOSAT_X_BXTASKRES=1` (BLASTX, speed at 2+ threads; in the strict set from the pull request that added it: on
  unless `=0`; no condition of its own: it acts on every subject flow, `x_subject_flow::run` and `run_chunks`, i.e.
  where `LOSAT_X_BXSUBJFLOW` runs): a flow's task ends after `k` subjects (as before) or before the subject that would
  take the task's residues past `k` times the mean subject length (`XTasks::Cut`), so a run of long subjects stored
  next to each other is spread over several tasks. The bounds come from the subjects' lengths (the records, the
  stored sequences, or with `-db` read in blocks the index, 16,384 lengths per read); with `-db` read in blocks they
  are made once per run for each `k` and shared (`XStreamDb::x_task_starts`). Where the lengths cannot be read (an
  index read error) the grid of `k` subjects is kept, so the reader's error comes where it did before. Exactness:
  only which thread computes which subjects changes; a task's value depends on its chunk and subjects alone and the
  loop's thread makes the sink calls in subject order (`x_subject_flow.rs` header); `LOSAT_X_BXSUBJFLOWSHADOW`
  compares every task (`tests/blastx_subject_flow.rs`: task counts of the cut and, with `=0`, of the grid; a unit
  test of the bounds). Why: the ordered flow ends when its last task ends. Fly 2L 50 kb x the 30,802 fly proteins:
  task 400 of each chunk (subjects 25,600-25,663, the isoforms of dumpy, about 20-23 kaa each) holds 248,369
  residues against a median task of 38,410 and took 42-60 ms against a median of 4 ms (local probe), and in the last
  chunk it ran alone for the last 21 ms of a 0.58 s run at 16 threads (about 2 of 16 CPUs busy). Cut by residues the
  largest task holds 42,344 residues (547 tasks per chunk instead of 482) and the tail is gone. Measured 2026-10-28
  on the workstation (node 1, the first n physical cores with their HT siblings idle, `--membind=1`, cores warmed
  3 s per case, rounds in rotated arm order, medians of the per-round ratios; every output the same bytes as dev's,
  the 50 kb slice also as the progress board's d26), dev 4403cc12 against this change: 50 kb slice x the fly
  proteins `-db` at 16 threads 0.580 -> 0.562 s (x0.965, 0.928-0.983, 8 of 8 rounds), 32 threads x0.910, 8
  threads x0.997, one thread x1.000 (S(16) in the same run 13.51 -> 13.95); 20 kb and 150 kb slices at 16
  threads x0.910 and x0.976; 800 cDNAs x the yeast proteins at 16 threads `-subject` x0.995, `-db` x0.998, `-db` at
  32 threads x0.991; with the switch off the same binary gives dev's times (x0.987-1.003). Max RSS at 16 / 32
  threads of the 50 kb slice 26 -> 20 / 28 -> 19 MiB. A second run with the gated binary (12 rounds): 16 threads
  0.571 -> 0.549 s (x0.962, 0.959-0.968, 12 of 12), S(16) 13.70 -> 14.29; the capped few-subject rows (800 cDNAs x
  2 / 8 / 50 yeast proteins, AP027131 x 10 / 50 proteins, 16 and 32 threads) x0.986-1.031, each within its
  per-round spread (`k` is 1 on the 2 and 8 protein rows, so the cut changes no task there).
- `LOSAT_X_BXKAPPAFIT=1` (BLASTX, memory at many threads; in the strict set from the pull request that added it: on
  unless `=0`; acts where a batch end at 2+ threads redoes composition for 2 or more matches, fewer than the
  threads): the batch end makes `min(threads, matches)` composition-redo workers (`runtime.rs`, `x_kappa_workers`)
  instead of one per pool thread; with fewer matches than threads one speculative batch holds every match and each
  match takes the worker of its position (with as many matches or more, the thread's own, as before). Why: each
  worker is a full redo state (`KappaState::x_for_thread`: the redo parameters, a gapped-alignment scratch and a
  composition workspace), made per batch end; a `-subject` search of 800 fly cDNAs against a few yeast proteins at
  32-64 threads made 32-64 of them at every batch end for a handful of matches, from several batch ends at a time
  (the capped batch flow), and the threads' arenas kept that memory. Exactness: a worker's state does not change a
  match's result (`x_speculate_match` redoes in stream order a result that read a predecessor's matrix); with fewer
  matches than threads `x_batch` (32 per worker, at least 64) holds every match, so positions name distinct workers;
  `tests/blastx_unit_cap.rs` compares reports with the switch on and off at 2, 4 and 8 threads for 3, 10 and 40
  subjects and the batch flow, and checks both sides of the reach. Measured 2026-10-29 on the workstation, dev
  ec91ea18 against this change, max RSS as the peak `VmHWM` sampled every 2 ms (`ru_maxrss` at 24-64 threads reads
  9-31 MiB below it); every output the same bytes in every arm. Against NCBI BLAST+ 2.17.0 (node 0, 3 rounds, NCBI
  `-subject` at one thread in the same window, the same bytes as NCBI's): 800 fly cDNAs x 8 / 50 / 200 yeast proteins
  `-subject` at 32 threads 82.2 / 87.2 / 97.9 -> 79.5 / 81.6 / 83.3 MiB (NCBI 81.3 / 82.0 / 83.7), at 64 threads 112.6
  / 131.5 / 161.4 -> 106.8 / 113.4 / 116.5, still above NCBI at 48 and 64 threads (the rest is per-thread live memory:
  kept diagonal tables, scratch, 48-64 arenas). Time and max RSS in the same runs (node 1, the first n physical cores
  with their HT siblings idle except at 64 threads, `--membind=1`, cores warmed 3 s per case, arms dev / new / `=0`
  rotated per round, medians of per-round ratios; run 3 12 rounds on 26 cases, run 4 16 rounds with an exact end time
  on 10): where the switch acts, max RSS 800 cDNAs x 200 proteins at 64 / 32 / 16 threads x0.716 / x0.855 / x0.968, x
  50 proteins x0.865 / x0.930 / x0.981, x 8 proteins at 64 / 32 threads x0.954 / x0.972 (lower in 12 of 12 and 14-16
  of 16 rounds), the proteome `-db` at 64 threads x0.933, AP027131 x 10 proteins at 16 threads x0.989; AP027131 x 50
  proteins at 64 threads and the proteome `-subject` at 16 threads x1.002 / x1.007 (new / `=0` x1.001 / x1.010: batch
  ends that are not the peak). Time x0.962-1.002 on every case (800 cDNAs x 2 / 8 / 50 / 200 proteins, AP027131 x 10 /
  50, the proteome `-subject` / `-db`, AP027131 x AP027133, fly 2L 50 / 150 kb `-db`; 1-64 threads), each within its
  per-round spread; the three run-3 readings above x1.01 (x 2 proteins at 16 threads x1.014, AP027131 x AP027133 at 32
  x1.017, x 8 proteins at 64 x1.011) were run 3's clock, which moved in 2.35 ms steps, and read x0.990-0.997 in run 4.
  The switch is not reached on 800 cDNAs x 2 proteins (no batch end with two matches) and takes the reference side on
  AP027131 x AP027133 and the fly rows.
- `LOSAT_X_BXDIAGKEEP=1` (BLASTX, memory; in the strict set from the pull request that added it, 2026-10-24: on
  unless `=0` or `LOSAT_X_DEFAULT=0`; acts where `LOSAT_X_BXSUBJFLOW` searches a chunk's subjects, and on the chunk
  subject tasks of `LOSAT_X_BXUNITCAP` below, where `-window_size` plus twice the longest subject is at most 2^28):
  each helper of the subject flow runs on the diagonal table its thread keeps (one of the chunk's length and window,
  else a new one, `search/x_bx_mem.rs` `XDiagSlot`) and gives it back when it ends, instead of one table per pool
  thread made by the searching thread for each chunk. A batch searched with the pool inside (the batch flow's last
  batches, a batch with a long query) made `threads` tables of next_pow2(letters + window) x 4 bytes on one thread,
  128 KiB each for a batch of 10,002 bases: 4 MiB per batch at 32 threads, kept in that thread's glibc arena.
  Exactness: a subject's computation depends on the subject and the chunk alone, not on the table's history, while
  no diagonal offset passes 2^30 (`x_subject_flow.rs` module comment, the bound the flow already requires).
- BXLUTPOOL (graduated 2026-10-27, in the set from the pull request that added it, 2026-10-24; no variable; always
  on; every word-size-3 table built in place, outside diagnostics mode; since `LOSAT_X_BXPOOLN1`, 2026-10-26, only in
  searches of two threads or more; the unpooled build of a one-thread search, `x_build_direct`, stays) (BLASTX,
  memory): the lookup tables
  are built with one set of neighbour lists shared by all threads (filled on first use as each thread's memo was,
  `tblastx/lookup/x_lut_direct.rs` `XSharedMemo`), and on vectors kept between builds: a table's backbone, overflow
  and pv vectors when the table is dropped (`blastx::seed::Lookup`'s drop), and the build's counts, layout and fill
  vectors and its exact backbone (chains and slabs) after the build, each cleared and zero-filled to the new table's
  sizes as calloc gives them, grown to n + n/4 where short. Each thread kept its own neighbour memo (about 1.1 MiB)
  and every thread that built a table kept a table's worth (3-6 MiB) in its arena after the table was freed. The
  pool holds as many sets as tables and builds were alive at once. The table equals `x_build_direct`'s (unit test
  `pooled_tables_equal_the_reference_on_kept_vectors`: smaller and larger tables in turn, two alive at once, a
  change of matrix).
- `LOSAT_X_BXUNITCAP=1` or `=K` (BLASTX, memory; in the strict set from the pull request that added it, 2026-10-24:
  on unless `=0` or `LOSAT_X_DEFAULT=0`; a value switch, K from 2 to 4,096, otherwise 8; acts with more than K
  threads, `LOSAT_X_BXDIAGKEEP` on, and K units with tasks for every thread): at most K query batches are searched
  at a time by the batch flow (`LOSAT_X_BXFLOW`), each with the pool inside (its subjects as the subject flow),
  instead of one batch per thread. The flow counts the batches being searched, not the batches not yet written:
  it reads up to 4K batches ahead of the writing (`x_bx_mem::flow_bound`), so a slow batch at the head of the
  output does not stop the others; a search that ends spawns the next batch itself instead of waiting for the
  writing thread; and the subject flow's helpers inside a capped batch leave once every task of the chunk is
  claimed (`xahead::Ahead::x_new_leaving`) instead of waiting, spinning and yielding, for the batch's thread. The
  chunks of a split query (`LOSAT_X_BXPAR`) are searched in waves of K, each chunk's subjects as tasks of 16 on the
  pool, instead of waves of one chunk per thread. The lookup tables alive at once are then at most K (and K more for
  a split query inside a batch of the flow), whatever the number of threads. Exactness: schedule only; the batch
  flow writes every batch where the waves wrote it and searches a capped batch with the function it uses for the
  batches it searches with the pool inside; a chunk's subject tasks keep the subjects with something for the sink
  in subject order and the replay gives the sink the calls of the serial chunk loop (`tests/blastx_unit_cap.rs`:
  reports equal to all of them `=0` with K = 2 and 3 at 4 threads, each part off, 16 batches with K = 2, a split
  query, `-db`). With `LOSAT_X_BXCAPFILL` (below) the capped units' subject tasks are smaller with few subjects
  (at most 64 in the batch flow, 16 in the chunk waves, as here) and the cap needs no task for every thread. Measured
  2026-10-24 on the workstation (node 1, CPUs 32..32+n-1 with their siblings and the rest of the node idle,
  `--membind=1`, 5 rounds after a warm-up in rotated order, medians; max RSS = max(`%M`, sampled VmHWM); every
  report the same bytes across arms; the pull request's first commit against dev 7164b5f7): bx_800 (the 800 fly
  2L cDNAs against the 6,600 yeast proteins, `-subject -evalue 1e-5`) at 1 / 8 / 16 / 32 threads x0.994 / 0.967 /
  0.976 / 0.953, max RSS 29 / 79 / 132 / 272 -> 29 / 71 / 83 / 108 MiB (NCBI 2.17.0 `-subject` 154 MiB); the same
  with `-db` at 32 threads x0.960, 265 -> 108 MiB; bx_ap (`AP027131.fasta` against the 543 proteins of
  `AP027133.faa`, `-query_gencode 4`: one query of 68 chunks) x1.000 / 0.973 / 0.944 / 0.964, 32 / 68 / 102 / 164
  -> 35 / 65 / 71 / 85 MiB (NCBI 95 MiB; at one thread the pool keeps a table and a build set that glibc unmapped
  before). The parts (prototypes, against all three `=0`): BXDIAGKEEP and BXLUTPOOL alone x0.974-0.994 with bx_800
  at 32 threads 274 -> 210 MiB (one table per thread still alive at once); BXUNITCAP with K + 1 unwritten batches,
  helpers waiting for the batch's thread and only the writing thread spawning batches cost bx_800 x1.07-1.14 at
  16-32 threads, which the three parts of the capped flow above remove.
- `LOSAT_X_BXCLUTCAP=1` or `=K` (BLASTX `-task blastx-fast` and `-word_size 5`, memory; in the strict set from the
  pull request that added it, 2026-10-31: on unless `=0` or `LOSAT_X_DEFAULT=0`; a value switch, K from 2 to 4,096,
  otherwise 3; acts where the compressed lookup table is used and `LOSAT_X_BXUNITCAP` is on): `LOSAT_X_BXUNITCAP`'s
  K is at most this K, or ceil(threads / subjects) where that is more, because K units of few subjects keep the
  threads busy only with K x subjects tasks (`blastx/x_bx_mem.rs`). The order of the report is the writing order of
  the units and of the subjects' tasks, as for `LOSAT_X_BXUNITCAP`, so reports do not change
  (`tests/blastx_fast_flat.rs`: 3, 4 and 8 threads with the switch on, off and `=3`, both sides of the rule).
  This pull request also gives BLASTX's compressed table the flat build of `LOSAT_X_BPCLUT` (`=0` keeps the
  reference table, 29 MiB of fixed cells per unit). Choosing K, binary f14f3ad5, ws node 1, 3 rounds, K through
  `LOSAT_X_BXUNITCAP`, 8 threads, time and peak `RssAnon` of K = 8 / 4 / 3 / 2 as parts of the reference table at
  K = 8: Trinity x Swiss-Prot 0.875 0.585 / 0.853 0.370 / 0.860 0.304 / 0.869 0.242 (231 MiB; NCBI 47);
  200 *D. melanogaster* cDNAs x Swiss-Prot 0.848 0.753 / 0.815 0.446 / 0.815 0.413 / 0.824 0.410; AMRProt with
  ABRicate's call 0.775 0.509 / 0.769 0.505 / 0.775 0.506 / 0.936 0.395; Trinity x four proteins (4 subjects)
  0.609 0.640 / 0.720 0.406 / 0.831 0.303 / 1.089 0.196; EDL933 x AMRProt 0.554 0.652 / 0.555 0.603 / 0.556 0.594 /
  0.556 0.690. K = 2 is slower than the reference table with 4 subjects; K = 3 is not, and is the smallest K of
  which that holds in all these cases, also at 16 threads (0.895 0.315, 0.763 0.501, 0.541 0.623).
  Final binary 0f22ecd0 against dev caf4c166 (reports equal to NCBI 2.17.0 in every arm), time and peak `RssAnon`
  as parts of dev, 3 pairs, ws node 1: Trinity x Swiss-Prot 4 threads 0.896 0.527, 8 threads 0.868 0.307, 16 threads
  0.906 0.422; cDNAs 8 threads 0.825 0.392; AMRProt 8 threads 0.746 0.516, 16 threads 0.727 0.517; rhodopsin x four
  proteins 8 threads 0.772 0.580; Trinity x four proteins 8 threads 0.856 0.360, 16 threads 0.717 0.385; EDL933 x
  AMRProt 8 threads 0.562 0.597, 16 threads 0.542 0.590. One thread (no cap): Trinity 0.928 0.851, AMRProt 0.843
  0.572, rhodopsin 0.773 0.578. Against `LOSAT_X_BXCLUTCAP=0` the cap costs time only with few subjects (Trinity x
  four proteins, 0.2-0.3 s, 8 threads 0.291 s vs 0.218 s, 16 threads 0.229 s vs 0.193 s, `RssAnon` 0.56 and 0.61 of
  that arm); cases of a second or more stay within about 2 %, and 0.1 s cases change sign between binaries (rhodopsin
  x four proteins, 8 threads: 0.119 s vs 0.110 s here, 0.095 s vs 0.114 s on 958e1293). Still above NCBI in `RssAnon`: Trinity x Swiss-Prot 8 threads 72 vs 47 MiB,
  EDL933 x AMRProt 16 threads 484 vs 124 MiB. BLASTX without `-word_size 5`: 8 pairs, cycles 0.9971 (byte copy of dev
  1.0005) on 60 Trinity transcripts x Swiss-Prot at one thread, instructions x1.00001, `Lookup::Compressed` keeps one
  variant so that the standard table's loops keep their match (a separate variant cost x1.003 to x1.006 there).
- `LOSAT_X_BXCLUTPOOL=1` (BLASTX `-task blastx-fast` and `-word_size 5`, memory; in the strict set from the pull
  request that added it, 2026-11-01: on unless `=0` or `LOSAT_X_DEFAULT=0`): the flat compressed lookup tables are
  built on vectors kept between builds: a table's starts, hits and pv when the table is dropped
  (`blastx::seed::Lookup`'s drop): the pv cleared and zero-filled as calloc gives it, the starts and hits only sized,
  since the bucket fill writes each of their entries before any is read (zero-filling them again cost the small
  searches 1-5 % of their user cycles); the hits grow to n + n/4 where short (`tblastx/lookup/x_clut_direct.rs`).
  A split query builds a table per chunk on whichever thread needs it, and each of those threads kept a table's worth
  and the build's record of the appends in its glibc arena after they were freed, about 20 MiB per thread
  (`MALLOC_ARENA_MAX=1` removed it). A set comes back only against one lent to a pooled build, so a table built
  without the pool (`LOSAT_X_BPFFILL=0`, a fill declined before it took a set) adds nothing, and the pool holds at
  most as many sets as pooled tables were alive at once (unit test `the_flat_pool_keeps_only_sets_it_lent`). The
  record of the appends is not kept (three builds side by side would hold three records through the search): on
  Linux with glibc its capacity is above glibc's largest mmap threshold, so it is always mapped and unmapped when
  dropped, and only the pages written are resident. Those pages are faulted in again on each build: EDL933 x AMRProt
  at one thread 873k minor faults against dev's 57k, sys +1.4 s, time x1.024 (user x0.988); keeping one record
  between builds removed the faults but raised ABRicate's call x AMRProt at 8 and 16 threads above dev (43.4 and
  45.6 against 40.0 and 39.5 MiB), so it was not taken. The table equals the unpooled one in every field (unit test
  `pooled_flat_tables_equal_the_unpooled_ones_on_kept_vectors`); reports: `tests/blastx_fast_flat.rs`.
- `LOSAT_X_BXENDMEM=1` (BLASTX, memory; in the strict set from the pull request that added it, 2026-11-01: on unless
  `=0` or `LOSAT_X_DEFAULT=0`): at the end of a query batch, no table of word hashes per context for the
  near-identical test of composition-based statistics (8 bytes per letter of six frames: 84 MiB for a 5.6 Mb query)
  and no copy of the batch's residues in its results, which the report does not read (21 MiB). The test reads the
  same words from the context's residues: a hash is its word's letters in 5 bits each and residues are below 32, so
  two hashes are equal exactly where the 8 letters are; the table's last entry stays 0 as NCBI's loop leaves it
  (`redo_alignment::x_test_near_identical_on_residues`; unit test
  `near_identical_on_residues_equals_the_table`, also at ranges not starting at 0, where NCBI reads the table at the
  range's offsets). Other programs keep the table: their queries are proteins of a few thousand letters. With
  `LOSAT_X_BXCLUTPOOL`, EDL933 (5.6 Mb) x AMRProt, ABRicate's call, peak `RssAnon` at 1 / 4 / 8 / 16 threads
  151 / 224 / 312 / 466 -> 71 / 102 / 127 / 168 MiB (NCBI 111 / 124 / 124 / 124), time x1.024 / 0.976 / 0.991 / 0.995
  (ws node 1, 3 rounds, binary c31d3d586 against dev 88bd1240d).
- `LOSAT_X_BXCAPFILL=1` or `=f` (BLASTX, memory with few subjects; in the strict set from the pull request that
  added it, 2026-10-25: on unless `=0` or `LOSAT_X_DEFAULT=0`; a value switch, f from 2 to 64 tasks per thread,
  otherwise 4; acts where `LOSAT_X_BXUNITCAP` caps): the capped units keep the threads busy with few subjects
  too (`blastx/x_bx_mem.rs`). A capped unit's subjects run as tasks of ceil(subjects / (f x ceil(threads / K)))
  subjects, at least 1 and at most the size without the cap (64 for the batch flow's subject flow, 16 for the
  chunk waves of a split query; a set keeps that size from 64 x f x ceil(threads / K) subjects, 512 at 16
  threads), and the cap holds wherever there are more threads than K: without this switch it also needed K units
  with a task for every thread, so with few subjects (fewer than 64 in the batch flow, where a batch had no subject
  flow) it stayed off and every thread searched its own batch with its own lookup table. And the capped batch flow
  prepares each query batch of at most two search units' letters (translation, query SEG, Karlin blocks:
  `reader_batch::prepare_lone_unit`, what `search_core` runs when a batch prepares itself) as a task of its own when
  the reader delivers it, as it prepares the units of a reader batch of several, so that the cap holds only the
  lookup table's build and the subject search to K at a time. Before, a batch of 64 subjects or fewer had no
  subject flow, so the cap left the batch flow alone and each thread searched its own batch with its own lookup
  table (about 6.5 MiB each): 800 cDNAs against 50 yeast proteins at 16 threads took 128 MiB, NCBI `-subject` (one
  thread) 80 MiB; with the cap forced on (`LOSAT_X_BXSUBJFLOW=4`) the preparation, about half of the work at one
  thread, ran on K threads: 68 MiB but x1.59. Exactness: schedule only; a batch prepared ahead is prepared by the
  statements its search runs otherwise, from the same records and options, and a preparation that fails is
  dropped, so the search prepares again and fails where it failed before (`tests/blastx_unit_cap.rs`: 40 proteins
  with K = 2 and 3, f = 2, 4 and 8, a split query of 10 proteins, all equal to every part `=0`). With
  `LOSAT_TIMING` the preparation ahead is not in the batch's `QuerySetup` clock. Measured 2026-10-25 on the
  workstation (node 0, CPUs 0-15 without siblings, untimed, 3 rounds in rotated order, medians, the same bytes in
  every arm; dev f24d161f against the pull request's first commit): 800 cDNAs x 50 proteins at 16 threads
  0.170 -> 0.167 s, 127.5 -> 70.9 MiB; x 200 proteins 0.310 -> 0.229 s (x0.739), 71 -> 74 MiB; AP027131 x 50 / 200
  proteins x0.958 / x0.997; the full sets (800 cDNAs x 6,600 proteins `-subject` and `-db`, AP027131 x 543) at 16
  threads x0.996-1.009, +2-3 MiB (the preparations held ahead); at 8 threads (no cap) x0.99-1.01. With the cap
  wherever threads > K, `LOSAT_X_BXBUILDPAR` and `LOSAT_X_BXCTXONCE` (workstation node 1, 2026-10-27, dev d40021ba
  against the pull request, 7-9 rounds in swapped order; peak = max(time %M, VmHWM), median of 3, untimed): peak
  MiB at 16 / 32 / 64 threads, 800 cDNAs x 2 proteins 97 / 164 / 284 -> 61 / 78 / 99, x 8 111 / 187 / 315 -> 64 /
  82 / 112, x 50 127 / 287 / 509 -> 67 / 90 / 131, x 200 71 / 93 / 304 -> 68 / 99 / 162; AP027131 x 2 / 10 proteins
  at 64 threads 222 / 229 -> 97 / 99; the full sets within 4 MiB. Time x0.48 (x 200 proteins at 32 threads) to
  x1.005, the full sets x0.991-0.998, except x 2 proteins at 32 threads x1.045 (0.100 -> 0.106 s).
- `LOSAT_X_BXBUILDPAR=1` (BLASTX, speed with few subjects at many threads; in the strict set from the pull request
  that added it: on unless `=0` or `LOSAT_X_DEFAULT=0`; acts where `LOSAT_X_BXUNITCAP` caps): the lookup table of a
  capped unit (a batch of the capped batch flow, a chunk of the capped chunk waves) is built in ceil(threads / K)
  parts, at most 8, on the threads the cap leaves it (`tblastx/lookup/x_lut_direct.rs` `x_partitioned_passes`,
  `x_counted_sources`, `x_claimed`). The present source cells (the exact backbone's cells, ascending, and their
  offsets in chain order) come from a stable counting sort of the words of `blast_lookup_index_query_exact_matches`'s
  loop instead of chains. They are cut into ranges of about equal appends. (a) Each range counts its appends per
  target cell. (b) Each stripe of target cells (a multiple of 32 cells) sums the counts, sets num_used and the pv
  bits, and adds up its overflow size and longest count. (c) After the overflow of the stripes before it, each
  stripe gives its overflow cells their cursors and each range its first position in every cell. (d) Each range
  writes its appends at its own positions. The entries of a cell come from the ranges in source order and, within a
  range, in the reference's order, so every cell holds the reference's entries in the reference's order (argument
  in the comments; `partitioned_tables_equal_the_reference` against the reference builder). The caller and the
  helpers it spawns claim the parts by an atomic counter, and the caller waits for parts others run without taking
  other work of the pool: with rayon's `join`, a waiting batch thread took other batches' searches and the build
  was no faster. Before, with K units at a time the builds were the part of a unit no other thread could take
  (2.7 ms of a 3.6 ms batch at 32 threads). Build alone (workstation node 0, a 10,002-base batch's six frames,
  threshold 11): 7.8 ms on one thread, 1.84 ms in 4 parts on 4 threads, 1.72 ms in 8 on 8 (a form by target stripes,
  where every stripe walks all source cells: x0.68 at 4). Measured (workstation node 1, against dev, 800 cDNAs x 8
  yeast proteins at 32 threads): the cap with serial builds x1.28 (2026-10-26, a probe), with this build x0.985
  (2026-10-27, 7 rounds in swapped order).
- `LOSAT_X_BXCTXONCE=1` (BLASTX, speed and memory; in the strict set from the pull request that added it: on unless
  `=0` or `LOSAT_X_DEFAULT=0`; acts where the `LOSAT_X_BXDIAG` fused loop runs): the fused loop's context table of a
  query batch (`XContextTable`, the binary search's index at every query offset) is made once per lookup table and
  kept with it (`XKeptPresence`), instead of once per thread and batch. In the capped batch flow each batch's
  subjects run on several threads, and every thread made the table again: 7 % of cycles and 2.9 MiB of page faults
  at 32 threads (800 cDNAs x 8 yeast proteins). Exactness: the table is a function of the batch's contexts alone.
- BXDBSTREAM (graduated 2026-10-27, in the set from the pull request that added it, 2026-10-22; no variable
  (it was a value switch, the block size in KiB, 2 to 1,048,576), no `LOSAT_X_BXDBSTREAMSHADOW` and no
  `[X_SHADOW] BXDBSTREAM` line; always on, blocks of 256 KiB) (BLASTX `-db`, memory): the search reads
  the database's sequences where it searches them, a block at a time, instead of the copy of the whole database that
  `BlastxDatabase::read_stores` made before the first query batch (deleted; `blastx/x_db_stream.rs`, exactness
  argument in its header). The database is checked once before the first batch with the reads and checks of the copy
  (`x_check_stores`: the same error fails the same batch with the same text) and nothing of it is kept but the length
  of its longest sequence. Each pass of the preliminary search over the database (a query chunk searched on one
  thread, a query chunk searched side by side) reads blocks of 256 KiB in ordinal order and drops each before the next; a task of
  `LOSAT_X_BXSUBJFLOW` reads its subjects with one read per volume, a wave of the pool reads its subjects before its
  fork. The sink reads the lengths of the subjects in its order from the index (16,384 per read) and, for `-ungapped`,
  their residues in blocks. The composition redo (speculated and serial), the traceback of `-comp_based_stats 0` and
  the reports' records read each subject they need with one read. Every read gives the stored bytes of the copy, so
  every report is the same. Before it graduated, a shadow mode (`LOSAT_X_BXDBSTREAMSHADOW=1`) made the copy too and
  compared every block, window of lengths and read with it, and the check's error with the copy's; it went with the
  copy (`tests/blastx_db_stream.rs` compares the reports across thread counts and schedules and the unit test of
  `x_db_stream.rs` compares each read with `SeqDb::read_sequence`).
  Measured 2026-10-22 (workstation, every report identical): peak RssAnon of `uniref50.03` x 201 reads at 8
  threads 49.3 MiB against 2,888.2 MiB with the copy (`=0`) (NCBI: 61 MiB at 1 thread), untimed. Timed on node 1, warmed 3 s,
  8 alternated pairs: against `=0` alone, Swiss-Prot x 201 reads x0.992 at 8 threads and x0.974 at 16,
  `uniref50.03` x 201 x0.996 and x0.992, but x1.0067 at 1 thread (7 of 8 pairs above 1: cycles:u +0.8 to +1.2 %
  at blocks of 16 KiB to 1 MiB, in the scan kernels, not in the new code). With `LOSAT_X_BXEVALEMPTY`, against the
  binary before both: Swiss-Prot x 201 x0.980 / x0.975 / x0.966 at 1 / 8 / 16 threads, `uniref50.03` x 201 x0.969 /
  x0.950 at 8 / 16, no pair above 1. Max RSS there: Swiss-Prot 15.6 / 31.3 / 38.4 MiB against 219.4 / 231.4 /
  239.9, `uniref50.03` 52.3 / 63.1 MiB against 2,825.7 / 2,835.7.
- BXEVALEMPTY (graduated 2026-10-27, in the set from the pull request that added it, with BXDBSTREAM, 2026-10-22; no
  variable, no `LOSAT_X_BXEVALEMPTYSHADOW` and no `[X_SHADOW] BXEVALEMPTY` line; always on) (BLASTX, speed): the gapped evaluation of a subject without
  preliminary HSPs (`statistics::evaluate_gapped`, which the sink calls for every subject, and the composition redo
  and the traceback call for theirs) returns at once the list that both of its ways give for no input
  (`link_uneven`'s first lines, and its loop: no HSP, best e-value 0), without the Karlin and Gumbel lookups and the
  per-context vectors it built first. Most subjects of a database have no HSP: on Swiss-Prot x 201 reads at one
  thread the evaluation was 2.1 % of the run's user cycles (`perf record`, workstation). Before it graduated, a
  shadow mode (`LOSAT_X_BXEVALEMPTYSHADOW=1`) also evaluated every empty input the full way and compared; it went
  with the full evaluation of no input. Measured 2026-10-22 (workstation, with
  `LOSAT_X_BXDBSTREAM`, against the binary before both): Swiss-Prot x 201 reads `-db` at one thread, instructions:u
  x0.9877 and cycles:u x0.9835 (`perf stat`, 3 runs each; x0.9873 and x0.9878 against both switches `=0`), wall
  x0.980 (node 1, 8 alternated pairs, none above 1); `-subject` (201 reads x PemoMJNVA) instructions:u x0.9993.
- TBNDBSTREAM (graduated 2026-10-27, in the set since 2026-10-24 (#443); no variable, no `LOSAT_X_TBNDBSTREAMSHADOW`
  and no `[X_SHADOW] TBNDBSTREAM` line; always on; every `-db` search) (TBLASTN `-db`, memory): the search reads each
  database sequence where it searches it instead of the copy of
  the whole database that `TblastnDatabase::read_stores` made before the first query batch (`read_stores` is deleted;
  `common::database::read_nucleotide_stores` stays for TBLASTX; `common/x_db_stream.rs`, exactness argument in its header). The database is checked once before the first batch
  with the blocks and checks of the copy (`common::database::x_check_nucleotide_stores`: the same error fails the
  same batch with the same text; a block of one sequence, a long one, is checked from its last 2-bit byte and its
  ambiguity data, `SeqDb::x_nucleotide_check`) and only the lengths of its sequences are kept. The preliminary search reads a
  subject when it makes the subject's plan and drops it once the plan is made (the plans of the short ones are made
  from the bytes of one block of 256 KiB per thread, without a store; a subject whose plan has no residues,
  `LOSAT_X_TBNUNITXLATE`, is read before its plan and held while its units translate it); the composition redo and the traceback of `-comp_based_stats 0` read each subject of their
  lists when they first align it and hold it until the stage has done its last match (`LOSAT_X_TBNSTAGEDROP` below;
  until the stage ends with `=0`); the reports read the sequences they show where a report reads their letters (outfmt 0, the XML family, the tabular rows), and only their lengths otherwise. A long
  sequence is read once while anything holds it: stages and batches that hold it at the same time share it, and a
  thread that asks for it while another thread reads it waits for that read. With more than one thread it is read
  into a buffer that glibc maps (taken from `utils::x_mapped::X_BYTES` when it holds one) and its store is that
  buffer cut to the 2-bit data; at one thread it is read into the heap. While a
  stream is open, the frames of 1 MiB or more of a dropped subject plan go to `X_BYTES` (idle bound
  `LOSAT_X_POOLIDLE`), and the next frames and long reads take them with their pages resident. A read that fails
  after the check fails its query batch with the check's error, after the reports of the batches before. Every
  read packs the stored bytes as the copy did, so every report is the same (the frozen NCBI `-db` fixtures hold the
  reports; `tests/tblastn_db_stream.rs` compares the reports across thread counts and query batches and the errors of a
  database that cannot be read). The shadow mode that compared the check's lengths and error and every read with the
  copy is gone with the copy. Measured 2026-10-24 (workstation, every report identical, against dev 4c9ba700):
  peak RssAnon (untimed, 2 runs) of 30 C. elegans proteins x the C. elegans genome (100 Mb in 7 sequences) 32 /
  87-89 / 139-140 MiB at 1 / 8 / 16 threads against 46 / 109-110 / 153-158 (NCBI 75 / 210 / 212), of 118 fly
  proteins x zebrafish cDNAs (57,776 sequences, `-max_target_seqs 5`) 22-23 / 55-56 / 86-87 against 52 / 81-82 / 112
  (NCBI 16 / 28 / 34), of one fly protein x the human genome (3.1 Gb) 563 against 832 MiB at one thread (NCBI 483).
  At one thread (`perf stat`, 3 runs) instructions:u x0.99990 / x1.00027 / x0.99960 on the three (x1.000002 /
  x1.000095 / x1.000005 with `=0`; on the cDNAs the check reads the database once more) and minor faults 12.7k /
  7.6k / 276k against 19.0k / 14.9k / 344k; cycles:u x0.9959 / x0.9983 / x0.9852. Timed against dev (node 1, warmed
  3 s, 8 alternated pairs): the C. elegans genome x0.991 / x0.986 (24 pairs at 8), the cDNAs x0.999 / x0.999,
  the human genome x0.790 / x0.673 at 8 / 16 threads.
- `LOSAT_X_TBNSTAGEDROP=1` (TBLASTN `-db`, memory; in the strict set from the pull request that added it: on unless
  `=0`): the composition redo and the traceback of `-comp_based_stats 0` of a `-db` stream let a subject's store go once
  the stage has done the subject's last match (redone, or passed by the early termination), instead of holding every
  subject they read until the stage ends (`tblastn/x_db_stream.rs` `XStageSubjects::done`; the number of matches of each
  subject is counted when the stage starts, so no subject is read again). With more than one thread a match redone by a
  worker is done when the worker ends it. It changes no value and no order of computation. NCBI's redo fetches each
  match's subject and releases it after the match (`blast_kappa.c`). Measured 2026-10-25 (workstation, every report
  identical, against the same binary with `=0`): peak RssAnon of one fly protein x the human genome (3.1 Gb) 112-115
  against 533-546 MiB at one thread (NCBI 483; the redo held the stores of 15 chromosomes, 2.3 Gb in 2 bits) and 421
  against 454-490 at 8 threads (697-702 against 691-700 at 16, where the peak is elsewhere); 30 C. elegans proteins x
  the C. elegans genome and 118 fly proteins x zebrafish cDNAs within the spread of their runs (their peaks are not in
  the redo). At one thread (`perf stat`, node 0) instructions:u x1.000000 / x0.999993 / x0.999920 on the three, cycles:u
  x0.9963 / x0.9991 / x0.9990 (8 and 3 rounds beside an A/A arm of `=0`: x0.9964 / x1.0035). Timed (node 1, warmed 3 s,
  alternating pairs): the human genome x0.984 / x0.980 at 8 / 16 threads (8 pairs, none above 1), the C. elegans genome
  x1.000 / x0.985 (16 / 8 pairs), the cDNAs x1.002 / x1.000 (8 pairs; x1.002 is one 10 ms step of 8.17 s).
- `LOSAT_X_BXDIAG=1` (BLASTX, speed; in the strict set since 2026-10-17, by the Owner's rule of that day
  that a speed-up lands on by default (DEFAULT_ON_POLICY_20261017 D4), with the rule below for few hits
  per cell: on by default where the cell test is `LOSAT_X_BXCROSS`'s (an x86-64 CPU with AVX2 and POPCNT;
  `LOSAT_X_BXCROSS=0` or no POPCNT leaves the direct loop on every table), `=0` turns it
  off; set to a value it runs with any cell test (`LOSAT_X_BXCROSS`'s where that runs, otherwise the masked AVX2
  batch on a CPU with AVX2 and the per-hit loop elsewhere); acts
  in the two-hit stage with the standard lookup table,
  `-word_size 3`, when no diagnostic stage traces the seeds)
  the subject scan and the two-hit diagonal test in one loop, with the hits of one lookup cell tested in
  batches of eight (`blastx/x_seed_diag.rs`, exactness argument in its header): blastx calls blastp's
  fused scan of BPSEED with the cell test of `LOSAT_X_BPDIAG` (the AVX2 batch on an x86-64 CPU
  with AVX2, in the form of `LOSAT_X_BXCROSS` where the CPU also has POPCNT, otherwise (with a value of
  `LOSAT_X_BXDIAG`) the masked batch of 2026-11-02's fold, below; the per-hit loop elsewhere), and the
  hits that
  pass the test go to the reference's context
  check, extension, save and entry update unchanged (the context from blastp's `XContextTable`, made once
  per query batch and thread, which gives the binary search's index). It rests on blastx's entries being
  NCBI's `DiagStruct`, as blastp's are. `LOSAT_X_BXDIAG=per-hit` and `LOSAT_X_BXDIAG=scalar` force a cell test
  (diagnostics; the tests use them). Shadow: `LOSAT_X_BXDIAGSHADOW=1` also runs the reference on a copy of
  the diagonal table for every subject and stops at a difference of the hit lists or the tables (the
  tables are not compared for a subject shorter than a word less one letter, whose one scan position the
  reference reads past the sentinel); at exit `[X_SHADOW] BXDIAG subjects_compared=N ...`. With
  `LOSAT_TIMING=1`: `[TIMING] blastx_diag: mode=<avx2-cross-defer|avx2-cross|avx2-defer-masked|scalar|per-hit> subjects=
  cells= batched_hits= candidates= min_hits_per_cell= tables_fused= tables_other= hits_per_cell=LEAST-MOST`
  (until 2026-11-02 it also gave `defer_min_hits_per_cell=` and `tables_avx2=`, `LOSAT_X_BXDEFER`'s rule), and
  with `mode=avx2-cross` `cross_splits=` (`LOSAT_X_BXCROSS`'s windows tested in parts). A query batch whose
  lookup table has fewer hits per used backbone cell than `LOSAT_X_BXDIAGMIN` (in hundredths; built-in 670 for
  every cell test since 2026-11-02, the fold of line SF S3 family bx-seed-loop, which had 670 with
  `LOSAT_X_BXCROSS`'s cell test and 900 with the others; a number set is the cutoff, `0` fuses every batch) is scanned by the direct loop instead, which is faster there; the
  choice is made once per thread and table and changes only the speed (the timing line counts the choices and
  gives the least and most hits per cell, in hundredths). The built-in value was 900 from the pull request
  perf/bx-diag-min (2026-10-23; 500 before, which `LOSAT_X_BXDIAGMIN=500` gives) to the fold of 2026-11-02: with `LOSAT_X_BXHITBF`'s
  branch-free hit test the direct loop is faster up to about 9 hits per used cell. The fused loop with the
  cell test it took then (the AVX2 batch below `LOSAT_X_BXDEFERMIN`, `LOSAT_X_BXDEFER`'s above) against the direct
  loop, user-space cycles, one core (workstation node 1, HT sibling idle, warmed 3 s, 8 rounds with the arm
  order reversed every other round, medians, dev cd8e0f89; every output the same bytes): 4.91 hits per used
  cell (the first 1,500 nt of a 2.1-kb fly cDNA x the fly proteome) x1.198, 5.79 (another cDNA) x1.187,
  6.77 (the whole cDNA) x1.133, then sets of two or three of 50 fly cDNAs: 7.34 x1.073, 8.38 x1.030, 8.63
  x1.006, 8.92 x0.994, 9.37 x0.981, 9.85 x0.975, 10.52 x0.965, 11.46 x0.946, 13.32 x0.899; 50 cDNAs x yeast
  (20-33 per batch) x0.781. Against Swiss-Prot (`-db`, 3 rounds) the loops break even later: 6.77 x1.178,
  8.92 x1.012, 9.81 x1.020, 11.46 x0.968, 13.32 x0.916; 900 lies between the fly proteome's break-even (about
  8.75) and Swiss-Prot's (about 10.7), within 2 % of the faster loop on either side. Against dev's default
  (the same runs; one thread): 5.79 x0.846 cycles (0.655 -> 0.555 s), 6.77 x0.883, 7.34 x0.931, 7.64
  x0.942, 8.38 x0.967, 8.63 x0.997, 8.92 x1.005, the cDNA x Swiss-Prot x0.847 (7.14 -> 6.06 s), 8.92 x
  Swiss-Prot x0.987, and x0.993-1.004 on the tables of 9 or more (the same loop as before).
  At 8 threads (node 1, 8 physical cores, 5 rounds): 6.77 x0.885 (0.125 -> 0.116 s), 8.38 x0.964, 8.92
  x1.007, 9.85 x0.999. With `LOSAT_X_BXCROSS`'s cell test (the default since 2026-10-23 on CPUs with
  POPCNT) the two loops break even lower, and the built-in value for it was 670 from the pull request
  perf/bx-diag-min-cross (2026-10-24; `LOSAT_X_BXDIAGMIN=900` gave the rule before; `LOSAT_X_BXCROSS=0` and
  CPUs without POPCNT kept 900 until the fold of 2026-11-02, which made 670 the one value). The fused loop with the cross cell test against the direct loop (the same
  setup, dev 0030a116, load 8-33; one cDNA or two of the 50 x the fly proteome, 8 rounds): 4.05 x1.135, 4.73
  x1.120, 5.34 x1.063, 5.79 x1.050, 6.33 x1.011, 6.63 x0.986, 6.77 x0.978, 7.34 x0.947, 7.64 x0.928, 8.38
  x0.902; against Swiss-Prot (4 rounds): 4.73 x1.157, 5.34 x1.109, 5.79 x1.079, 6.33 x1.040, 6.77 x1.013,
  7.64 x0.951. They break even at about 6.46 against the fly proteome and 6.95 against Swiss-Prot; 670 lies
  between, within 2 % of the faster loop on either side. Against dev's default (0030a116, the same setup, load
  7-28; every output the same bytes), 1 thread, cycles: 6.77 x0.979 (0.578 -> 0.565 s), 7.34 x0.949, 7.64
  x0.927, 8.38 x0.898 (0.738 -> 0.663 s), 8.63 x0.885 (0.886 -> 0.784 s); against Swiss-Prot 6.77 x1.013, 7.64
  x0.953, 8.38 x0.920 (7.17 -> 6.61 s); at 8 threads (node 1, 8 physical cores, 5 rounds) 6.77 x0.994, 7.34
  x0.960, 8.38 x0.911. Tables that take the same loop as before: 6.33 x0.999 and 6.63 x1.001 (the direct loop),
  9.81 x0.999 and 50 cDNAs x yeast x0.991 (fused; one node-0 core), the same instructions (x1.0000). The
  measurements below came before `LOSAT_X_BXHITBF`. Where the two loops break even (one thread, `-db`
  UniRef50 2026_03 volume 00, both with `LOSAT_X_BXSCANBF=1`, every batch fused against the direct loop, 2
  rounds each; pieces of one 2.1-kb fly cDNA): hits per used cell 1.76 (300 nt) x1.20 the direct loop's
  time, 2.34 (600 nt; 41.48-41.66 s against 34.49-34.74 s) x1.20, 3.56 (1,000 nt; 68.62-68.86 / 63.55-64.06)
  x1.08, 4.91 (1,500 nt; 90.08-90.69 / 90.53-90.54) x1.00, 6.77 (2.1 kb) x0.94; 800 cDNAs x yeast (4.52-68.09
  per batch) x0.72. With the rule (dev 7a5a7f87 and the rule, `LOSAT_X_BXDIAG=1` against the direct loop, both
  with `LOSAT_X_BXSCANBF=1`, 2 rounds): 300-nt case 21.28 / 21.42 s against 21.20 / 21.68 s (no batch fused;
  x1.20 before the rule), 2.1-kb cDNA 115.83 / 116.79 s against 124.72 / 125.10 s (x0.93), 800 cDNAs x yeast
  62.59 / 62.65 s against 87.55 / 88.18 s (x0.71; 185 of the 186 tables fused at 2 threads, the other one had
  4.52 hits per cell). As the default (the same pairs, the default against `LOSAT_X_BXDIAG=0`, the build before the
  short-subject fix): 300-nt case x1.00, 2.1-kb cDNA x0.93, 800 cDNAs x0.71-0.73, max RSS equal; at 8 threads
  (node 1 CPUs 56-63, rotated, 3 rounds, with `BXFLOW`) 800 cDNAs 7.93-8.18 s against 11.05-11.22 s (0.71-0.73), 50
  cDNAs 0.430-0.435 s against 0.582-0.585 s (0.74), 5 cDNAs 0.086-0.089 s against 0.097-0.098 s (0.89-0.91). Measured 2026-10-17 on the workstation (node 1, one pinned core per arm with its HT sibling
  idle, the two arms of a pair at once on two cores, swapped every round; 800 fly 2L cDNAs x yeast proteome
  `-db`, `-evalue 1e-5`, one thread; every output the same bytes as NCBI BLAST+ 2.17.0): off 87.21-87.42 s,
  on 64.83-65.41 s (x0.742-0.750; user-space cycles x0.743-0.750, instructions 723 -> 493 G; 4 rounds);
  NCBI `-db` 102.86 / 103.61 s beside on 66.55 / 64.99 s (2 rounds). The context table came after these
  rounds: against the same loop with the binary search, 63.75-64.43 s against 64.84-65.27 s (x0.977-0.994,
  cycles x0.98-0.99; 4 rounds). The shadow compared 1,214,400 subject scans with `-db` and with `-subject`
  (6 threads): identical. At 16 and 8 threads (same input, node 1 physical CPUs 32-47, 5 rounds rotated, with
  `LOSAT_X_BXFLOW=1`): 4.28 s against 5.80 s without it and NCBI's 7.06 s at 16 threads (0.61 of NCBI), 8.25 s
  against 11.33 s and 13.19 s at 8 threads (0.63). Where a lookup cell has few hits it is slower than the
  default loop: a 300-nt cDNA piece against the first UniRef50 volume (1.5 hits per cell) took 38.08-38.25 s
  against 35.76-36.00 s (one thread, 3 rounds, the build with LOSAT_X_BXSCANBF's const parameter; instructions
  250 G against 185 G); one 2.1-kb cDNA (4.6 hits per cell) 123.9-125.4 s against 140.4-140.7 s.
- `LOSAT_X_BXDEFER` and `LOSAT_X_BXDEFERMIN` (BLASTX; folded on 2026-11-02, line SF S3 family bx-seed-loop: both
  variables are gone, the timing line lost `defer_min_hits_per_cell=` and `tables_avx2=`, and
  `LOSAT_X_BXDIAGMIN` has one built-in value, 670, for every cell test; default results and paths did not
  change). Why: since the built-in `LOSAT_X_BXDIAGMIN` became 900 (and 670 with `LOSAT_X_BXCROSS`'s cell test) the
  places this cell test acted fused tables from 9 hits per used cell, above its own rule's 7, so no default table
  reached the rule, and every CPU with AVX2 has POPCNT, so `LOSAT_X_BXCROSS`'s cell test took every default
  table. The masked AVX2 batch (`XDiagMode::Avx2DeferMasked`, `mode=avx2-defer-masked`) stays: with a value of
  `LOSAT_X_BXDIAG` and `LOSAT_X_BXCROSS=0` (or no POPCNT) the fused loop takes it on every table it fuses; with
  `LOSAT_X_BXDIAG` not set and `LOSAT_X_BXCROSS=0` the direct loop scans every table. The plain AVX2 batch
  (`mode=avx2`, what `LOSAT_X_BXDEFER=0` gave) is no longer reachable from blastx (blastp has it). History of the
  switch (in the strict set from the pull request that added it, 2026-10-22; acted only where
  `LOSAT_X_BXDIAG`'s fused loop ran the AVX2 batch: not with `LOSAT_X_BXDIAG=per-hit` or `scalar`, not in the
  query batches `LOSAT_X_BXDIAGMIN` sent to the direct loop, not without AVX2 (aarch64, WebAssembly); since
  2026-10-23 only where `LOSAT_X_BXCROSS`'s cell test did not run (`LOSAT_X_BXCROSS=0`, or an x86-64 CPU with AVX2
  but without POPCNT); and only on lookup tables with at least `LOSAT_X_BXDEFERMIN` hits per used backbone cell,
  in hundredths, built-in 700, `0` every table) the AVX2 cell test of blastx's fused seed loop took
  `LOSAT_X_BPDEFER`'s form (a cell's candidates packed into a buffer and the closure called for the buffer after
  the cell's batches) with two changes
  (`blastp/x_seed_fused.rs`, `x_cell_avx2_defer_masked`, exactness argument in its header): a batch of eight
  computes its lanes' entry coordinates from scalar loads of its hits instead of taking them out of the hits'
  register, and a cell's last batch of fewer than eight hits is one masked AVX2 batch (the missing lanes repeat
  the first hit, their candidate bits are dropped) whose candidates join the buffer, instead of the scalar batch
  after the calls. Exact because the coordinates are the same values, a repeated lane reads its entry before any
  lane of the batch stores and so stores what lane 0 stores, and the last batch's hits lie on entries no other
  hit or call of the cell touches, so the calls (the batches', then the last batch's lanes', each with the entry
  its test read) and the table are those of the AVX2 batch. The two changes come from line R's probes of
  blastp's cell test (v4-R R2); a tested blastx cell holds about 19 hits on the dense tables below (a blastp
  cell about 106), so its last batch weighs more than in blastp. On a table with few hits per cell the masked
  batch and the buffer cost more than they save, so a fused table below `LOSAT_X_BXDEFERMIN` took the AVX2
  batch (`XDiagMode::Avx2`); the choice was made with `LOSAT_X_BXDIAGMIN`'s, once per thread and table. Since
  `LOSAT_X_BXDIAGMIN`'s built-in value became 900 (above 700; perf/bx-diag-min; 900 still where `LOSAT_X_BXCROSS`'s
  cell test did not run, 670 where it did) every fused table that took
  `LOSAT_X_BXDEFER`'s choice took its masked batch, and the rule acted where `LOSAT_X_BXDIAGMIN` was set below 700. With the
  fused loop and the cell test forced on every table (`LOSAT_X_BXDIAGMIN=0` on the binary before the rule,
  the same as `LOSAT_X_BXDEFERMIN=0`; one core, medians of 3-8 per-round ratios against `=0`) the cycles were
  x1.256 at 1.76 hits per used cell (a 300-nt
  cDNA piece x UniRef50 volume 00), x1.071 at 3.56 (the first 1,000 nt of a 2.1-kb fly cDNA x the fly proteome),
  x1.033 at 4.91 (its first 1,500 nt), x0.996 at 6.77 (the whole cDNA), x0.946 at 9.85 and x0.918 at 13.32 (the
  first two and three of 50 fly cDNAs x the fly proteome), and x0.88-0.96 on the tables of 17 or more below.
  `LOSAT_X_BPDEFER`'s form alone was x1.040 at 6.77 and x0.984 on 201 reads x Swiss-Prot (v4-A14 probe 1, ratios
  of the medians of 3 rounds). blastp's hot loops keep their code (the new cell test has a scan function of its
  own, never inlined, and calls none of blastp's batch helpers): blastp's instructions equal those of the binary
  before the switch (x1.0000 on 300 yeast proteins x the fly proteome and on WSSV x PajaWSV; since 2026-10-24
  blastp's fused loop takes the same cell test by default, `LOSAT_X_BPMASKED`). With `LOSAT_TIMING=1` the
  `blastx_diag` line named `mode=avx2-defer-masked` (`mode=avx2` with `=0`) and gave
  `defer_min_hits_per_cell=` and `tables_avx2=` (the fused tables that took the AVX2 batch). Shadow:
  `LOSAT_X_BXDIAGSHADOW=1` compares the cell test each table takes with the reference for every subject (every
  table through the fused loop); `tests/blastx_diag.rs` runs the masked batch (a value of `LOSAT_X_BXDIAG` with
  `LOSAT_X_BXCROSS=0`) with every search and the shadow; its cases of the switch and the rule's two ends went
  with the fold.
  Measured (history) 2026-10-22/23 on the workstation (perf stat, user-space cycles; one physical core of node 0 for 1 thread
  (load 29-83 on the other cores), 16 physical cores of node 1 for 8 and 16 threads; HT siblings idle, the booked
  cores warmed 3 s before every run, arm order reversed every other round, one warm-up round and 8 rounds (fewer
  where noted), medians of the per-round ratios; the arms the binary before the switch (dev 7cfec941), the
  default and `=0`, before the rule (every fused table of these cases has at least 7 hits per used cell: with
  the rule `tables_avx2=0` on all of them); every report identical between the arms, and those of the reads
  the same bytes as NCBI BLAST+
  2.17.0 on the same inputs (201 reads `-num_threads 1`; 1,000 reads the 16-thread runs of line E); max RSS
  equal within noise), the default against `=0`, wall time (cycles, instructions, mispredicted branches): at 1
  thread 201 150-nt reads (ENA DRR1003646) x Swiss-Prot (`-db`) x0.891 (59.19 -> 52.26 s; x0.881, x0.933,
  x0.730; 5 rounds), a 50-kb piece of the fly genome x the fly proteome (`-db`, `-evalue 1e-5`) x0.894 (11.44 ->
  10.23 s; x0.892, x0.948, x0.726), 800 cDNAs x the yeast proteome x0.899 (64.27 -> 56.53 s; cycles x0.884; 3
  rounds), 50 fly cDNAs x the yeast proteome x0.885, AP027131 x AP027133 proteins (`-subject`,
  `-query_gencode 4`) x0.964 (cycles x0.959); at 16 threads 1,000 reads x Swiss-Prot x0.887 (17.24 -> 15.28 s),
  1,000 reads x UniRef50 2026_03 volume 03 x0.882 against the binary before the switch (256.3 -> 226.0 s; 2
  rounds), the 50-kb piece x0.905; AP027131 x AP027133 at 8 threads x0.969. NCBI BLAST+ 2.17.0 on the same
  machine (v4-E): 201 reads x Swiss-Prot 79.6 s at 1 thread, 1,000 reads x Swiss-Prot 27.59 s and x UniRef50
  volume 03 385.94 s at 16 threads. The `=0` arm has the instructions of the binary before the switch (x1.000)
  and cycles x0.998-1.013 (the loaded node's spread). Query batches the direct loop takes by default (the
  first 1,000 nt of the 2.1-kb cDNA, 3.56 hits per used cell; the 300-nt piece, 1.76, no hit at `-evalue 1e-5`)
  run as before (instructions x1.000, cycles x1.001-1.003). With the rule (the same setup on one core, load
  9-38, 3-8 rounds, against the binary before the switch): the 2.1-kb cDNA (6.77 hits per used cell, the AVX2
  batch by the rule) cycles x1.001, instructions x1.000; cycles x0.939 and x0.922 for the first two and three
  of the 50 cDNAs, x0.963 for AP027131 x AP027133, x0.897 for the 50-kb piece and x0.883 for 201 reads x
  Swiss-Prot (wall 57.33 -> 50.10 s, x0.874); WSSV x PajaWSV (blastp) instructions x1.000.
- `LOSAT_X_BXCROSS=1` (BLASTX, speed; in the strict set since 2026-10-23, from the pull request that made it the
  default (it was an opt-in prototype the same day; v4-A14 DECISIONS 14): `=0` turns it off; acts where
  `LOSAT_X_BXDIAG`'s fused loop runs an AVX2 cell test on an x86-64 CPU with AVX2 and POPCNT, on every table
  that loop takes; off, or without POPCNT, `LOSAT_X_BXDIAG` not set leaves the direct loop on every table, and
  with a value the fused loop takes the masked AVX2 batch (`mode=avx2-defer-masked`; until 2026-11-02 the default
  took `LOSAT_X_BXDEFER`'s choice there).
  Which query batches the fused loop takes is `LOSAT_X_BXDIAGMIN`'s rule, a separate cutoff against the direct
  loop, whose built-in value is 670 (900 without this cell test until 2026-11-02)) the cell test of blastx's fused seed loop cuts its batches of eight
  across the cells' borders (`blastp/x_seed_fused.rs`, `x_scan_cells_avx2_cross`, exactness argument in its
  header): the scan's call for a cell only stages the cell's hits (coordinate and subject offset, in two arrays
  of one allocation; three masked 8-lane loads per 24 hits), and every 1,024 staged hits and at the end of the
  subject the stage is tested in windows of eight (the lanes of the masked batch, once `LOSAT_X_BXDEFER`'s, with a subject offset
  per lane, each batch's candidates called before the next batch loads). A window on eight distinct diagonal
  entries is one batch; otherwise (out of line) it is tested in parts, each the hits before the first one that
  repeats an entry of the part (the batch's other lanes repeat the part's first hit and their candidate bits are
  dropped); the last fewer than eight hits of a subject take the hit-by-hit test, and those left by a chunk move
  to the front of the stage. The query offset is not staged: it is `(coord + s) & mask`, exact because the
  diagonal table has at least the query's length plus the window entries (`Diagonals::new`), so the mask is at
  least every query offset of the lookup (asserted once per subject scan). Exact because the stream is the scan
  order and a hit's test and call touch only its own entry. blastp's loops keep their code: the cross-cell test
  is a cell test of blastx's own (`XCellTest::Cross`), which blastx's entry runs, not one of the choice that
  blastp's loops share, and every function of it is its own (in the release build every blastp function
  disassembles as before but for the `x_reach!` indices of the switches after this one in the set). With
  `LOSAT_TIMING=1` the `blastx_diag` line names `mode=avx2-cross` and ends with
  `cross_splits=N` (the windows with a repeated entry, counted per subject scan with the other counts).
  Shadow: `LOSAT_X_BXDIAGSHADOW=1` compares it with the reference for every subject. Tests: `tests/blastx_diag.rs`
  runs the default with every search, the shadow and `=0` (the direct loop with `LOSAT_X_BXDIAG` not set, the masked batch with a value); the unit tests compare the program's chunk of 1,024
  and a test-only instance of 64 (drains inside a subject and inside long cells, the hits a drain leaves,
  windows in parts across drains) with the per-hit test. Measured 2026-10-23 on the
  workstation (`perf stat` user-space cycles, one physical core of node 0 with its HT sibling idle, warmed 3 s
  before every run, load 8-38 on the other cores; paired arms, arm order reversed every other round; medians of
  the per-round ratios, 3 rounds for the reads and 5 for the others), against `dev`'s default 7d7de94e
  (then `LOSAT_X_BXDEFER` with its rule, `LOSAT_X_BXUNGAPSENT`, `LOSAT_X_BXDIAGMIN` 900): 201 reads x Swiss-Prot
  (`-db`; 26-28 hits per used cell) x0.880, a 50-kb piece of the fly genome x the fly proteome (`-db`,
  `-evalue 1e-5`; 25-31) x0.899, AP027131 x AP027133 proteins (`-subject`, `-query_gencode 4`; 17-26) x0.966,
  50 fly cDNAs x the yeast proteome (20-33) x0.861, the first two of the 50 cDNAs x the fly proteome (9.85)
  x0.830; one 2.1-kb fly cDNA x the fly proteome (6.77) takes the direct loop in both arms (x1.002); `=0`
  x0.998-1.001 of `dev` (the same instructions); every report identical between the arms, max RSS unchanged.
  Against the `dev` before `LOSAT_X_BXDIAGMIN` 900 and `LOSAT_X_BXUNGAPSENT` (32a63700, the same setting),
  where every one of these tables was fused: x0.886, x0.916, x0.973, x0.882, the one cDNA x0.846, the two
  x0.847. `LOSAT_X_BXDIAGMIN`'s built-in 900 (perf/bx-diag-min; gone with the fold of 2026-11-02) was set
  against `LOSAT_X_BXDEFER`'s cell test, so in these runs a table below 9 hits per used cell, such as the one cDNA's (6.77), still took the
  direct loop; with this cell test the fused loop breaks even with the direct loop lower, and the pull
  request perf/bx-diag-min-cross (2026-10-24, v4-A10 bdm5) set the built-in value to 670 for it (the
  `LOSAT_X_BXDIAG` item). Timings: v4-A14 results/fab2
  (and results/fab against 32a63700). Design and timings: v4-A14 results/p2 (DESIGN_cross.md), results/p3
  (DESIGN_s3.md: round 2) and results/p4 (DESIGN_s4.md: round 3, r1's loops with a constant chunk, two arrays,
  the windows with a repeated entry out of line).
- `LOSAT_X_BXCROSSDEFER=1` (BLASTX, speed; in the strict set since 2026-10-25, from the pull request that
  added it (perf/bx-cross-defer; v4-A14 DECISIONS 26-27): `=0` turns it off; acts where `LOSAT_X_BXCROSS`'s cell
  test runs, on every table blastx's fused loop takes on an x86-64 CPU with AVX2 and POPCNT) `LOSAT_X_BXCROSS`'s
  drain defers the candidate calls of each unit (a batch, or a part of a window with a repeated entry) to the
  unit that next reads one of their entries, or to the end of the drain (`blastp/x_seed_fused.rs`,
  `x_cross_drain` with `DEFER`, exactness argument in its header). A unit's candidate lanes store a marker
  (raw `0xFFFFFFFF`: flag 1, last_hit -1) in place of the entry they read, and the unit's entries, first
  staged hit and candidate lanes go to a record past the stage (written without a branch). A unit that reads the
  marker stops before it stores, the deferred calls are made in order (each entry first given back the value
  its unit read) and the unit runs again from the loop's head, so no value of the batch is live across the
  calls; the calls left are made after the drain's last window, before the last hits are tested or moved. A
  stored value equal to the marker (a flagged last_hit of -1, which a call writes only near offsets of 2^31)
  met with no call deferred is taken as a value: the unit runs with its calls at once, so any table gives the
  same calls. Why it pays: the candidate branch of each batch was the drain's main misprediction (28 of 127
  issued uops per batch squashed, v4-A14 s7); on dense tables a later unit reads a deferred candidate's entry
  0.017-0.029 times per batch and 5.5-7.2 candidates are deferred between two such reads (v4-A14 s8 probe,
  results/s8c). With `LOSAT_TIMING=1` the `blastx_diag` line names `mode=avx2-cross-defer` (`=0`:
  `mode=avx2-cross`). Tests: `tests/blastx_diag.rs` compares the default's and `=0`'s reports with the
  reference's (every search, the built-in `LOSAT_X_BXDIAGMIN` and every table fused) and checks the mode the
  timing line names; the unit tests compare both arms with the per-hit test (the chunks of 1,024 and 64,
  random and edge tables, which hold the marker's value) and count the units that met the marker with calls
  deferred and with none. Measured 2026-10-25 on the workstation (`perf stat` user-space cycles, one physical
  core of node 0 with its HT sibling idle, warmed 3 s before every run, load 13-109 on the other cores; paired
  arms, arm order reversed every other round, medians of the per-round ratios), against `dev` f24d161f: a 50-kb
  piece of the fly genome x the fly proteome (`-db`, `-evalue 1e-5`) x0.966 (0.956-0.970, 8 rounds), 201 reads x
  Swiss-Prot (`-db`) x0.968 (0.962-0.974, 6), 50 fly cDNAs x the yeast proteome x0.981 (0.973-0.991, 10), AP027131
  x AP027133 proteins (`-subject`, `-query_gencode 4`) x0.992 (0.984-0.999, 10), one 2.1-kb fly cDNA x the fly
  proteome (6.77 hits per used cell, fused by `LOSAT_X_BXDIAGMIN` 670) x0.985 (3), and 300 UniRef50 reads (the
  direct loop) x1.002 with the same instructions; instructions x1.003-1.009 and branch misses x0.83-0.94. `=0`
  x0.999-1.004 of `dev` with the same instructions (x1.000): its drain, window and scan functions disassemble as
  `dev`'s but for rip-relative displacements. Every report and every `blastx_diag` count identical between the
  arms, max RSS unchanged; the record space is 41 KiB per thread (the stage's 8 KiB without it). The prototype
  (v4-A14 s8ab2, against `dev` 4c9ba700): x0.961, x0.979, x0.983, x0.994, x0.982. Design and timings: v4-A14
  results/s8c, s8ab1 (a first form with the calls inside the loop: spills, x0.99-1.02), s8ab2 and s9ab1.
- BXSCANBF (graduated 2026-10-26; no variable; always on) (BLASTX, speed; in the strict set from 2026-10-17, by
  the Owner's decision of that day (perf-si-20261012 DECISIONS 175) on the timings below; the scan on the
  presence test went with the variable; acts where the standard lookup table is scanned by the
  direct loop or `LOSAT_X_BXDIAG`'s fused loop; the compressed table and the reference word finder keep their
  scans) the subject words with lookup hits found without a branch per word, by tblastn's loop of
  tblastn's TBNSCANBF (graduated; `tblastn/x_scan_bf.rs`): the direct loop calls `x_scan_one_range_bf` (NCBI's pairs,
  stop and resume, and `scan_range` values, as for tblastn), and the fused loop lists the positions of 256 words
  with `x_block_candidates` and passes the listed cells to the cell test in position order
  (`blastp/x_seed_fused.rs` `x_scan_cells_bf`; blastp's own loop keeps the presence test). Tests: the fused
  scan's cell calls against the presence-test loop (word sizes 2-4, subjects across the block edges, letters
  1-27), blastp's batched-scan test, and the blastx direct and fused tests (the direct loop's arms, the fused loop's cell tests). The
  shadow switches check it (`LOSAT_X_BXDIAGSHADOW=1` compares the fused loop against `word_finder`; the
  direct loop's shadow, `LOSAT_X_BXSCANSHADOW`, went with BXSCAN): 1,214,400 subject scans (800 fly cDNAs x yeast) and 2,178,834
  (a 300-nt piece x the first UniRef50 volume), identical. Measured 2026-10-17 on the workstation (node 1, one
  pinned core per arm with its HT sibling idle, the two arms of a pair at once, swapped every round; `-db`,
  `-evalue 1e-5`, one thread; every output the same bytes in both arms, and NCBI BLAST+ 2.17.0's where timed
  beside), time with the switch / without, the share of subject words with hits in brackets: a 300-nt fly cDNA
  piece x UniRef50 2026_03 volume 00 (2,178,834 sequences, 3.0 G residues) [34 %] default loop 21.73-22.01 s /
  35.63-36.41 s (x0.60-0.61; branch misses 1.84 G -> 0.46 G, cycles 108.9 G -> 64.3 G; 4 rounds), with BXDIAG
  26.08-26.32 / 38.09-38.61 (x0.68), NCBI `-db` 32.20 / 32.53 s beside 21.77 / 22.12 s (0.68); one 2.1-kb cDNA x
  the same volume [84 %] 125.64 / 126.05 against 140.41 / 140.72 (x0.89-0.90), with BXDIAG 117.68 / 118.63
  against 123.88 / 125.37 (x0.95), NCBI 145.69 / 146.47 beside 125.45 / 125.52 (0.86); 800 fly cDNAs x yeast
  proteome [98 %] 87.70-87.90 / 87.95-88.31 (x0.997), with BXDIAG 63.22 / 63.44 against 63.77 / 63.84 (x0.99).
  The switch is a const parameter of the fused loop's instances (`x_scan_two_hit_in_with::<BRANCH_FREE>`): with a
  run-time flag, blastx's instance had two scan loops calling the cell test, which the compiler then kept as a
  function called per cell (BXDIAG 41.3-42.3 s on the 300-nt case). The timings above are of the const form;
  BXDIAG without this switch then takes 0.96 of the time before this change on the 300-nt case (37.70-38.10 s
  against 39.20-40.01 s) and 0.99 on 800 cDNAs (63.04 / 63.05 s against 63.45 / 63.72 s). Timed again for the
  default set on the same day (dev 7a5a7f87, whose candidate loop has tblastn's constant shifts of #273; one build,
  the default against `LOSAT_X_BXSCANBF=0`, the pairs as above): 300-nt case 21.34-21.36 s against 35.46-35.54 s
  (x0.60, 3 rounds; it reports no hit at `-evalue 1e-5`, and at the default `-evalue` both give the same rows),
  NCBI `-db` 32.01 / 32.34 s beside 21.10 / 21.21 s (0.66); one 2.1-kb cDNA 125.24 / 129.20 s against 140.98 /
  144.30 s (x0.89), NCBI 146.17 / 145.90 s beside 124.73 / 124.64 s (0.85); 800 cDNAs x yeast 87.18 / 87.26 s against
  87.76 / 88.12 s (x0.99-1.00).
- `LOSAT_X_BXSCANBYTE=1` (BLASTX, speed; in the strict set from the pull request that added it, 2026-10-24;
  `=0` turns it off; `blastx/x_seed_scan.rs` module header, the scan is tblastn's `x_scan_one_range_byte` of
  `LOSAT_X_TBNSCANBYTE`, `tblastn/x_scan_byte.rs`, exactness argument in its header) the direct seed loop scans the
  standard lookup table on presence bytes: the table indexes of a block of 256 words made from the letters at once
  (16 at a time with AVX2), the candidates listed from one byte per index (32 KB, built the first time the direct
  loop scans with the table and kept with it, in `Lookup::Standard`), the hits of a block copied with one capacity
  test, NCBI's test where the block could fill the offset array; every call returns NCBI's pairs and scan range
  (unit test `x_scan_byte_gives_ncbi_scan_calls_on_blastx_tables` against `s_blast_aa_scan_subject_one_range` on
  `Lookup::new` tables, both forms; the direct-loop test with its arms). The offset array has four entries more.
  CPU policy of `LOSAT_X_TBNSCANBYTE`: the AVX2 forms where the CPU has AVX2, the plain loops for sparse tables (at
  most three tenths of the backbone cells with hits) on x86-64 without AVX2 and in WebAssembly, nothing on aarch64
  (not timed); `=2` the plain loops for every table, `=3` for sparse tables, on any CPU. Elsewhere the scans before
  run (the branch-free scan of BXSCANBF, always on). The fused loop of `LOSAT_X_BXDIAG` keeps its candidate loop, so on x86-64
  with AVX2 the switch acts on batches with few hits per used cell (`LOSAT_X_BXDIAGMIN`: short queries, one to a
  few; with `LOSAT_X_BXCROSS=0` on every batch); without AVX2 the direct loop scans every batch, and the switch acts on its sparse tables on x86-64 and in
  WebAssembly (on aarch64 not at all). `LOSAT_X_REACH` records per
  subject of the direct loop whether the scan on presence bytes ran (`tests/blastx_scan_byte.rs`: reports equal
  `=0`'s for a 300-nt query (sparse table) and six pieces of 1.5-3.9 kb (not sparse), `-subject` and `-db`, 1 and 3
  threads). Timed on the workstation (node 1, one booked
  core with its HT sibling idle (8 for 8 threads), warmed 3 s, perf stat; dev 4c9ba700 against this change, the same
  output in every arm; medians of the per-round ratios): one 300-nt read x UniRef50 2026_03 volume 00 `-db`
  (2,178,834 sequences), one thread x0.828 (14.96 -> 12.40 s; cycles x0.815, instructions x0.741), `-evalue 1e-5`
  x0.833, `=3` (the plain loops) x0.853, 8 threads x0.871 (2.53 -> 2.20 s); batches the fused loop takes keep their
  time: one 2.1-kb cDNA x fly proteins `-db` x1.002 (instructions x1.000), 200 and 1,000 short reads x SwissProt
  `-db` x0.998 (one thread) and x0.996 (8 threads); the policy without AVX2 (`=3` with `LOSAT_X_BXDIAG=0` against dev
  with `LOSAT_X_BXDIAG=0`) x1.004 on the cDNA (its table is not sparse; rounds x0.984-1.019) and x0.994 on 200 reads.
  Again on dev bf92a834 (the pull request's base): the read x0.822 at one thread and x0.869 at 8, the cDNA x0.999, 200
  reads x0.983. In WebAssembly (the serial WASI command under Node 26, one booked core, the default and `=0`
  alternating in one process, 12 rounds, process CPU): the 300-nt read x the first 3,000 UniRef50 records (43.7 M
  letters) `-db` x0.882 (`-subject` x0.981: there the FASTA parse dominates); 600-nt, 1-kb and 2.1-kb queries
  x0.997-1.001 (their tables are not sparse, the scan before runs); 50 cDNAs x yeast proteins x1.002, as much as an A/A
  arm (`=1`) against the default (arm order). With `=2` the plain loops take x0.963 on the 600-nt query's table and
  x1.026 on the 2.1-kb one there, so blastx keeps tblastn's boundary of three tenths (a boundary of its own is not
  measured yet).
  Probe: `exp/a7-probe` bf65e9ba (v4-A s7, `DESIGN_A7.md`).
- `LOSAT_X_BXHITBF=1` (BLASTX, speed; in the strict set from the pull request that added it, 2026-10-20; `=0`
  turns it off) the direct seed loop (`blastx/x_seed_scan.rs`, exactness argument in its
  header) takes the two outcomes most hits of a flag-0 cell end in, too far (`diff >= window`: `last_hit = s +
  diag_offset`) and overlap (`diff < wordsize`: nothing), without a branch, as TBNHITBF (graduated) does for
  tblastn: with the packed cell `raw`, `d2 = 2 (s + diag_offset) - (raw << 1)` is twice the reference's `diff`
  (the shift brings the 31-bit field's sign; a cleared cell holds `-window`), a candidate is `(d2 - 2 wordsize)
  as u32 < 2 (window - wordsize)`, and the cell's next value is selected. A flag-1 cell or a candidate goes on
  to the reference's statements. Acts per subject where `1 <= wordsize < window`, `0 <= diag_offset` and
  `diag_offset + length + 2 window < 2^30` (otherwise the loop without the test runs: the one-hit search,
  and offsets near 2^30 at a `-window_size` near 2^29). The direct loop scans every subject where the CPU
  has no AVX2 and POPCNT (aarch64, WebAssembly: `LOSAT_X_BXDIAG` is off there) or with `LOSAT_X_BXCROSS=0`, and on x86-64 the query batches whose
  lookup table has fewer hits per used cell than `LOSAT_X_BXDIAGMIN` (short queries); the fused loop of
  `LOSAT_X_BXDIAG` keeps its own tests. With `LOSAT_X_BXDIAG=0` 79-81 % of the hits of the three cases
  below are too far, 14.5-17.4 % overlap, 1.2-1.4 % meet a flag-1 cell and 2.1-2.8 % are candidates
  (v4-A10 probe bx1). The shadow (`LOSAT_X_BXSCANSHADOW=1`, which also printed `[X_SHADOW] BXHITBF
  subjects_compared=N`) went with BXSCAN's graduation; `LOSAT/tests/blastx_hitbf.rs` compares whole reports.
  Measured 2026-10-20 on the workstation (node 0: one physical core for 1 thread, 8 physical cores for 8
  threads, HT siblings idle, the booked cores warmed 3 s before every run, arm order reversed every other
  round, one warm-up round + 6 rounds (4 for the last two cases), medians of the per-round ratios against
  dev 7e73ec9c, the binary before the switch; every report identical between the arms, max RSS equal; a gate
  ran on other node-0 cores meanwhile). With `LOSAT_X_BXDIAG=0` (the direct loop on every subject, as
  without AVX2): AP027131 x AP027133 proteins (`-subject`, `-query_gencode 4`) 1 thread x0.912 (4.62 ->
  4.22 s) and 8 threads x0.949; 50 fly cDNAs x the yeast proteome (`-db`, `-evalue 1e-5`) x0.812 (4.06 ->
  3.30 s) and x0.824; one 2.1-kb fly cDNA x the fly proteome (`-db`, `-evalue 1e-5`) x0.812 and x0.886.
  By default with AVX2 the fused loop takes these three (x1.00-1.01, instructions x1.000). Where the direct
  loop takes a table by default: a 300-nt fly cDNA piece x UniRef50 2026_03 volume 00 (`-db`, 1.76 hits per
  used cell; no row at `-evalue 1e-5`) x0.839 (20.19 -> 17.19 s) and x0.887 at 8 threads (4.18 -> 3.73 s);
  the first 1,000 nt of the 2.1-kb cDNA x the fly proteome x0.833 and x0.916. Cycles x0.81-0.92,
  instructions x0.90-0.96, mispredicted branches x0.38-0.73. With `=0` the instruction counts equal those of
  the binary before the switch (x1.000), cycles x0.997-1.005. On WebAssembly (the serial WASI command,
  wasm32-wasip1, under Node 26 in one process: the module compiled once, a new instance per search, the
  default and `=0` alternating, one warm-up round + 8 rounds; the development machine's i9 at load 16-59, so
  the spread is wide): AP027131 x AP027133 x0.907 (6.11 -> 5.20 s, 7 of 8 rounds below 1), 50 cDNAs x the
  yeast proteome (`-subject`) x0.848 (CPU time x0.789, 8 of 8), the same output bytes. (blastx is not in
  LOSAT Web's ABI v2 yet.)
- `LOSAT_X_BXUNGAPSENT=1` (BLASTX, speed; in the strict set from the pull request that added it, 2026-10-23;
  `=0` turns it off) the two-hit ungapped extension of BLASTX's seed loops, the fused loop of `LOSAT_X_BXDIAG`
  and the direct loop (BXSCAN, graduated), for word size 3, with the walks of blastp's `LOSAT_X_BPUNGAPSENT`
  (`blastx/x_ungap_sent.rs`, exactness argument in its header; the extension is
  `blastp/x_ungapped_sent.rs`'s `two_hit_in`): no bound test in the walks, the extension inlined into a
  closure (fused loop) or a loop instance (direct loop) of its own. The six-frame query buffer has a sentinel
  before each frame and after the last, and each subject is read with a sentinel at each end (a database's
  bytes as they are stored, the sequence before ending with one; a FASTA subject encoded with one added
  before it), so a step on a sentinel (-128 in the byte table) stops each walk where the reference's
  sentinel or bound stops it (the left walk at the bound one step later, on the buffer's first byte), with
  the same best, length and stop. The checks of BPUNGAP's walks are made once: the word size, the window,
  residue codes below 32 and the sentinels per query batch, codes and sentinels per subject; the X-drop
  (1..=127) per candidate (the X-drops are the contexts' parameters). Off the path (the reference
  extension): a word size other than 3 (`-word_size 5`), `-window_size 0` (one-hit), a query batch of more
  than `i32::MAX / 16` letters, a candidate whose context's X-drop is outside 1..=127, and the searches
  where `word_finder` runs (the diagnostic stages that trace the seeds).
  Measured 2026-10-23 on the workstation (node 1, CPU 40 with its HT sibling and CPUs 41-47 idle, warmed 3 s
  before every run, arm order reversed every other round, one warm-up round + 5 rounds (4 for the UniRef50
  case, 3 for SwissProt), medians of the per-round ratios against dev cd8e0f89, the binary before the switch;
  every report identical between the arms). One thread: a 50-kb fly genomic piece x the fly proteome (`-db`,
  `-evalue 1e-5`) x0.887 (10.08 -> 8.94 s, 5 of 5 rounds below 1; cycles x0.887, instructions x0.799);
  AP027131 x AP027133 proteins (`-subject`, `-query_gencode 4`) x0.965 (3.89 -> 3.76 s); one 2.1-kb fly cDNA
  x the fly proteome x0.940; 201 reads x SwissProt (`-db`) x0.944 (48.58 -> 45.85 s); a 300-nt fly cDNA piece
  x UniRef50 2026_03 volume 00 (`-db`, `-evalue 1e-5`; the direct loop, 1.76 hits per used cell) x0.960. In
  the 50-kb case the extension with the candidate closure went from 250 to 188 cycles per candidate (46.2
  million candidates; 36.0 % of dev's cycles in `perf record`). With `=0` the instruction counts equal dev's
  (x1.0000-1.0002) and the machine code of the loops without the switch is dev's (normalized `objdump`);
  cycles x0.996-1.008 (the same code at other addresses). Max RSS within the spread of repeated runs of one
  binary (14-15 MB on the 50-kb case). On WebAssembly (the serial WASI command, wasm32-wasip1, under Node 26
  in one process: the module compiled once, a new instance per search, the default and `=0` alternating, one
  warm-up round + 8 rounds, on one booked workstation core; the direct loop, as there is no AVX2): AP027131 x
  AP027133 x0.987 (5.63 -> 5.56 s, 8 of 8 rounds below 1), 50 fly cDNAs x the yeast proteome (`-subject`)
  x0.974 (4.28 -> 4.17 s, 8 of 8), the same output bytes.
- `LOSAT_X_BXUNGAPSENTSHADOW=1` shadow of `LOSAT_X_BXUNGAPSENT` (turns the path on): every extension of the
  path also runs `extend_two_hit_blosum62` and the results (score, offsets, length, `right_extend`,
  `s_last_off`, and the early return) must be equal; any difference panics with the inputs. One line at exit,
  after a subject took the path: `[X_SHADOW] BXUNGAPSENT subjects= compared=` (`subjects=`: the pairs of a
  query batch and a subject made ready for the path).
- `LOSAT_X_BXKAONCE=1` (BLASTX, speed; in the strict set from the pull request that added it, 2026-10-26; `=0`
  turns it off): the Karlin blocks of a query batch that is not split are computed once, not twice, and the
  BLOSUM62 ideal block once per process. `search_core` takes the batch's blocks from `score_block` (in the query
  setup), runs `effective_lengths` on them and sets each context's `is_valid` to the block's `valid`;
  `split_queries` then returns, for a batch that is not split (one chunk, also without gapped), a chunk whose
  `prepared` is a clone of that batch, and the chunk loop called `score_block` on it again: the same compositions,
  so the same ungapped blocks (`valid` = `is_valid && block`, and `is_valid` is the first block's `valid`), and
  `effective_lengths` writes only `length_adjustment` and `search_space`, never `valid` or `ungapped` (it reads
  `ungapped` in ungapped mode). With the switch the serial chunk loop uses the blocks `score_block` gave, kept
  before `effective_lengths` (`search.rs`, dispatch point); the chunks of a split batch have sequences of their
  own and keep `score_block`. NCBI computes the score block once in the query setup and creates a chunk's own only
  for a split query (`prelim_stage.cpp:243-246`). The ideal block (`Blast_ScoreBlkKbpIdealCalc`: the standard
  composition against itself) takes no input, and `score_block` computed it on each call; with the switch it is a
  `OnceLock` (`parameters.rs`; a constant because BLASTX accepts only BLOSUM62), which acts on every BLASTX search.
  Reach: the own side is recorded where a batch that is not split runs through the serial chunk loop; the
  reference side (`score_block` on a chunk) only where a split query runs through that loop (1 thread, or
  `LOSAT_X_BXCHUNKFLOW`), not in the chunk waves at 2+ threads (`x_compute_chunk`, which keeps `score_block`).
  Exactness: a unit test (`kept_blocks_equal_score_block_of_the_unsplit_chunk`) compares the kept blocks with
  `score_block` of the chunk field by field, bit for bit, with invalid contexts among valid ones, for a batch of
  several queries and for one long query without gapped; `tests/blastx_unit_cap.rs` compares the reports of
  several queries with it on and off, at 1 and 4 threads.
  Measured on ws (2x Xeon 6338, node 1, one binary, default against `=0`, medians of paired rounds, byte-identical
  outputs): 800 fly cDNAs against 8 / 50 / 200 yeast proteins (`-subject`) at 1 thread x0.771 / 0.800 / 0.875
  (1.382 -> 1.066 s, 1.579 -> 1.264 s, 2.583 -> 2.261 s; instructions x0.690 / 0.723 / 0.820); against 50 and 200
  at 16 threads x0.863 / 0.860 (50 with `-seg no` x0.863); against 8 / 50 at 32 threads x0.875 / 0.918; the whole
  yeast proteome at 16 threads `-subject` / `-db` x0.993 / 0.995; a 660 kb genome (AP027131, split) against 50
  proteins x0.997 at 1 thread, x0.985 at 16. The gain is largest with many queries and few subjects; a split
  query saves only the ideal block.
- `LOSAT_X_BXPOOLN1=1` (BLASTX, memory; in the strict set from the pull request that added it, 2026-10-26; `=0`
  turns it off): a search on one thread (`-num_threads 1`) builds its word-size-3 lookup tables without the pool of
  the pooled build: `blastx::seed::Lookup::new` asks `x_bx_pooled_for(threads)` (`tblastx/lookup/x_lut_direct.rs`),
  builds with `x_build_direct` (its own neighbour memo, `vec!` vectors), as NCBI does, and
  marks the table unpooled (`XKeptPresence::pooled`), so the table's vectors are freed when the lookup drops; only
  tables built pooled go back to the pool. Searches of two threads or more build pooled as before. Why: the pool moves
  vectors between threads (no thread's glibc arena keeps a table's worth); at one thread it held a table and its
  build's temporaries through each query batch's end (the composition-based statistics and the traceback), whose own
  allocations came on top. Results do not change: both builders run `x_build_direct_in` (unit tests compare both with
  the reference builder); `=0` restores the previous choice. Reach: own = a search on one thread, ref = two threads or
  more (both required; the blastx fixtures run at 1 and 4 threads). Test: `tests/blastx_unit_cap.rs`
  (`queries_give_the_report_without_pool_n1`: 4 x 4 kb queries and one 45 kb split query, report equal to `=0` at one
  thread, reach sides at 1 and 4 threads).
  Measured on ws (2x Xeon 6338; dev f5326a59 against e85c95c3 and its `=0`; outputs byte-identical): peak `RssAnon`
  sampled every 50 ms (as v4-MEM), AP027131 (660 kb genome, 68 chunk tables) x 543 proteins (AP027133) `-subject` at
  one thread 26.0 / 22.7 / 26.0 MiB (dev / default / `=0`, 4 runs each; the MF floor is 22.5); sampled every 2 ms
  (one node-0 core): the same row 25.7-26.0 -> 24.1-24.6 (a 1 ms peak at the very end of the run, in every arm, is
  caught by some runs; 22.6-22.8 otherwise), `-db` 25.9-26.0 -> 23.9-24.8, t10 against Swiss-Prot `-db` 7.0 -> 5.7;
  the price is the `=0` form for many small batches: 800 fly cDNAs against 8 / 50 / 200 yeast proteins 12.8-13.2 ->
  13.8-14.1, against the whole proteome `-subject` / `-db` 21.1 -> 22.0 / 13.5 -> 13.9 (NCBI `RssAnon` 17 and more).
  Time at one thread (one node-0 core, TIMED under a node load of 30-95; default against `=0`, medians of paired
  rounds): wall x0.971-1.027 over eight cases (every round above 1.0 only for the proteome `-subject` x1.002 and
  `-db` x1.008, 3 rounds each; every round below for 200 proteins x0.971), cycles:u x1.001-1.006, instructions:u
  x0.999-1.000. 16 threads (node 1): AP027131 x 543
  proteins x1.004, 800 cDNAs x proteome `-subject` x1.000, peaks equal (65-66, 76-78 MiB).
  Forms not taken (line A DECISIONS 53, local): emptying the pool when the last batch leaves its chunk loop raised
  800 cDNAs x 50 proteins at one thread to 15.0-16.7 MiB (after the first free of a mapped table glibc's dynamic mmap
  threshold rises and each later table comes from the heap); freeing a table that lived alone raised AP027131
  `-subject` at 16 threads by 3-5 MiB (`LOSAT_X_BXCHUNKFLOW` builds chunk tables one at a time while other threads
  allocate).

Counters and tests: `LOSAT_X_STATS=1` (counters on stderr; DP and Newton counts need the
`xstats` cargo feature), `LOSAT_X_ADJMEMO_STATS=1` (diagnostic build; distinct composition-adjustment inputs),
`LOSAT_FUZZ_CASES=<n>` (random groups of the linking fuzz test; random problems of the X-drop DP
kernel tests).

X-drop DP capture and replay (diagnostics, not part of the strict set):

- `LOSAT_X_DPCAPTURE=<dir>` (diagnostic build) records every Nth call of each thread that reaches the vector X-drop
  DP as a self-contained problem (residues, matrix, gap costs, x_drop, result
  and edit script) in `<dir>/xdp-<pid>-t<k>.xdp`, one file per thread, appended
  (`utils/x_xdrop_capture.rs`, format in its header). `LOSAT_X_DPCAPTURE_EVERY=<n>` (diagnostic build) sets N
  (default 25). The search output does not change.
- `LOSAT_X_DPREPLAY=<dir>` input of the ignored test `replay_captured_problems`
  (`cargo test --release replay_captured_problems -- --ignored --nocapture`): every captured
  problem through the scalar loops, the 128-bit, the DPW128 and the 256-bit kernels and
  the dispatch, compared with the recorded result, then the three kernels timed (median ns per cell
  and per call, traceback and score-only, 8- and 16-bit lanes). `LOSAT_X_DPREPLAY_REPS` (default 7),
  `LOSAT_X_DPREPLAY_OUT=<file>` (append the table to a file). The ignored test
  `w128_replay_captured_problems` does the same for the 128-bit and the DPW128 kernels
  only and also runs in the wasm32-wasip1 test build (the libtest `.wasm` under a WASI host such
  as Node's `wasi` module with a root preopen), which times the simd128 kernels.
