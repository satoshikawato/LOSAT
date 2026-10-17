# LOSAT debug and diagnostics environment variables

Moved from `AGENTS.md` on 2026-10-08 so that it is read only when needed.

- `LOSAT_TRACE_HSP="qstart,qend,sstart,send"` trace a specific TBLASTX HSP.
- `LOSAT_TRACE_HSP_MASKS=1` print mask coverage for the traced TBLASTX HSP.
- `LOSAT_TRACE_CHAIN_HSP="qstart,qend,sstart,send"` trace TBLASTX chain
  selection for a specific HSP.
- `LOSAT_TRACE_LINK_SELECTIONS=1` print TBLASTX link-selection details.
- `LOSAT_LINKING_LEGACY=1` TBLASTX sum-statistics linking by the literal port
  `link_hsp_group_ncbi` instead of the default kernel (the incremental kernel of
  `tblastx/sum_stats_linking/linking_incr.rs`; same result). The linking traces
  (`LOSAT_TRACE_HSP`, `LOSAT_TRACE_CHAIN_HSP`, `LOSAT_TRACE_LINK_SELECTIONS`,
  `LOSAT_DEBUG_CHAINING`) also run the literal port, which prints them.
- Checks of the default TBLASTX linking kernel: `LOSAT_LINK_FAST_VERIFY=1` (every
  pass against a pass computed from scratch by NCBI's scans; slow on large
  groups), `LOSAT_LINK_FAST_SHADOW=1` (every group also linked by the literal
  port, all fields compared), `LOSAT_LINK_WIDE_KEYS=1` (u128 keys for every
  group, to measure their cost), `LOSAT_LINK_STATS=1` (per-group counters,
  `[LINK_INCR_STATS]` lines).
- `LOSAT_DUMP_TBLASTX_STAGE=<dir>` append TBLASTX stage snapshots as TSV files.
- `LOSAT_TRACE_BLASTN_HSP="qstart,qend,sstart,send"` trace a specific BLASTN
  HSP.
- `LOSAT_TRACE_BLASTN_SEED="q,s"` trace a specific BLASTN seed.
- `LOSAT_TRACE_BLASTN_CONTEXT=<context_idx>` restrict BLASTN tracing by context.
- `LOSAT_TRACE_BLASTN_SUBJECT=<subject_id_or_index>` restrict BLASTN tracing by
  subject.
- `LOSAT_TRACE_BLASTN_STAGE=<seed|ungapped|prelim|traceback|purge|hitlist|all>`
  restrict BLASTN tracing by stage.
- `LOSAT_DEBUG_CUTOFFS=1` cutoff calculations (tblastx + blastn).
- `LOSAT_DEBUG_CUTOFFS_ALL=1` verbose TBLASTX cutoff diagnostics.
- `LOSAT_DEBUG_CHAINING=1` chaining debug (legacy; tblastx).
- `LOSAT_DEBUG_EXTENSION=1` tblastx extension debug.
- `LOSAT_DEBUG_HSP_SAVING=1` TBLASTX HSP-save diagnostics.
- `LOSAT_DEBUG_OUTPUT_FILTER=1` TBLASTX output filter diagnostics.
- `LOSAT_DEBUG_BLASTN=1` blastn hit loss diagnostics.
- `LOSAT_DEBUG_COORDS=1` blastn coordinate transforms.
- `LOSAT_DEBUG_COORDS_START=<int>` narrow selected BLASTN coordinate diagnostics.
- `LOSAT_DEBUG_SCAN_SOFF=<int>` tblastx scan debug center subject offset.
- `LOSAT_DEBUG_SCAN_WINDOW=<int>` tblastx scan debug window size.
- `LOSAT_TIMING=1` timing breakdown. TBLASTX prints, per query batch, after its other lines:
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
  cannot go below; `largest_group` is the number of HSPs in the largest group. A linking call
  that takes the reference body (`LOSAT_X_HSPINPLACE=0`, `LOSAT_X_DEFAULT=0`,
  `LOSAT_DIAGNOSTICS`, `LOSAT_DEBUG_CUTOFFS`) is not timed there: `calls=0`;
  `tblastx/x_link_timing.rs`). The clock is not read without the variable.
  BLASTX prints four kinds of lines at the end of a run (`blastx/x_timing.rs`):
  `[TIMING] blastx_run: startup .. prolog .. read_stores .. loop .. epilog .. total .. (search .. render .. titles .. input ..)`
  (wall times of the thread that runs the program; inside `loop`, the time it waits for searches, writes reports (of
  which `titles` reads the reported sequences' titles with `-db`), and the rest, `input`);
  `[TIMING] blastx_waves: waves=N batches=N wall .. busy .. threads=N idle .. tail=N tail_wall ..`
  (with `LOSAT_X_BXBATCH`: the waves of query batches searched side by side, the searches' own times added over the
  threads, `idle = wall x threads - busy`, and the batches searched one at a time with the whole pool: the end of the
  input, or a long query);
  `[TIMING] blastx_batch_side:` / `blastx_batch_alone: batches=N query_setup .. setup .. chunks_par .. lookup .. scan_par .. scan_inline .. replay .. chunk_end .. batch_end .. sum ..`
  (the stages of a query batch in the order they run; `side` = batches searched on one thread beside others, times
  added over the threads; `alone` = a batch searched by itself, with the pool if there is one, wall times; `replay`
  hands the subjects to the sink in subject order and includes `scan_inline`) and
  `[TIMING] blastx_batch_end_side:` / `blastx_batch_end_alone: init .. read .. kappa_spec .. kappa_serial .. loop .. finish .. (matches=N reused=N serial=N)`
  (inside the end of a batch: the composition-based redo made ahead on the pool, the redo done in the in-order loop,
  and the hit lists after it).
- `LOSAT_DIAGNOSTICS=1` general diagnostics counters.
- `LOSAT_STARTUP_TRACE=1` startup trace.
- `LOSAT_WASI_THREADS_DEBUG=1` threaded-WASI scheduling diagnostics.
- `LOSAT_WASI_PRESPAWN=auto|N` (Node threaded runner, opt-in) boot N thread
  workers (`auto`: `-num_threads` - 1) before `_start`; same output, extra
  `standby_*` debug events. See `LOSAT/tests/README.md`, "Threading contract gates".
- `LOSAT_WASI_GUARD_CACHE=<dir>` (Node threaded runner, opt-in) reuse the
  shared-memory guard's output from `<dir>`, checked by SHA-256 and rebuilt
  when it does not match.
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
strict set are on by default (43 then, 72 since `LOSAT_X_PAIRLUT`, `LOSAT_X_HSPINPLACE`,
`LOSAT_X_SEEDBUCKET`, `LOSAT_X_TBNOIDCHUNK`, `LOSAT_X_SEEDRELEASE`, `LOSAT_X_GREEDYFUSE`,
`LOSAT_X_LINKSORT`, `LOSAT_X_NEWTONLANES8`, `LOSAT_X_TBNSTATE`, `LOSAT_X_TBNXLATE`,
`LOSAT_X_TBNOIDNEXT`, `LOSAT_X_TBNOIDUNITS`, `LOSAT_X_LINKSLICES`, `LOSAT_X_BPSUBJRUN`,
`LOSAT_X_TBNSTORE`, `LOSAT_X_DPNEON`, `LOSAT_X_TBNBATCHRULE` and the eight blastp switches
`LOSAT_X_BPSEED`, `LOSAT_X_BPUNGAP`, `LOSAT_X_BPGAP`, `LOSAT_X_BPLUT`, `LOSAT_X_BPPFREE`,
`LOSAT_X_BPPCOMPACT`, `LOSAT_X_BPREDOFIT` and `LOSAT_X_BPQBLOCK`, `LOSAT_X_BNSTREAM`,
`LOSAT_X_TBXNEXT` and `LOSAT_X_TBXFRAMES`, and `LOSAT_X_BPDIAG` joined), and their variables mean
(`LOSAT/src/utils/xswitch.rs`):

| `LOSAT_X_<NAME>` | `LOSAT_X_DEFAULT` | the switch is |
| --- | --- | --- |
| `0` | any | off |
| any other non-empty value | any | on (a value switch takes the number if it is in range, otherwise its built-in value) |
| not set, or empty | `0` | off, the default before the change |
| not set, or empty | not set, empty, or anything but `0` | on |

- `LOSAT_X_DEFAULT=0` runs the previous (reference) paths everywhere. With some variables set
  as well, it gives the run those variables gave before the default changed.
- `LOSAT_X_<NAME>=0` turns one switch off, to find which path a problem comes from. Single
  switches off are a diagnostic aid: `DPAVX2`, `DPW128` and `DPNEON` act only with `DPFAST`; `TBNSSIDE`,
  `TBNQSIDE` and `TBNBATCH` only with `TBNPAR` (`TBNEVENTS` only without it); `TBNBATCH` needs
  `TBNSSIDE`; `BXBATCH` needs `BXPOOL`; `AHEAD` needs `GREEDYSPEC` for megablast; `MBBATCH` runs
  where `MBDELAY` does not; `PAIRLUT` acts only while `PAIRPAR` is on; `SEEDBUCKET` acts only
  in the scan loop, which TBLASTX runs where `PAIRPAR` does not take the search;
  `NEWTONLANES8` acts only with `NEWTONLANES`, and only on aarch64; `TBNSTATE` and `TBNXLATE`
  act only with `TBNPAR`; `TBNBATCHRULE` only with `TBNBATCH`; `TBNOIDNEXT` and `TBNOIDUNITS` act
  only with `TBNOIDCHUNK`; `BPSUBJRUN` acts only in the parallel `-subject` loop of blastp (2 or
  more threads, several subjects); `BPSEED` acts on the two-hit seed stage only; `BPDIAG` acts only
  in `BPSEED`'s fused loop and batches only with AVX2 (else the per-hit loop); `BPGAP` runs the
  vector alignment only with AVX2, BMI1 and LZCNT (else the reference function); `BPLUT` declines in
  diagnostics mode; `BPQBLOCK` takes the block path only for one `-outfmt 6` report with more than
  one block (see its entry); `BNSTREAM` only in blastn's `-subject` search at 2 or more threads;
  `TBXNEXT` only in TBLASTX's parallel subject loop (2 or more threads, several subjects, native
  builds).
- White space around a value is ignored (`set LOSAT_X_DEFAULT=0 && ...` in cmd.exe stores `0 `).
  Only `0` means off; `00`, `false` or `off` are "any other value".
- The variables are read once per process (once per Wasm instance), except `GREEDYSPEC`,
  `MBDELAY`, `MBBATCH` and `TBNBATCH`, which are read once per search or lookup build. The Web
  reactors get an empty WASI environment, so they always run the default.
- The `...SHADOW` / `...VERIFY` switches run the new and the reference computation side by side
  and stop with an assertion on the first difference. They, the counters, the timing aids and the
  tuning variables are not part of the set: they are off unless set to any value, and a shadow
  switch turns its switch on whatever the two variables above say (except
  `LOSAT_X_AHEADSHADOW`, which compares only while `AHEAD` is on, and
  `LOSAT_X_TBNXLATESHADOW`, which compares only where `TBNPAR` makes the subject plans). A
  shadow compares the kernel the default selects (for `LOSAT_X_DPSHADOW` the 256-bit, `DPW128`
  or base kernel that takes the call).
- The names are provisional (`x_` prefixes, `EXPERIMENT` comments). The 72 single variables are
  to be removed once the default has run on every CI environment and through a few releases;
  `LOSAT_X_DEFAULT` stays. In the descriptions below `LOSAT_X_<NAME>=1` shows the switch, not a
  setting that is needed, and "strict-set candidate" marks the round-5 switches, which are in
  the set.

The strict set (the variables that are on; `STRICT_SET` in `xswitch.rs`), written as the
environment that turned it on before the change. The three numbers are the built-in values:

```bash
export LOSAT_X_DPFAST=1 LOSAT_X_COMPFAST=1 LOSAT_X_PURGEFAST=1 LOSAT_X_ITREEFAST=1 LOSAT_X_CTXFAST=1
export LOSAT_X_GREEDYFAST=1 LOSAT_X_GREEDYSPEC=1 LOSAT_X_AHEAD=16 LOSAT_X_ENVCACHE=1
export LOSAT_X_DUSTFAST=1 LOSAT_X_DUSTRING=1 LOSAT_X_MBDELAY=1 LOSAT_X_MBBATCH=64
export LOSAT_X_SEGFAST=1 LOSAT_X_SEGMEMO=1 LOSAT_X_SEGSHARE=1 LOSAT_X_KARLINFAST=1 LOSAT_X_IDEALMEMO=1
export LOSAT_X_ERFMEMO=1 LOSAT_X_DEKKER=1 LOSAT_X_LUTARENA=1 LOSAT_X_LUTSPLIT=1
export LOSAT_X_TBNPAR=1 LOSAT_X_TBNEVENTS=1 LOSAT_X_TBNQSIDE=1 LOSAT_X_TBNSSIDE=1 LOSAT_X_TBNBATCH=1
export LOSAT_X_BXPAR=1 LOSAT_X_BXLAZYCTX=1 LOSAT_X_BXLEAN=1 LOSAT_X_BXPOOL=1 LOSAT_X_BXCHUNK=64 LOSAT_X_BXBATCH=1
export LOSAT_X_CODONFAST=1
export LOSAT_X_NEWTONEXACT=1
export LOSAT_X_DPAVX2=1 LOSAT_X_DPW128=1
export LOSAT_X_BXLUT=1 LOSAT_X_BXSCAN=1 LOSAT_X_BXSEGMEMO=1
export LOSAT_X_PAIRPAR=1 LOSAT_X_NEWTONLANES=1 LOSAT_X_GREEDYPRELIM=1
export LOSAT_X_PAIRLUT=1 LOSAT_X_HSPINPLACE=1
export LOSAT_X_SEEDBUCKET=1
export LOSAT_X_TBNOIDCHUNK=1
export LOSAT_X_SEEDRELEASE=1
export LOSAT_X_GREEDYFUSE=1
export LOSAT_X_LINKSORT=1
export LOSAT_X_NEWTONLANES8=1
export LOSAT_X_TBNSTATE=1 LOSAT_X_TBNXLATE=1
export LOSAT_X_TBNOIDNEXT=1 LOSAT_X_TBNOIDUNITS=1
export LOSAT_X_LINKSLICES=1
export LOSAT_X_BPSUBJRUN=1
export LOSAT_X_TBNSTORE=1
export LOSAT_X_DPNEON=1
export LOSAT_X_TBNBATCHRULE=1
export LOSAT_X_BPSEED=1 LOSAT_X_BPUNGAP=1 LOSAT_X_BPGAP=1 LOSAT_X_BPLUT=1
export LOSAT_X_BPPFREE=1 LOSAT_X_BPPCOMPACT=1 LOSAT_X_BPREDOFIT=1 LOSAT_X_BPQBLOCK=50000
export LOSAT_X_BNSTREAM=1
export LOSAT_X_TBXNEXT=1 LOSAT_X_TBXFRAMES=1
export LOSAT_X_BPDIAG=1
```

Value switches: `LOSAT_X_AHEAD` (2 to 65,536, built-in 16), `LOSAT_X_MBBATCH` (2 to 4,096,
built-in 64), `LOSAT_X_BXCHUNK` (1 or more, built-in 64), `LOSAT_X_BPQBLOCK` (1 or more, built-in
50,000: query letters per block; a value that is not a positive number is the built-in value).

`LOSAT_X_SEEDBUCKET` was in the block of the round 3-4 set and was left out when the set became
the default: it numbered the lookup hits of one chunk scan in a 32-bit counter, exact only up to
2^32 hits per scan, which a query and a subject of a few megabases pass. It now numbers the hits
per flush, which is exact for any number of hits, and joined the set on 2026-10-12. With the set
on it acts on the TBLASTX searches that `LOSAT_X_PAIRPAR` leaves to the scan loop and whose
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
the unit path of `LOSAT_X_PAIRPAR` then builds only the tables it scans with, and the HSP list
of a subject is sorted and linked without a second whole copy. `LOSAT_X_PAIRLUT` acts only
while `LOSAT_X_PAIRPAR` is on.

`LOSAT_X_TBNOIDCHUNK` was added after the default changed, as a switch that was off unless set,
and joined the set on 2026-10-13. It changes only which subjects one entry into the tblastn
search pool takes (2 or more threads, several subjects) and needs no other switch. Since
2026-10-14 its chunks are bounded by residues, as NCBI's database reader cuts them for its
search threads.

`LOSAT_X_SEEDRELEASE` was added after the default changed, as a switch that was off unless set,
and joined the set on 2026-10-13. It lowers the peak memory of TBLASTX further: the tables that
only the seed stage reads are freed before the HSP list of a query batch's last subject is
merged and linked. It changes no value and no order of computation, acts on every TBLASTX path
(the unit path of `LOSAT_X_PAIRPAR` and the scan loop) and needs no other switch.

`LOSAT_X_GREEDYFUSE` was added after the default changed, as a switch that was off unless set,
and joined the set on 2026-10-14. It runs the preliminary score-only greedy alignment of megablast
as one fused loop over distances. It changes no value and no order of computation, acts on
the blastn searches whose preliminary stage takes the score-only greedy alignment on an unpacked
subject (`LOSAT_X_GREEDYPRELIM`, which is in the set), and its shadow
`LOSAT_X_GREEDYFUSESHADOW` stays outside the set.

`LOSAT_X_LINKSORT` was added after the default changed, as a switch that was off unless set, and
joined the set on 2026-10-14. It speeds up TBLASTX at two or more threads: the linking call sorts
its output order once where it sorted twice, and a long index sort runs as a stable parallel sort
on the search pool. It changes only how an index order is computed, not the order, and acts only
where `LOSAT_X_HSPINPLACE` does (the reference bodies have their own sorts). The shadow
`LOSAT_X_LINKSORTSHADOW` stays opt-in.

`LOSAT_X_NEWTONLANES8` was added after the default changed, as a switch that was off unless set,
and joined the set on 2026-10-15. On aarch64 the Newton lanes of `LOSAT_X_NEWTONLANES` (blastp and
tblastn) solve eight problems per round instead of four. It changes only which independent
problems share a vector, not the operations of a problem, and acts only where
`LOSAT_X_NEWTONLANES` runs. On x86-64 and Wasm it does nothing: the lanes stay four wide there.
`LOSAT_X_DEFAULT=0` runs no lanes (the reference Newton port), so on the arm64 runners
platform-checks compares eight lanes with no lanes.

`LOSAT_X_TBNSTATE` and `LOSAT_X_TBNXLATE` were added after the default changed, as switches that
were off unless set, and joined the set together on 2026-10-15. They cut the work tblastn does per
subject on the unit path (`LOSAT_X_TBNPAR`): a thread keeps the scan state of the units from one
subject to the next instead of making and clearing it per subject and per unit, and the six frames
of a subject are read from codon tables instead of through letters. Neither changes a value; they
are timed together because `TBNSTATE` alone was 3 % slower on the fly genome at 8 threads. The
shadow `LOSAT_X_TBNXLATESHADOW` stays opt-in; `LOSAT_X_TBNPARSHADOW` compares the unit path, with
both switches, with the serial reference.

`LOSAT_X_TBNOIDNEXT` and `LOSAT_X_TBNOIDUNITS` were added after the default changed, as switches
that were off unless set, and joined the set together on 2026-10-15. They change how tblastn hands
the chunks of `LOSAT_X_TBNOIDCHUNK` to its threads (2 or more threads, several subjects): the
threads take the next chunk when they are free instead of waiting at the end of each entry, and
the units of a subject longer than a chunk's budget (a chromosome) run on the pool instead of on
the thread that took it. Neither changes a value or the order in which the results are entered;
they are timed together because `TBNOIDUNITS` alone was slower on the human genome at 8 and 16
threads. On main's path before `TBNSTATE` joined the set they cost 5-11 % more CPU time on genomes;
with it the CPU time is the same.

`LOSAT_X_LINKSLICES` was added after the default changed, as a switch that was off unless set, and
joined the set on 2026-10-15. It lowers the memory of TBLASTX and, at two or more threads, its time:
the linking call links the groups of a subject's HSP list where they lie in the one sorted list,
where it cut the list into one Vec per group and copied the linked groups into a new list, so the
list no longer exists twice. It changes where the groups are linked, not a value or an order, and
acts only where `LOSAT_X_HSPINPLACE` does. The shadow `LOSAT_X_LINKSLICESSHADOW` stays opt-in.

`LOSAT_X_BPSUBJRUN` was added after the default changed, as a switch that was off unless set, and
joined the set on 2026-10-16. With two or more threads blastp's `-subject` loop searches the
subjects in runs of consecutive subjects, each run with one diagonal table that is cleared once and
one search scratch per thread, instead of a new scratch for every piece rayon splits the subjects
into and a clearing of the table for every subject. It changes no value; a run's later subjects
start at the diagonal offset the one-thread loop has. On the workstation (16 physical cores, medians
of alternating pairs, every output hash equal, max RSS equal) `blastp-fast` of the human proteome
against the yeast proteome took 0.89 of the time at 16 threads (6 pairs; 0.88 with
`LOSAT_X_BPQBLOCK=20000`), 0.90 and 0.92 at 8 and 0.92 at 4, exact `blastp` 0.998 at 16 (8 pairs)
and 0.995 at 8 (4 pairs). One thread and `-db` are not affected. The shadow
`LOSAT_X_BPSUBJRUNSHADOW` stays opt-in.

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
joined the set on 2026-10-15. On aarch64 it lets `LOSAT_X_DPFAST` (and `LOSAT_X_DPW128`) run the
vector X-drop kernels built over NEON, the kernel bodies x86_64 and wasm32 run. It changes no
value and no order of computation, acts on aarch64 only (it does nothing on other targets) and
only where `LOSAT_X_DPFAST` is on; `LOSAT_X_DPSHADOW` compares the NEON kernels whatever it says.

`LOSAT_X_TBNBATCHRULE` was added after the default changed, as a switch that was off unless set,
and joined the set on 2026-10-16. With `LOSAT_X_TBNBATCH` it chooses which tblastn query batches
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
50,000 query letters in more than one query batch now runs in query blocks. Its shadow
variables (`LOSAT_X_BPSEEDSHADOW`, `LOSAT_X_BPUNGAPSHADOW`, `LOSAT_X_BPGAPSHADOW`,
`LOSAT_X_BPLUTSHADOW`) and `LOSAT_X_BPPCOMPACT=shadow` stay outside the set, and so do
`LOSAT_X_BPCLUT` and `LOSAT_X_BPSEED_PREFETCH`.

`LOSAT_X_BNSTREAM` was added after the default changed, as a switch that was off unless set, and
joined the set on 2026-10-16. With two or more threads blastn hands each subject's HSPs of a
subject set (`-subject`) to the hit lists in OID order from a window of at most 64 subjects per
thread while the threads search the next ones, instead of holding the HSPs of every subject of a
query batch until all are searched (one thread always hands them over at once). It changes no
value and not the order in which the hit lists take the subjects; it lowers the peak memory and,
on the workstation, the time (0.94 at 8 threads and 0.88 at 16 on 100 16S queries against 28,036
16S sequences). A database search (`-db`) has its own drivers and is not affected.

`LOSAT_X_TBXNEXT` and `LOSAT_X_TBXFRAMES` were added after the default changed, as switches that
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

What runs per target:

- x86_64: the vector X-drop DP (`DPFAST`: SSE4.1 or AVX 128-bit; `DPAVX2`: 256-bit with AVX2;
  `DPW128`), the AVX2 builds of the greedy rows and of the Newton kernels where the CPU has AVX2,
  `MBDELAY` with the prefetch instruction.
- wasm32 (built with simd128): the simd128 X-drop DP (`DPW128` replaces the first kernel),
  `MBBATCH` instead of `MBDELAY`, `DEKKER`; everything else is the same Rust as on x86_64
  without AVX2.
- aarch64: the NEON X-drop DP (`DPFAST` with `DPNEON`, and `DPW128` where the 128-bit kernel
  would run; `DPAVX2` does nothing), `MBBATCH` instead of `MBDELAY`; `DEKKER` does nothing (it
  is compiled for wasm32 only); `NEWTONLANES8` gives the Newton lanes eight problems per round;
  everything else as on wasm32.
- Any other target: no vector X-drop DP is compiled, so `DPFAST`, `DPAVX2`, `DPW128` and
  `DPNEON` do nothing (the scalar port runs), and the Newton lanes stay four wide; otherwise as
  on aarch64.
- One thread (`-num_threads 1`, a build without the `parallel` feature, wasm32 without
  `wasm-threads`): the scheduling switches (`AHEAD`, `GREEDYSPEC`, `PAIRPAR` for
  blastn/megablast, `SEGSHARE`, `TBNBATCH`, `TBNOIDCHUNK`, `TBNOIDNEXT`, `TBNOIDUNITS`, `BXPAR`,
  `BXPOOL`, `BXBATCH`, `BXCHUNK`, `BNSTREAM`) leave the reference code running. The unit paths have no thread-count gate: the tblastn `TBNPAR` units
  (with `TBNSSIDE` / `TBNQSIDE`; only its parallel redo needs a pool) and the TBLASTX `PAIRPAR`
  units run in order on the calling thread.

Same values, faster computation:

- `LOSAT_X_DPFAST=1` X-drop gapped DP (score-only and traceback) in 128-bit vectors (16 x i8 or
  8 x i16 lanes; AVX2/SSE on x86_64, simd128 on Wasm, NEON on aarch64 with `LOSAT_X_DPNEON`,
  which is in the set),
  blastn task blastn and the protein programs; calls that do not fit 8/16 bits use the original
  loop (`utils/xdrop_simd.rs`).
  Shadow: `LOSAT_X_DPSHADOW=1`. Timing aids: `LOSAT_X_DPNO8=1` (no 8-bit kernel),
  `LOSAT_X_DPNOAVX=1` (no AVX2).
- `LOSAT_X_DPAVX2=1` (with `LOSAT_X_DPFAST`; strict-set candidate) the same X-drop DP in 256-bit
  AVX2 vectors (32 x i8 or 16 x i16 lanes, x86_64 with AVX2 only): the same values, lane width
  and band steps as the 128-bit kernels, so the same result and edit script; a call the 256-bit
  kernel does not take runs on the 128-bit kernel (`utils/x_xdrop_avx2.rs`). Covered by
  `LOSAT_X_DPSHADOW=1`; `LOSAT_X_DPNOAVX=1` turns it off as well.
- `LOSAT_X_DPW128=1` (with `LOSAT_X_DPFAST`; strict-set candidate) the 128-bit X-drop DP kernels
  with the per-row savings of `LOSAT_X_DPAVX2` (first and last survivor kept while the row runs,
  row gap as a prefix maximum with a two-operation carry, the first lane's missing diagonal as a
  lane mask, column residues copied 16 at a time, row tables of `Rows28` / `Flat` matrices built
  with vector instructions, byte shifts in the form V8 compiles to `palignr`): the same values,
  lane width, blocks and result as the 128-bit kernels; a call it does not take (8-bit lanes with
  `16 * gap_extend - (gap_open + gap_extend) > 128`) runs on the 128-bit kernel
  (`utils/x_xdrop_w128.rs`). On x86_64 it runs where the 128-bit kernel would (after
  `LOSAT_X_DPAVX2` when both are set on an AVX2 machine, which keeps priority), on the SSE4.1 or
  AVX-128 build; on wasm32 it replaces the simd128 128-bit kernel (the only DP kernel there).
  Covered by `LOSAT_X_DPSHADOW=1`.
- `LOSAT_X_DPNEON=1` (aarch64 only, speed; in the strict set since 2026-10-15, `=0` turns it off)
  lets `LOSAT_X_DPFAST` run on aarch64: the 128-bit kernels and the `LOSAT_X_DPW128` kernels
  built over a NEON `arch` / `ext` layer (`mod arch` in `utils/xdrop_simd.rs`, `mod ext` in
  `utils/x_xdrop_w128.rs`; the kernel bodies are the ones x86_64 and wasm32 run). Off, aarch64
  runs the scalar port, as before the NEON layer. `LOSAT_X_DPSHADOW=1` runs and compares the NEON
  kernels on aarch64 whatever it says. Platform checks hash the default (NEON) and, on the Linux
  and macOS arm64 runners, `LOSAT_X_DPNEON=0` (mode `noneon`: the other switches on, the scalar
  X-drop).
- `LOSAT_X_DPNEON2=1` (aarch64 only, speed; not in the strict set: off unless set, to any value, and
  `LOSAT_X_DEFAULT` does not touch it) the X-drop DP kernels of `LOSAT_X_DPAVX2` on NEON, each
  256-bit vector in two registers (32 x i8 or 16 x i16 lanes per block, the same band steps, the
  row gap in the form of `LOSAT_X_DPW128`): the same values and result as the 128-bit kernels; a
  call they do not take (8-bit lanes with `32 * gap_extend - (gap_open + gap_extend) > 128`, 16-bit
  lanes with `16 * gap_extend - (gap_open + gap_extend) > 32768`) runs on the `LOSAT_X_DPW128` /
  128-bit kernels (`utils/x_xdrop_neon2.rs`). It turns the NEON kernels on as `LOSAT_X_DPNEON`
  does and takes the place `LOSAT_X_DPAVX2` has on x86_64 (first in the dispatch). Covered by
  `LOSAT_X_DPSHADOW=1` when set. Platform checks hash it on the Linux and macOS arm64 runners
  (mode `dpneon2`).
- `LOSAT_X_GREEDYFAST=1` megablast greedy alignment, one distance in three passes of 256 cells,
  match extension 8 bases at a time. Shadow: `LOSAT_X_GREEDYSHADOW=1`.
- `LOSAT_X_GREEDYPRELIM=1` (strict-set candidate, S-G, round 5) the preliminary (score-only)
  megablast greedy alignment reads a one-byte-per-base copy of the packed ncbi2na subject (the
  NCBI-resolved codes) with `compressed_subject = false` instead of the packed bytes: the same
  comparisons, bounds and results (`algorithm/blastn/alignment/x_greedy_prelim.rs`). Only with
  `LOSAT_X_GREEDYFAST` does the preliminary stage then use the fast rows; without it the ported
  loop runs on the copy. Memory: one byte per subject base while the subject is searched (made on
  the first chunk with ungapped hits, shared by its chunks and threads). Shadow:
  `LOSAT_X_GREEDYPRELIMSHADOW=1` (implies the switch) also runs the packed call and stops on the
  first difference; at exit `[X_SHADOW] GREEDYPRELIM calls_compared=N`. With `LOSAT_X_STATS=1`:
  `[X_STATS] GREEDYPRELIM calls= subjects_unpacked= bytes_unpacked=`.
- `LOSAT_X_GREEDYFUSE=1` (in the strict set since 2026-10-14; `=0` turns it off) the score-only megablast greedy
  alignment on an unpacked subject (the preliminary stage with `LOSAT_X_GREEDYPRELIM`): the whole
  loop over distances of `BLAST_GreedyAlign` in one function, one pass over the diagonals of a
  distance, eight bases per slide step, the same cells of the two rolling rows and of `max_score`
  written in the C order (`algorithm/blastn/alignment/x_greedy_fused.rs`, exactness argument in
  its module doc). It replaces the per-distance `LOSAT_X_GREEDYFAST` row there; traceback, packed
  subjects and the affine alignment keep their paths. x86_64 with AVX2/BMI runs a build with
  `tzcnt`/`lzcnt`, other CPUs and targets the plain build of the same code. Shadow:
  `LOSAT_X_GREEDYFUSESHADOW=1` (implies the switch) runs it on copies, runs the reference call as
  well and stops on the first difference in the return value, align lengths, seed, fence flag,
  rows or `max_score`; at exit `[X_SHADOW] GREEDYFUSE compared=N`.
- `LOSAT_X_COMPFAST=1` composition adjustment (Newton) loops reordered with the same element
  operations in the same order. Shadow: `LOSAT_X_COMPSHADOW=1`.
- `LOSAT_X_NEWTONEXACT=1` the Newton iteration of the composition adjustment with fixed-size
  storage, a structured right-looking column Cholesky and AVX2, same element operations
  (`core/composition_adjustment/x_newton_exact.rs`). Shadow: `LOSAT_X_NEWTONEXACTSHADOW=1`
  (convergence state and the 400 values compared bitwise). `LOSAT_X_NEWTONEXACT_LEVEL=0|2`
  (scalar | try AVX-512), `LOSAT_X_NEWTONEXACT_PROF=1` (cycle counters per phase).
- `LOSAT_X_NEWTONLANES=1` (blastp and tblastn) Newton problems of the composition adjustment
  solved ahead, four problems side by side (one per AVX2 lane, each lane with the operations of
  `x_newton_exact`; `core/composition_adjustment/x_newton_lanes.rs`): blastp the first adjustment of
  the next 16 matches of a query (serial and query-parallel redo; the single-query match-parallel
  redo is unchanged), tblastn the adjustments of the coming windows of a match (at least 32 problems;
  the window ranges, translation + SEG, computed ahead are handed to the redo loop). A call uses a
  stored result only when its input (q, row and column probabilities, relative entropy) has the same
  bits; otherwise the existing path runs. Shadow: `LOSAT_X_NEWTONLANESSHADOW=1` (every stored result
  that is used is compared with the reference port, status and the 400 values bitwise; every kept
  tblastn range with a fresh `get_range`). With `LOSAT_X_STATS=1` it prints
  `[X_STATS] NEWTONLANES prefetched= hits= misses= unused= ...` (no `xstats` feature needed).
  Strict-set candidate (S-F, round 5).
- `LOSAT_X_NEWTONLANES8=1` (aarch64 only; in the strict set since 2026-10-15) the Newton lanes of
  `LOSAT_X_NEWTONLANES` solve eight problems per round instead of four (`F8`, four 128-bit NEON
  operations per value) in batches of five or more problems; once the queue is empty and at most
  four lanes hold a problem, those problems move to four lanes (bit copies of the values kept
  between rounds) and finish there. Each problem keeps its operation sequence, so the results are
  the same bits (`x_newton_lanes.rs`, module header). It only widens the lanes: with
  `LOSAT_X_NEWTONLANES` off it does nothing, and `LOSAT_X_NEWTONLANES8=0` runs four lanes. On other
  targets it is ignored: on x86-64 the two 256-bit operations per value gained nothing (0 to 3 %
  slower on blastp / tblastn), and Wasm is not measured. Measured 2026-10-15 on the GitHub runners,
  bench5's blastp_wssv / tblastn_small at one thread, wall relative to four lanes: Linux arm64
  (Neoverse-N2) 0.981 / 0.959, macOS arm64 0.974 / 0.893 (macOS rounds are noisy); on Neoverse-N2
  `LOSAT_X_NEWTONLANES=0` gave 0.967 / 0.941. Shadow: `LOSAT_X_NEWTONLANESSHADOW=1` compares the
  eight-lane results too. With `LOSAT_X_STATS=1` it adds
  `[X_STATS] NEWTONLANES8 batches8= rounds8= tail_rounds=` (rounds on eight lanes and after the move).
- `LOSAT_X_SEGFAST=1`, `LOSAT_X_SEGMEMO=1` SEG window as a residue histogram updated in O(1),
  entropy memoised per state vector. Shadow (`SEGFAST` path): `LOSAT_X_SEGSHADOW=1`.
- `LOSAT_X_DUSTFAST=1`, `LOSAT_X_DUSTRING=1` DUST perfect-interval list built in one merge,
  window kept in a ring buffer. Shadow: `LOSAT_X_DUSTSHADOW=1`.
- `LOSAT_X_KARLINFAST=1` Karlin K convolution in eight lanes (same per-lane order).
- `LOSAT_X_ERFMEMO=1` Spouge ErfC memoised per thread by the argument's bits.
- `LOSAT_X_DEKKER=1` (Wasm only) ErfC product error by Dekker's product instead of `fma`.
- `LOSAT_X_CODONFAST=1` an unambiguous codon translated by one table read.
- `LOSAT_X_PURGEFAST=1` blastn common-endpoint purge without shifting the array per removal.
- `LOSAT_X_ITREEFAST=1` blastn interval tree with side indexes (grid, endpoint bitmaps) that
  answer "none contained" / "no such endpoint" without walking the tree. Shadow:
  `LOSAT_X_ITREESHADOW=1`.
- `LOSAT_X_CTXFAST=1` blastn query-offset-to-context lookup by binary search instead of a per-base
  table. Shadow: `LOSAT_X_CTXSHADOW=1`.
- `LOSAT_X_LUTARENA=1` protein lookup-table chains carved from one arena per table.
- `LOSAT_X_MBDELAY=1`, `LOSAT_X_MBBATCH=64` megablast lookup construction touches (prefetches)
  hash cells ahead in batches of 64; same update order.
- `LOSAT_X_ENVCACHE=1` blastn debug variables read once instead of per extension.
- `LOSAT_X_BXLAZYCTX=1` blastx seed context lookup moved after the diagonal test.
- `LOSAT_X_BXLEAN=1`, `LOSAT_X_TBNEVENTS=1` do not build records nobody reads (candidate copies
  without an observer, test events); no tree or DP space for a (chunk, subject) without HSPs.
- `LOSAT_X_BXLUT=1` blastx lookup table filled in place from per-cell counts (same append order as
  `s_AddNeighboringWords` and `BlastAaLookupFinalize`), neighbour-word lists memoised per thread
  (`tblastx/lookup/x_lut_direct.rs`). Shadow: `LOSAT_X_BXLUTSHADOW=1` (every table compared field by
  field). Strict-set candidate (S-C, round 5).
- `LOSAT_X_BXSCAN=1` blastx seed loop reads the scanner's per-thread offset-pair buffer in place instead
  of a fresh buffer and a copy per subject (`blastx/x_seed_scan.rs`); not in the diagnostic stages, and
  only with the reference diagonal order (no `LOSAT_X_BXSEEDBUCKET`, no `LOSAT_X_SEEDBUCKETSHADOW`). Shadow: `LOSAT_X_BXSCANSHADOW=1`
  (hit lists and diagonal table). Strict-set candidate (S-C, round 5).
- `LOSAT_X_SEEDBUCKET=1` (in the strict set since 2026-10-12) TBLASTX two-hit stage: hits buffered
  per diagonal range and processed
  range by range, HSPs restored to scan order (`tblastx/x_seed_bucket.rs`); only for diagonal tables
  of at least `LOSAT_X_SEEDBUCKET_MIN_CELLS` cells (default 2^21). Shadow:
  `LOSAT_X_SEEDBUCKETSHADOW=1` (HSP lists and the whole diagonal table; its summary line gives
  `max_scan_hits`, the most lookup hits one scan went through). Tuning:
  `LOSAT_X_SEEDBUCKET_CELLS` (cells per bucket), `LOSAT_X_SEEDBUCKET_BUDGET` (hits per flush;
  values above 2^31 are read as 2^31, and memory, 12 bytes per hit reserved up front, limits a
  real setting long before that). A shadow run compares something only where the scan loop runs:
  with the default set use `LOSAT_X_PAIRPAR=0` and check that `chunks_compared` is not 0. With
  `LOSAT_X_STATS=1` the line `[X_STATS] SEEDBUCKET bucketed_scans=... flushes=...` tells how many
  scans ran in the bucketed order (`LOSAT/tests/tblastx_seed_bucket.rs` uses it).
  BLASTX runs the bucketed order only with `LOSAT_X_BXSEEDBUCKET=1` as well, which is not in the
  set (off unless set to any value; no gain measured), so BLASTX keeps the reference order by default.

Rebuilding less:

- `LOSAT_X_SEGSHARE=1` blastp with 2+ threads: the subject SEG masks shared by the whole search,
  as in the one-thread path. Shadow: `LOSAT_X_SEGSHARESHADOW=1`.
- `LOSAT_X_TBNQSIDE=1` tblastn: query SEG and lookup table once per query set.
- `LOSAT_X_TBNSSIDE=1` tblastn: subject six-frame translations once per search of four or more
  query batches (20,000 query residues per batch), up to 256 MiB of translated residues. A
  search of one to three batches keeps none and translates a subject again for each batch, as
  NCBI does (`blast_engine.c:627-631`): one batch has no later reader, and two or three give
  the time the kept translations would save for their memory. Timed with two and three batches
  at 4 and 8 threads, with other work on the machine: on 28,036 16S subjects (40.6 Mb) the
  search without them has 115 MiB less (max RSS 218-226 -> 103-113 MiB) and takes a median of
  1.03-1.11 of the time, on one 606 kb genome the same memory and 1.05-1.08 of the time
  (0.01-0.06 s). Four batches at 8 threads gave 1.06 on the 16S set, so the saving does not
  start at four batches: the number is where the memory is preferred to it. `TBNBATCH` needs
  the kept translations, so two or three batches with at least one batch per thread are now
  searched one after another instead of side by side: 1.02-1.07 of the time on the 16S set and
  1.12-1.24 on the genome (0.04-0.09 s) with one batch per thread, 0.75-0.86 with three
  batches on two threads.
- `LOSAT_X_TBNBATCHRULE=1` tblastn, with `LOSAT_X_TBNBATCH` (speed; in the strict set since
  2026-10-16, on by default, `=0` turns it off; `tblastn/stage_d_pipeline.rs`
  `x_batches_side_by_side`): which query batches `TBNBATCH` searches side by side. With `=0`
  every batch goes side by side once there is one batch per thread (`LOSAT_X_DEFAULT=0` alone
  turns `TBNBATCH` off too, so nothing goes side by side). By default, only the batches of the
  full waves (a multiple of the threads) go side by side and the others one after another with the pool inside each, and
  only when the subjects hold at most 24 Mb or there are fewer than two subjects per thread
  (otherwise every batch one after another). The report does not change (every batch is searched
  as before, the reports keep the batch order). With `LOSAT_WASI_THREADS_DEBUG=1` a line
  `[losat-thread-stage] program=tblastn stage=query_batches batches= threads= side_by_side=
  one_after_another= own_pools=` gives the split (printed whenever `TBNBATCH`'s block runs, with
  or without the rule; `own_pools` is 0 unless `LOSAT_X_TBNOWNPOOLS` is set). Timed on the workstation (perf-si-20261012 session 17, node 1, one hash per
  case; time of each arm / time of main's side by side): 16S (28,036 records, 40.6 Mb) x 10
  batches at 8 threads 0.68 (13.2 -> 8.9 s), 5 batches at 4 threads 0.69 with
  `LOSAT_X_TBNBATCH=0`; one 0.6 Mb genome x 10 batches at 8 threads 0.77 with the rule (8 side
  by side, 2 after), 1.07 with `TBNBATCH=0`; three E. coli genomes (15.7 Mb) x 16 batches at 12
  threads 0.80 with the rule; three fly chromosomes (76.9 Mb) x 16 batches at 16 threads 1.29
  with `TBNBATCH=0` (the rule keeps them side by side), the same bases in 79 pieces 0.95.
  Where all batches are a multiple of the threads and the set is small (0.3-15.7 Mb, one to
  5,000 subjects) side by side took 0.26-0.96 of the time of one after another.
- `LOSAT_X_TBNOWNPOOLS=1` tblastn, with `LOSAT_X_TBNBATCH` and `LOSAT_X_TBNBATCHRULE` (speed;
  not in the strict set: off unless set, to any value; `tblastn/stage_d_pipeline.rs`
  `x_batches_on_own_pools`): the query batches `TBNBATCHRULE` leaves one after another (all of
  them when there are fewer batches than threads, else those after the full waves) are searched
  at the same time instead, each on a search pool of its own with its share of the threads
  (`threads / k`, the first `threads % k` one more), after the shared pool has ended. Only with
  fewer than two subjects per thread, and for two to `threads - 1` batches when the subjects
  hold at most 24 Mb (a pool may then have one thread), otherwise two to `threads / 2`. The
  report does not change (every batch is the same search, the reports keep the batch order;
  `tests/tblastn_batch_rule.rs`); the `query_batches` line above gives `own_pools=`, and one
  `[losat-thread-pool]` line is printed per pool. Timed on the workstation (perf-si-20261012
  session 19, node 0, one hash per case, outputs equal; time with the switch / main's default):
  one 4.7 Mb genome (E. coli MG1655) x 16 batches at 32 threads 0.20 (`-subject`, 3.86 -> 0.79 s)
  and 0.19 (`-db`), x 10 batches at 32 threads 0.24, x 5 batches at 32 / 16 / 8 threads 0.27 /
  0.45 / 0.97, x 2 batches at 8 threads 0.89; one 0.6 Mb genome x 16 batches at 32 threads 0.15,
  x 5 batches at 16 / 8 threads 0.45 / 0.61; three fly chromosomes (76.9 Mb) x 16 batches at 32
  threads 0.81, x 5 batches at 16 threads 0.89. Left out by the conditions, as measured with all
  of them on pools of their own: three fly chromosomes x 5 batches at 8 threads 1.44 (two pools
  of one thread), 5,000 16S records (7.7 Mb) x 5 batches at 8 / 16 threads 1.38 / 1.09, 28,036
  16S records x 5 batches at 8 / 16 threads 1.67 / 0.98, the fly transcriptome x 5 batches at
  16 threads 1.16. Without the switch, 16 batches on the 4.7 Mb genome at 32 threads took 1.60
  (`-subject`) and 2.86 (`-db`) times NCBI BLAST+ 2.17.0 `-db` at 32 threads (1.59 s).
- `LOSAT_X_LUTSPLIT=1` tblastn: build only the half (contexts or table) a caller uses.
- `LOSAT_X_IDEALMEMO=1` BLOSUM62 ideal Karlin block once per process.
- `LOSAT_X_BXSEGMEMO=1` blastx redo stage: the subject SEG once per subject instead of once per window
  (per-thread memo of the last subject, byte-equal input; `blastx/x_seg_memo.rs`). Shadow:
  `LOSAT_X_BXSEGMEMOSHADOW=1` (the masker also runs on every call). Strict-set candidate (S-C, round 5).
- `LOSAT_X_HSPINPLACE=1` (TBLASTX, memory; in the strict set) the HSP list of a subject is sorted, cut
  into linking groups, put together and replayed into NCBI's output order without a second whole copy
  alive (`tblastx/x_hsp_inplace.rs`, exactness argument in its header). The index sorts of
  `sort_ungapped_hits_by_index_replay` (score order, link order) move the records inside the list, cycle by
  cycle, after the same `sort_by`; the linking call cuts its groups off the end of the sorted list as exactly
  sized lists, links them with the same per-group functions (serially or on the pool, as before), joins them
  in one allocation of the exact size, and applies the output replay as a permutation in place. A replay that
  leaves an HSP out or repeats one returns the reference's cloned list (`LOSAT_X_STATS=1` counts both:
  `[X_STATS] HSPINPLACE sorts= linkings= hits= groups= replay_inplace= replay_fallbacks= gate_reference=`).
  The per-group scratch of the literal port is not reserved unless that port can run (`LOSAT_LINKING_LEGACY`,
  a linking trace, `LOSAT_LINK_FAST_SHADOW`). A linking call with `LOSAT_DEBUG_CUTOFFS` or
  `LOSAT_DIAGNOSTICS` runs the reference function (`gate_reference`).
- `LOSAT_X_HSPINPLACESHADOW=1` shadow of `LOSAT_X_HSPINPLACE` (implies nothing else): every sort call and
  every linking call runs both ways (the reference on a clone of the input) and the results are compared HSP
  by HSP, every field, E-values by bits; the first difference stops the run and the search continues with the
  reference's result. One line at exit:
  `[X_SHADOW] HSPINPLACE sorts= linkings= hits_compared= replay_fallbacks=`.
- `LOSAT_X_BPLUT=1` (BLASTP, memory; in the strict set since 2026-10-14; `=0` turns it off) the lookup table
  of the default (non-compressed) BLASTP lookup is
  filled in place from per-cell counts by `LOSAT_X_BXLUT`'s builder (same append order as
  `s_AddNeighboringWords` and `BlastAaLookupFinalize`; `tblastx/lookup/x_lut_direct.rs`, exactness
  argument in its header), so the per-cell chains of the reference builder, about 1.5 times the
  finished table, are never held. In diagnostics mode the reference builder runs. Shadow:
  `LOSAT_X_BPLUTSHADOW=1` (both builders, every table compared field by field, the reference's table
  used); one line at exit: `[X_SHADOW] BPLUT tables_compared=`.
- `LOSAT_X_BPREDOFIT=1` (BLASTP, memory; in the strict set since 2026-10-14; `=0` turns it off) the
  composition-based redo builds each redone HSP list with room for
  exactly its alignments, and a list from which HSPs were removed afterwards gives back their room
  (`blastp/x_redo_fit.rs`): without it the list's vector grows by `push`, whose first allocation has room for
  four HSPs of 168 bytes, while most lists hold one. The lists hold the same HSPs in the same order, so there
  is no shadow. With `LOSAT_TIMING=1` each search of a query set prints
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
- `LOSAT_X_BPSUBJRUNSHADOW=1` shadow of `LOSAT_X_BPSUBJRUN` (implies it): before each later subject of a run
  the run's diagonal offset is compared with the one-thread order's and every table entry is checked to read
  as a cleared one at that offset (flag 0: `last_hit <= offset - window`; flag 1: `last_hit <= offset`), which
  a wrong offset, a reset that does not clear the table or an entry left too high fails; then each subject
  is searched again with a second scratch of the task from a cleared table at the same offset (the loop
  without the switch) and the two results are compared (every HSP list of every query, every field, by their
  `Debug` text). The first difference stops the run (panic) after printing the line that is otherwise printed
  at exit: `[X_SHADOW] BPSUBJRUN searches= runs= subjects= fills_saved= resets_in_runs= subjects_compared=
  tables_checked= identical` (`fills_saved`: subjects searched without a fill of their own; `resets_in_runs`:
  INT4_MAX / 4 resets inside a run; `tables_checked`: tables checked before a later subject of a run).
- `LOSAT_X_BPDBCUT=1` (BLASTP, speed; not in the strict set: off unless set, to any value, and
  `LOSAT_X_DEFAULT` does not touch it; acts only where a `-db` search runs its preliminary stage on the pool:
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
- `LOSAT_X_BPDBCUTSHADOW=1` shadow of `LOSAT_X_BPDBCUT` (implies it): after each `-db` search with the cut,
  the search runs again with blocks of 1 MiB (the same pool, rounds and kernel) into separate hit lists, and
  the two are compared (every preliminary and compact hit list of every query, every field, by their `Debug`
  text; or the texts of the two errors). The first difference stops the run (panic) after printing the line
  that is otherwise printed at exit: `[X_SHADOW] BPDBCUT searches= blocks= blocks_ref= lists_compared=
  identical` (`blocks`: blocks of the cut; `blocks_ref`: blocks of 1 MiB; `lists_compared`: hit lists
  compared that are not empty on either side). The preliminary stage runs twice, so its `LOSAT_TIMING`
  counts double.
- `LOSAT_X_BPSUBJMETA=1` (BLASTP, speed; not in the strict set: off unless set, to any value, and
  `LOSAT_X_DEFAULT` does not touch it; acts only in the `-subject` block search of `LOSAT_X_BPQBLOCK`, for
  blastp and blastp-fast) the subjects' tabular IDs and titles are made once for all blocks instead of once
  per block (`blastp/blast_engine/x_subject_meta.rs`, exactness argument in its header). Without it each
  block's search makes one title `Arc` per subject before its lookup table and every subject's tabular ID
  (`tabular_subject_ids`: the title's first word, or for a `Subject_` first word the first word of NCBI's
  generated definition line) in its output stage, and drops them at its end, on the block's thread while the
  pool's other threads wait. With it the block search's subjects (`x_query_blocks::XBlockSubjects`, made with
  the subjects' encoding) also hold each subject's ID result (`None` where the ID cannot be written: not an
  error until a row needs it), the shown IDs and the titles, and every block borrows them; the check that
  stops a block whose rows name a subject without a writable ID runs over that block's hit lists as before,
  with the stored results (in the block search it never fails: the block search is taken only when every
  subject's ID can be written). The values are functions of each subject record only, so reports, errors
  and their order do not change. The search of all queries, `-db` (no subject records; its titles are read
  per block from the database) and one-thread or many-thread searches without `LOSAT_X_BPQBLOCK` are not
  changed. Memory: the IDs and titles of the subject set, once for the search (they were made per block). They now live
  for the whole block search, alongside the search's peak, about 80 B per subject (about 0.5 MiB for 6,600
  subjects, 8 MiB per 100k), and are freed when the block search ends.
  Measured on the workstation (16 threads, 4 alternating pairs, `LOSAT_X_BPSUBJRUN` and the other blastp
  switches on both sides, same outputs): human x yeast blastp-fast `-subject` with `LOSAT_X_BPQBLOCK=50000`
  (254 blocks) 7.82 -> 7.45 s (0.952), with `LOSAT_X_BPQBLOCK=20000` (633 blocks) 10.58 -> 9.81 s (0.927);
  exact blastp of 3,000 human proteins (23 blocks) unchanged (6.09 / 6.11 s); max RSS within the runs' spread.
- `LOSAT_X_BPSUBJMETASHADOW=1` shadow of `LOSAT_X_BPSUBJMETA` (implies it): in each block whose output takes
  the stored values, each subject's ID result and title and the result of `tabular_subject_ids` over the
  block's hit lists (the shown IDs, or the error's text) are made again as without the switch and compared
  element for element with the stored ones. The first difference stops the run (panic) after printing the
  line that is otherwise printed at exit: `[X_SHADOW] BPSUBJMETA blocks= subjects_compared= identical`
  (`blocks`: blocks whose output took the stored values; `subjects_compared`: the subjects compared in them;
  both 0 where the block search of `-subject` is not taken).
- `LOSAT_X_BPREDOPAR=1` (BLASTP, speed; not in the strict set: off unless set, to any value, and
  `LOSAT_X_DEFAULT` does not touch it) the composition-based redo of a search with more than one thread and
  at least two matches splits its matches between the threads, as NCBI's `omp for` over its array of
  matches, in jobs of up to four consecutive matches of one query, instead of one job per query (or the
  matches of the only query with matches); the early termination test, the heap insertions and the error
  that ends the redo are then replayed on the calling thread in the order of the one-thread loop, on the
  same heaps (`blastp/x_redo_par.rs`, exactness argument in its header). The jobs run in windows of
  consecutive jobs, one `pool.install` each, that may end inside a query: at most 64 jobs per thread
  (`LOSAT_X_BPREDOPARWINDOW` below), so at most 4 x 64 x threads redone results are held before their
  replay. A job builds its query's workspace once per window (`OnceLock`; a query that spans windows has it
  built again from the same inputs, a small cost that grows with the query's length) and its own
  parameters, prefetches the first Newton problems of its own matches with `LOSAT_X_NEWTONLANES` (one vector
  of four, against the heaps as the window's start left them), and redoes its matches with per-worker
  scratch and the search's `LOSAT_X_SEGSHARE` store. A
  query's preliminary list is released after the replay of the window that holds its last match, as the
  one-thread loop releases it (`LOSAT_X_BPPFREE`, `LOSAT_X_BPPCOMPACT`). The early termination test is
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
  `LOSAT_X_NEWTONLANES8` (aarch64) a job's prefetch of at most four problems never fills eight lanes.
  Reason: the query-level redo waits for the largest query of a search (in the blocks of the human x yeast
  blastp-fast job with `LOSAT_X_BPQBLOCK=20000`, 33 % of a block's redo work at the median, 54 % at p90).
- `LOSAT_X_BPREDOPARWINDOW=<jobs per thread>` (tuning variable of `LOSAT_X_BPREDOPAR`; not in the strict set,
  read once when the path first runs) the job budget of a window per thread, a positive integer. Unset or
  any other value (0, a negative number, not a number) gives 64, without a message, as the value switches of
  `utils/xswitch.rs` give their built-in value. The tests set 1 to make windows that end inside queries. The
  output does not depend on it.
- `LOSAT_X_BPREDOLPT=1` (BLASTP, speed; not in the strict set: off unless set, to any value, and
  `LOSAT_X_DEFAULT` does not touch it; acts only where `LOSAT_X_BPREDOPAR` runs, its shadow mode included,
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
  one thread have none). `LOSAT_X_BPREDOPARSHADOW` checks it (every match the replay redoes against the
  one-thread computation, every hit list against the one-thread loop's), and then the shadow line ends with
  ` lpt=1` after `diffs=`; that says only that the switch is set (the fields before it and their order do
  not change; in the repository only `tests/blastp_redo_par.rs` reads that line, by field name), and
  `lpt_windows=` on the `[TIMING]` line counts the windows that ran with it. Memory: one scratch per task
  instead of one per rayon work unit.
- `LOSAT_X_BPREDOPARSHADOW=1` shadow of `LOSAT_X_BPREDOPAR` (implies it): the new path runs first into
  separate heaps without releasing the preliminary lists, and every match that its replay redoes is
  computed again on the calling thread with one set of scratch, as the one-thread loop holds it, and the
  two results are compared (error or not, best score, E-value bits, every HSP field); then the one-thread
  loop runs for the output and the hit lists of both sets of heaps (`fill_results_from_compo_heaps`) are
  compared field by field, floating-point values by their bits. When the one-thread loop ends with an error
  instead, its error is compared with the new path's (both must end, with the same message). The first
  difference stops the run (panic) after printing the line that is otherwise printed at exit:
  `[X_SHADOW] BPREDOPAR searches= queries= matches_compared= lists_compared= diffs=`. With `LOSAT_TIMING=1`
  the new path's `[TIMING] x_bpredopar:` line is printed as well.
- `LOSAT_X_BPSEED=1` (BLASTP, speed; in the strict set since 2026-10-14; `=0` turns it off) the two-hit seed
  stage with an amino-acid lookup table
  (`blastp/x_seed_fused.rs`, exactness argument in its header): the subject scan and the diagonal test in one
  loop that reads each lookup cell's list in place, in the reference's order, and hands only the hits that need
  the context check and the ungapped extension to the reference's code; the context of a candidate comes from a
  two-level table (a coarse table of at most 16,384 entries and a walk over the frame bases) instead of the
  binary search; the diagonal entry of a later hit of the cell is prefetched `LOSAT_X_BPSEED_PREFETCH` hits
  ahead (default 16, 0 = none; tuning only). Measured 2026-10-13 on the workstation (OrthoFinder-style job,
  fly x human, every 10th protein, one thread, medians of 6 rotating rounds under load): 51.4 s -> 34.2 s,
  output unchanged. The switch also covers the compressed lookup table of blastp-fast
  (`blastp/x_seed_fused_compressed.rs`, exactness argument in its header): the scan (NCBI's
  `s_BlastCompressedAaScanSubject` statements) lists the subject offsets whose word passes the PV test, then
  reads those cells in scan order, the backbone cell 8 cells ahead prefetched, and hands each cell's hits in
  NCBI's order (the backbone's offsets, then the overflow list from the cell filled last) to the same test and
  closure; the offset-pair array's rounds are not repeated (a resumed round reads the index the rolling scan has,
  as NCBI's reciprocal division is exact below the table size). Measured 2026-10-14 on the WSL machine (blastp-fast,
  2M query letters x 1M subject letters, `LOSAT_X_BPQBLOCK=20000` on both sides, one thread, 8 alternating pairs,
  load 7.5-10.9): `blastp_scan_ungapped` summed over the blocks 2.05 -> 1.40 s (median ratio 0.67), wall
  8.94 -> 8.22 s (0.89), output unchanged.
- `LOSAT_X_BPSEEDSHADOW=1` shadow of `LOSAT_X_BPSEED` (turns the fused loop on), for both lookup tables: for
  every subject the fused loop and the reference loop run from the same diagonal table, and the initial HSPs
  (every field) and the whole table must be equal; the first difference stops the run (panic). For a subject
  shorter than a word less one letter (amino-acid table, one residue with word size 3) only the initial HSPs are
  compared: the reference's only word there ends one byte past the subject's buffer (NCBI reads the same byte),
  whose value can name a cell with hits; such a hit is never extended and leaves an entry that later subjects read
  as no previous hit (`x_shadow_compare`). The compressed scan reads letters past the buffer as the sentinel in
  both loops, so its tables are always compared.
  One line at exit: `[X_SHADOW] BPSEED subjects= initial_hsps= identical`, with
  `(tables not compared for N subjects shorter than a word)` when there were such subjects.
- `LOSAT_X_BPCLUT=1` (BLASTP `-task blastp-fast`, speed and memory; not in the strict set: off unless set, to
  any value, and `LOSAT_X_DEFAULT` does not touch it) the compressed lookup table (`CompressedAaLookupTable`) is
  built flat in three passes (`tblastx/lookup/x_clut_direct.rs`, exactness argument in its header): a count pass
  makes the reference builder's appends in its order (the words of `add_word_hits`, the neighbors' last letter
  taken as a prefix of the sorted row) and records the cell of each append; a layout pass places each cell's hits
  in one slice of one array (an array of starts, and the reference's PV array and longest chain from the counts);
  a fill pass writes the hits in the order of their appends and then moves the hits of each cell of more than 5
  to NCBI's scan order (`x_emission_position`: the first three, the overflow cell filled last, then the older
  cells newest first). The test for NCBI's overflow banks (the per-batch refusal) is the reference's, with the
  same message. Each block of `LOSAT_X_BPQBLOCK` builds its own flat table. The scans read a cell as its slice:
  the offset-pair scan (`scan_subject`'s statements) when `LOSAT_X_BPSEED` is off, the fused loop
  (`blastp/x_seed_fused_clut.rs`) when it is on. The table holds 4 bytes per hit plus 4 per backbone cell (the
  reference: 32 bytes per backbone cell plus 8 per overflowed hit), and its build holds another 4 bytes per hit.
  Measured 2026-10-14 on the WSL machine (blastp-fast, 2M query letters x 1M subject letters,
  `LOSAT_X_BPQBLOCK=20000` and `LOSAT_X_BPSEED=1` on both sides, one thread, 8 alternating pairs, load 5.9-8.2):
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
  other sessions; same input, `LOSAT_X_BPQBLOCK=20000 LOSAT_X_BPSEED=1 LOSAT_X_BPCLUT=1` on both sides,
  with this change and the per-block set-up on the pool of `LOSAT_X_BPQBLOCK` below; alternating pairs,
  medians; load 15-32): wall at 1 / 4 / 8 / 16 threads 7.36 / 5.00 / 4.55 / 7.92 -> 7.19 / 2.74 / 2.33 /
  3.97 s (median pair ratios 0.98 / 0.55 / 0.51 / 0.49), `blastp_build_lookup` 3.04 / 3.07 / 3.00 / 4.91
  -> 3.03 / 1.03 / 0.86 / 1.41 s, max RSS 39 / 45 / 46 / 53 -> 42 / 42 / 45 / 63 MB, output unchanged.
- `LOSAT_X_BPCLUTSHADOW=1` shadow of `LOSAT_X_BPCLUT` (for blastp-fast): both builders run, the refusal and the
  tables are compared (every cell's hits in the reference's scan order, the PV array, `longest_chain` and the
  scan's parameters) and the first difference stops the run (panic); the search uses the reference's table.
  With more than one thread the pool's build runs too and its outcome and table (every field, slot for slot)
  are compared with the serial flat build's. One line at exit: `[X_SHADOW] BPCLUT tables_compared=
  refusals_compared= hits= parallel_tables_compared= parallel_refusals_compared= identical`.
- `LOSAT_X_BPLKPFLAT=1` (BLASTP `-task blastp-fast`, speed; not in the strict set: off unless set, to any value,
  and `LOSAT_X_DEFAULT` does not touch it; acts only where the search's pool builds the flat table of
  `LOSAT_X_BPCLUT` or `LOSAT_X_BPCLUTSHADOW`, with more than one thread) changes how that build is scheduled
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
- `LOSAT_X_BPLKPFLATSHADOW=1` shadow of `LOSAT_X_BPLKPFLAT` (implies it): after each build of the switch (the
  table the search uses, on the recycled vectors) the table is built again with the `par_iter` passes on fresh
  zeroed vectors, and the outcomes (table, refusal, `TooLarge`) and the tables (every field, slot for slot) are
  compared; the switch's build fills its vectors with `u32::MAX` / `i32::MIN` first, so that an entry no
  stripe writes always differs from the fresh build's; the contexts set after a build with Karlin tasks are compared with those of
  `x_query_blocks::x_bp_lookup_query` (the query, the segments, every context field, the Karlin blocks bit for
  bit). The first difference stops the run (panic). One line at exit: `[X_SHADOW] BPLKPFLAT tables_compared=
  outcomes_compared= recycled= recycled_hits= context_sets_compared= contexts_compared= contexts_from_tasks=
  identical` (`recycled` / `recycled_hits`: builds whose `starts` / whose `hits` are the buffer of an earlier
  table, `hits` only where it had the room; with `LOSAT_X_BPLKPAHEAD` on, also builds ahead on the second
  spare's vectors, which `x_make_room` allocated or moved there, since the reuse test sees only capacity; `contexts_from_tasks`: contexts computed by the builds' tasks). With `LOSAT_X_BPCLUTSHADOW` the pool's build is the switch's (on the vectors that the
  previous block's pool build gave back) and is compared with the serial flat build as well.
- `LOSAT_X_BPLKPAHEAD=1` (BLASTP `-task blastp-fast`, speed; not in the strict set: off unless set, to any value,
  and `LOSAT_X_DEFAULT` does not touch it; acts only in the block search of `LOSAT_X_BPQBLOCK`, with `-subject`
  or `-db`, where `LOSAT_X_BPLKPFLAT` acts on the pool's flat table of `LOSAT_X_BPCLUT`: more than one thread;
  not with `LOSAT_X_BPCLUTSHADOW`, not for exact blastp) builds the next query block's query frames, lookup
  table and contexts during the block before it, as one task on the search's pool, instead of at the start
  of the next block's own search (`blastp/blast_engine/x_lookup_ahead.rs`, exactness argument in its
  header): the block loop runs each block's search in a `rayon::in_place_scope`, and when the block's own
  table is ready (at once when it was itself built ahead) the next block's build is spawned there; the next
  block's search then takes the values instead of computing them. They are the values of the same calls
  (`blastp_query_frames`, `blastp_compressed_lookup`) on the same records, batches and pool, so the reports
  do not change. A build ahead that fails, refuses, would need the reference builder (`TooLarge`) or lacks
  the room in its vectors gives nothing, and the block builds inline, so an error is raised where it is
  without the switch; a panic of the build is kept and resumed where the block's search would start (with
  `panic = "abort"`, the release profile, it ends the process when it happens). Memory: one more table
  (`starts`, `hits`, PV array) and one block's frames and contexts are alive while the next block's table
  is built; every vector of that table is allocated by the block search's thread, never by a pool task:
  that thread keeps a second spare (`x_query_blocks.rs` `XBlockSubjects`) with the room for the largest
  table so far plus a half (new vectors of twice its hits where they are short; untouched capacity is not
  resident), which the build ahead takes whole (`XFlatArgs::x_ahead`). At one thread nothing changes.
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
- `LOSAT_X_BPLKPAHEADSHADOW=1` shadow of `LOSAT_X_BPLKPAHEAD` (implies it): each block that takes values
  built ahead also builds them inline as without the switch (on the block spare's vectors, given back
  after) and compares them: the query frames field for field, the flat tables field for field and slot for
  slot, the contexts with their Karlin blocks bit for bit. The first difference prints the line below with
  `DIFFERENT` and stops the run (panic). One line at exit: `[X_SHADOW] BPLKPAHEAD blocks= prebuilt= dropped=
  compared= identical` (`blocks`: blocks searched where the switch acts; `prebuilt`: blocks that took
  values built ahead; `dropped`: builds ahead that gave none, whose blocks built inline; `compared`: blocks
  compared; every block but the last starts one build once its own table is ready, if that table is a flat
  one, as it is unless the reference builder made it, so `prebuilt + dropped = blocks - 1`).
- `LOSAT_X_BPDIAG=1` (BLASTP, speed; in the strict set since 2026-10-17; `=0` turns it off; acts only in
  the fused loop of `LOSAT_X_BPSEED` or `LOSAT_X_BPSEEDSHADOW`) the two-hit
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
  nothing (`LOSAT_X_BPSEED_PREFETCH` is the per-hit loop's). Shadow: `LOSAT_X_BPSEEDSHADOW=1` with
  `LOSAT_X_BPDIAG=1` compares the batched loop with the reference for every subject. With `LOSAT_TIMING=1`
  each search of a query set prints `[TIMING] x_bpdiag: mode=<avx2|scalar|per-hit> cells=<n> batched_hits=<m>
  candidates=<k>` (the path that ran, cells tested in batches, their hits, the hits passed on; zeros with
  `per-hit`). Measured 2026-10-14 on the workstation (one pinned core, alternating, the other 8 exact switches
  and `LOSAT_X_BPQBLOCK=50000`, x human proteome by `-subject`, rows identical in every arm; `blastp_scan_ungapped`
  (wall) in s): human block 1 per-hit 11.59 (20.61), gather 10.17 (18.90), scalar loads 9.13 (17.86), scalar
  batch 12.56 (21.53); C. elegans block 5 per-hit 32.68 (70.68), gather 28.55 (66.16), scalar loads 26.45
  (64.62).
- `LOSAT_X_BPUNGAP=1` (BLASTP and BLASTX one-hit and two-hit BLOSUM62 ungapped extension, speed; in the
  strict set since 2026-10-14; `=0` turns it off) the word loop and the
  X-drop walks of `s_BlastAaExtendLeft` / `s_BlastAaExtendRight` read the scores from one flat 32 x 32 byte
  table (`q << 5 | s`; the sentinel's `BLAST_SCORE_MIN` stored as -128) with raw-pointer walks over the
  positions the reference indexes, NCBI's statements and exit tests unchanged (`blastp/x_ungapped_fast.rs`,
  exactness argument in its header). Acts in `extend_one_hit_blosum62` / `extend_two_hit_blosum62` (TBLASTN
  calls the matrix-generic functions and is not affected); a dropoff outside 1..=127, a word size outside
  1..=11 or a residue code >= 32 run the reference.
- `LOSAT_X_BPUNGAPSHADOW=1` shadow of `LOSAT_X_BPUNGAP` (turns the path on): every extension runs both ways and
  the results (score, offsets, length, `right_extend`, `s_last_off`, and the early return) must be equal; any
  difference panics with the inputs. One line at exit: `[X_SHADOW] BPUNGAP compared=`.
- `LOSAT_X_LINKSORT=1` (TBLASTX, speed; in the strict set since 2026-10-14) the index sorts of the
  in-place HSP list path (`tblastx/x_link_sort.rs`, exactness argument in its header; it acts only where
  `LOSAT_X_HSPINPLACE` does).
  Two changes, both of how an index order is computed, not of the order: (1) the linking call sorts its
  output order once, stably, with `s_FwdCompareHSPsTransl` and ties by `s_RevCompareHSPsTransl`, where it
  sorts twice (a stable sort by B followed by a stable sort by A is one stable sort by "A, then B"); this
  also saves a sort at one thread. (2) An index sort of 32,768 indices or more runs as rayon's stable
  parallel sort (`par_sort_by`) when its call may use the search pool: the two sorts of a linking call
  under the call's `allow_parallel`, as its per-group linking; the other HSP sorts (the score sorts after
  each frame, after the first linking call and after the re-evaluation) when the calling thread is a worker
  of a rayon pool, which in a search started by the command line or the Web adapter it is exactly when
  `-num_threads` is 2 or more (the header names the one case this repository's code does not show: a
  program that calls a one-thread search from a thread of its own rayon pool). A stable sort with a
  comparator that is a total preorder has one result, so the index order is the same at any number of
  threads; the comparators are NCBI's four, each a lexicographic comparison of integer fields. With one thread, in a
  build without the `parallel` feature and on wasm32 without `wasm-threads` every sort is the serial
  statement. With the switch off the statements of `x_hsp_inplace.rs` run as before it existed.
  `LOSAT_X_STATS=1` prints, only with the switch or its shadow on,
  `[X_STATS] LINKSORT sorts= indices= parallel_sorts= parallel_indices= composite_sorts=`
  (index sorts through the path and their lengths added; those that ran on the pool; final sorts made as
  one). With several subjects on a pool, a thread that waits for the parts of its sort can run another
  subject meanwhile, as it can while it waits for its linking groups, so the `sort`, `final_sorts`,
  `merge` and `reap` times of `LOSAT_TIMING` can then hold time of another subject. Measured once on
  2026-10-13 (medians of three alternating runs without and with the switch, one binary, a loaded
  machine; `[TIMING] x_link1` and `x_subject`, seconds, without -> with): AP027131 x AP027133 (genetic
  code 4, 338,859 HSPs in the first linking call) at 8 threads `sort` 0.085 -> 0.039, `final_sorts`
  0.067 -> 0.009, `reap` 0.073 -> 0.045, `merge` 0.167 -> 0.156, and at one thread `final_sorts`
  0.060 -> 0.035 (the other sorts are today's statements there); NZ_CP006932 against itself (516,395
  HSPs) at 8 threads `sort` 0.139 -> 0.075, `final_sorts` 0.100 -> 0.014, `reap` 0.098 -> 0.068, and at
  one thread `final_sorts` 0.105 -> 0.047. Max RSS did not change (164 -> 161 and 242 -> 243 MiB at 8
  threads, one run each). What is left of `sort` and `reap` includes the move of the records into the
  sorted order, which stays on one thread. In a search of 200 short subjects no list reached 32,768 HSPs
  and no sort went to the pool.
- `LOSAT_X_LINKSORTSHADOW=1` shadow of `LOSAT_X_LINKSORT` (turns the path on): every index sort runs both
  ways, the path's sort on a copy of the index vector and today's statement or statements on the vector
  itself, and the two vectors must be equal element for element; the first difference stops the run with
  the sort (first, final or dispatch), the length and the position, and the search continues with today's
  result. One line at exit: `[X_SHADOW] LINKSORT sorts= indices_compared=`. Together with
  `LOSAT_X_HSPINPLACESHADOW=1` the reference bodies are compared with the path's result end to end
  (`LOSAT/tests/tblastx_link_sort.rs` runs both).
- `LOSAT_X_LINKSLICES=1` (TBLASTX, memory; in the strict set since 2026-10-15) the linking call of
  the in-place HSP list path links its groups as disjoint sub-slices of the one sorted list
  (`tblastx/x_link_slices.rs`, exactness argument in its header; it acts only where
  `LOSAT_X_HSPINPLACE` does). Without it the call cuts the sorted list into
  one Vec per group, links each and copies the linked groups into a new list of the whole size, so the
  list exists twice while the new one is filled (AP027131 x AP027133, genetic code 4: 49 MiB twice in
  the first linking call, where the heap of that search peaks). With it the groups are the same runs
  (the same start test, `group_starts`), each kernel writes fields of its own records in place and
  changes no order, and the list after the kernels is today's assembled list, record for record; no
  second list is made. The kernels, their arguments, the per-thread scratch and the pool condition
  (`allow_parallel` and more than one group) are today's; `LOSAT_LINKING_LEGACY` and the linking traces
  link each group's copy with the literal port and write it back, since that port returns its group
  in the same order and length. With the switch off the statements of `x_hsp_inplace.rs` run as
  before it existed, after one read of the mode per linking call. `LOSAT_TIMING=1`: `assemble` of
  the `x_link1` / `x_link2` lines is then about zero. `LOSAT_X_STATS=1` prints, only with the
  switch or its shadow on,
  `[X_STATS] LINKSLICES linkings= groups= parallel_linkings= legacy_fallbacks=` (linking calls
  linked as slices, their groups, those on the pool, and calls handed back to the cut, which never
  happens: always 0). Measured once on 2026-10-13 (one binary, a loaded machine): massif of AP027131 x
  AP027133 at one thread, heap requested, 107.3 MiB at the peak in the first linking call without the
  switch; with it that call reaches 77.3 MiB and the search peaks at 102.4 MiB at the end of the seed
  stage, as without it. Max RSS (medians of three alternating runs, without -> with): that pair 110 ->
  111 MiB at one thread and 162 -> 163 at eight (its RSS peak is in the seed stage); NZ_CP006932
  against itself 142 -> 132 at one thread and 242 -> 222 at eight. `groups` 0.011 -> 0.002 s and
  `assemble` 0.010 -> 0.000 s in the first linking call of the pair at one thread.
- `LOSAT_X_LINKSLICESSHADOW=1` shadow of `LOSAT_X_LINKSLICES` (turns the path on): every linking call's
  group step runs both ways, today's cut and assembly on the list and the slices on a copy of the
  sorted list, and the two lists are compared record by record, every field, E-values by bits; the
  first difference stops the run with a message that starts with the variable's name, and the search
  continues with today's list. One line at exit: `[X_SHADOW] LINKSLICES linkings= hits_compared=`
  (`LOSAT/tests/tblastx_link_slices.rs` runs it with the other linking variables).
- `LOSAT_X_BPQBLOCK=<letters>` (BLASTP, memory and speed; in the strict set since 2026-10-14, a value switch:
  a positive number of letters, built-in 50,000; `=0` turns it off) the queries are searched in blocks of whole query batches
  (`blastp/x_query_blocks.rs`, exactness argument in its header): a block closes at the first batch end where
  its letters reach the value, and each block is a search of its own (its lookup table, its diagonal table,
  its hit lists) whose rows are written before the next block is searched, as NCBI searches and reports one
  batch at a time. Both lookup tables take it: the default one, and the compressed one of `-task blastp-fast`
  (word size 5), whose table each block builds for its own batches and tests for NCBI's overflow banks with
  them. Taken only where a block's report is its rows: one output in tabular format 6, no `-query_loc`, every
  batch read, more than one block, no subject whose tabular ID the reference path rejects, and none of the
  reference path's other refusals that come before its rows (a `-window_size` whose diagonal table with all
  queries passes 2^30, a query split that LOSAT rejects, or with the compressed table a query batch whose
  offsets add more than four times NCBI's overflow cells, about 859 million words, counted before the search);
  otherwise the search of all queries. With `-db` the same, and the blocks run only when every definition line
  of the database reads and none of LOSAT's refusals of a reported title applies (`x_titles_allow_blocks`, one
  pass over the header files in the search pool before the first block); a query whose effective search space
  could pass NCBI's `Int8` (with `-db` or `-dbsize`) or an empty database also takes the search of all queries.
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
- `LOSAT_X_BPPFREE=1` (BLASTP, memory; in the strict set since 2026-10-14; `=0` turns it off) the preliminary
  hit lists hold no more than their HSPs need
  (`blastp/x_prelim_free.rs`, exactness argument in its header): a subject's HSP list for a query gives back
  the room of the HSPs that the preliminary E-value reaped, and the serial composition-based redo (one thread)
  frees each query's preliminary list after its redo, as NCBI's redo frees each list it has redone. The
  parallel redo paths keep the lists, as before. No value changes, so there is no shadow. With
  `LOSAT_TIMING=1` each search of a query set prints `[TIMING] x_bppfree: freed_lists=<n>`. Measured once
  (C. elegans block of 150 proteins, 50,869 letters, x human proteome, one thread, max RSS): 153,744 ->
  117,968 kB (NCBI -db 109,080 in the same run).
- `LOSAT_X_BPPCOMPACT=1` (BLASTP, memory; in the strict set since 2026-10-14; `=0` turns it off) the subject
  scan keeps each preliminary HSP as a 64-byte element (the
  `BlastpPreliminaryHsp` and its E-value) instead of a 168-byte `BlastpHsp`, in the same hit-list helpers
  made generic over the element (`blastp/x_prelim_compact.rs`, exactness argument in its header); the redo
  reads lists rebuilt as the off path holds them: the one-thread redo rebuilds each query's list when it
  reaches it and drops it after its redo, a search with more than one thread rebuilds all lists before the
  redo, query chunks hand on rebuilt lists; with `-db` the compact lists are listed and renamed to ranks with the
  others (`ListedSubjects::read`). `LOSAT_X_BPPCOMPACT=shadow` builds both lists, compares every
  rebuilt list field by field with the off path's (it panics on a difference) and prints
  `[X_SHADOW] BPPCOMPACT compared=<lists> hsps=<n>` at exit. With `LOSAT_TIMING=1` each search of a query
  set prints `[TIMING] x_bppcompact: expanded_lists=<n>`.
- `LOSAT_X_BNSTREAM=1` (BLASTN, all tasks, `-subject` with two or more threads, memory; in the strict set
  since 2026-10-16, on by default) the threads claim the subjects
  of a query batch one at a time in OID order and one of them (the owner) hands the results to the hit lists in
  OID order while the others go on with the next subjects; a thread claims a subject only while it is fewer
  than `threads x 64` subjects past the last one handed over, so at most that many subjects' results wait
  (`x_bnstream::search_in_order`, called from `search_query_batch` in `blastn/blast_engine/run.rs`; each
  thread runs its loop as one `rayon::broadcast` job, so that a thread waiting inside a nested parallel stage
  of a search never starts a second loop), in both passes: the preliminary HSPs to the preliminary hit lists
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
- `LOSAT_X_BNSTORE=1` (BLASTN, all tasks, `-subject` without `-subject_loc`, memory; not in the strict set: off
  unless set, to any value, and `LOSAT_X_DEFAULT` does not touch it) each subject record's letters are packed as
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
- `LOSAT_X_SEEDRELEASE=1` (TBLASTX, memory; in the strict set since 2026-10-13) the tables that only the
  scan reads are freed before the HSP list of the batch's last subject is merged and linked, which is where
  the heap of a large search peaks
  (`tblastx/x_seed_release.rs`, exactness argument in its header). A query batch keeps its tables in two
  holders, each reached only through one stage per subject: the unit path of `LOSAT_X_PAIRPAR` with its group
  tables, whose stage ends when the subject's units have run, before their lists are merged (the unit path
  runs as two calls for this, `run_subject_units` and `merge_subject` in `x_pair_par.rs`, the same statements
  in the same order); and the batch's lookup table, whose stage ends after the seed stage. The end of the
  batch's last stage of a holder frees it. The worker's diagonal table, which only the reference frame loop
  reads, is not created when the unit path takes the batch without its shadow, and is otherwise dropped by
  the worker that ended the last lookup stage. With the switch off the same holders keep the tables until the
  batch has been searched, as before it existed. Measured once on AP027131 x AP027133 (genetic code 4, one thread, bytes
  requested under massif), without -> with the switch: the end of the units 110.4 -> 102.4 MiB, the merge of
  their lists 155.8 -> 81.0 MiB, the first linking call 167.6 -> 107.3 MiB, which is the peak of the run in
  both; max RSS in one run each 167 -> 110 MiB at one thread and 221 -> 163 MiB at eight. Known limits: the
  gain is for the last subject of a batch (the only one of a pairwise search), since earlier subjects are
  merged and linked while later ones still need the tables; and with `LOSAT_X_PAIRLUT=0` the batch's full
  lookup table, which the frame loop names, stays alive through the merge, so that path keeps the higher
  seed-stage peak. `LOSAT_X_STATS=1` prints, only with the switch on,
  `[X_STATS] SEEDRELEASE batches= tables_released= unit_tables_released= batches_without_worker_table= worker_tables_released=`
  (batches searched; batches whose lookup table a stage end released; batches whose unit path a stage end
  released; batches whose workers held no diagonal table; worker tables dropped at the release). There is no
  shadow switch: the switch changes no value and no order of computation, only when three allocations are
  freed, so there is no second result to compare. A read after a release is excluded by ownership instead:
  a stage is the only way to its tables and ending it consumes it (the compiler rejects a later use), the
  merge after the units stage is a function that cannot name a group table, an open stage keeps its tables
  alive, a subject that arrives after a release stops the run, and an assertion that runs with and without
  the switch stops a frame loop that would run without the worker's table
  (`LOSAT/tests/tblastx_seed_release.rs` reads the counters).
- `LOSAT_X_SEEDCOMPACT=1` (TBLASTX, memory; not in the strict set: off unless set, to any value, and
  `LOSAT_X_DEFAULT` does not touch it) the seed stage of the unit path of `LOSAT_X_PAIRPAR` and its merge
  hold fewer and exactly sized lists (`tblastx/x_pair_par.rs`, exactness argument in the section above
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
  scan loop and `LOSAT_X_PAIRPAR=0` are untouched). There is no shadow variable of its own:
  `LOSAT_X_PAIRPARSHADOW=1` with the switch compares the unit path's `combined_ungapped_hits`, made from
  the compact lists, with the reference loop's for every subject. `LOSAT_X_STATS=1` prints, whenever the
  switch is on, `[X_STATS] SEEDCOMPACT subjects= units= records= bytes_saved=` (subjects and units merged
  from compact lists, records expanded, and records x (size of `InitHSP` - 28): the bytes the unit
  outputs no longer request; all 0 where no unit path runs). `LOSAT/tests/tblastx_seed_compact.rs` reads
  the counters. Measured once on 2026-10-17 (WSL machine under load, one binary, max RSS from `wait4`, means
  of three alternating runs without -> with the switch, every report the frozen one): AP027131 x AP027133
  (genetic code 4) 111.0 -> 88.7 MiB at one thread and 164.0 -> 154.5 at eight; NZ_CP006932 against itself
  132.2 -> 123.9 and 226.3 -> 215.1 (at eight threads the peak is then the first linking call). The medians
  of `units`, `merge` (`[TIMING] x_pairpar`) and `search_total` moved within the run-to-run spread of that
  machine; the probe build measured the same memory on the workstation (-22.6, -7.5, -8.1 and -9.4 MiB).

- `LOSAT_X_TBXFRAMES=1` (TBLASTX, memory; in the strict set since 2026-10-17, on by default) a
  subject's frames are made when they are read and only as long
  as they are read (`Subjects::x_preliminary_frames` / `x_letter_frames` in `tblastx/subjects.rs`,
  `x_generate_frames_in_windows` in `tblastx/translation.rs`). Without it `process_subject` makes, before
  the seed stage, the six frames of the subject's letters and, for a subject with ambiguous bases, the six
  frames of the bases the preliminary stage translates, each from a whole copy of its sequence and of its
  reverse complement, and keeps both sets to the end of the subject (a stored sequence with ambiguity codes:
  about 7 bytes per base at once). With it only the preliminary frames are made before the seed stage, in
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

Scheduling (each value computed by one thread with the original function; order of output kept):

- `LOSAT_X_AHEAD=16` blastn/megablast ordered gapped-extension and traceback loops: helpers compute
  the next elements' alignments ahead, the tree is read and written in the original order by one
  thread (`utils/xahead.rs`). Shadow: `LOSAT_X_AHEADSHADOW=1`. `LOSAT_X_GREEDYSPEC=1` adds the
  megablast greedy traceback; `LOSAT_X_SPECBATCH` tunes the batch.
- `LOSAT_X_PAIRPAR=1` (strict-set candidate; S-E, round 5) blastn/megablast, one query-subject
  pair on a pool of 2+ threads. (1) The megablast lookup table is filled by the pool threads: every
  thread reads the whole query in order and applies the update only to the words whose cell lies in
  its blocks of `2^max(11, pv_array_bts)` cells, so each cell gets its words in the original order
  (`blastn/x_mb_lookup_par.rs`; not with database word counts or `BLEMIR_DEBUG`). (2) The one-hit
  seed stage of a subject chunk with the diagonal hash (query block over 8000 letters,
  `window_size` 0, `scan_range` 0, contiguous megablast table, unmasked subject, one subject, no
  small-query word extension, no debug or trace): the subject scan grid is cut into strips that
  the pool threads scan (lookup chains and word extension), and the owner thread folds the strips
  in scan order through a copy of the diagonal-hash code (`blastn/blast_engine/x_pair_seed.rs`).
  Anything else runs the original code. Shadow: `LOSAT_X_PAIRPARSHADOW=1` (lookup tables and, per
  chunk, the hit list and the whole diagonal hash compared with the original code; one
  `[X_PAIRPARSHADOW]` summary line on stderr at the end). With `LOSAT_TIMING` a
  `[TIMING] x_pairpar_seed` line gives the strips, seeds and the owner's fold, scan and wait time.
- `LOSAT_X_DISCLUTPAR=1` (BLASTN `-task dc-megablast`, and megablast given `-template_type`,
  `-template_length` and `-word_size` 11 or 12; `-subject` and `-db`; not in the strict set: off unless set
  to any value, whatever `LOSAT_X_DEFAULT` says) the discontiguous megablast lookup tables
  (`lookup.rs` `build_disc_mb_lookup`, NCBI `s_FillDiscMBTable`: `hashtable`, `next_pos` and the counters
  of the longest chain, and with two templates `hashtable2`, `next_pos2` and their counters, both
  templates sharing one presence vector) are filled by the threads of the search pool
  (`blastn/x_disc_lookup_par.rs`). The words are read in chunks of 524,288 positions, each in three
  phases on the pool: (A) the chunk's positions are cut into one span per thread, and each thread
  reads the words of its span with the loop's statements (entering a range up to `template_length - 1`
  positions early, so the words and their template indices are the loop's) and sorts them into
  buckets by the owner of their cell (blocks of `2^max(11, pv_array_bts)` cells dealt round robin, the
  plan of `LOSAT_X_PAIRPAR`'s megablast table); (B) each owner applies the loop's update to the words of
  its cells span by span, so each cell gets its words in the original order, and keeps each word's
  previous head; (C) each thread writes `next_pos` (and `next_pos2`) of the words of its span from the
  kept heads. Each presence-vector word and counter is written by one owner (a presence bit is only
  ever set), and `longest_chain` is computed afterwards by the original code. The buckets hold one
  8-byte entry per word of the current chunk and template; before each chunk a bucket gives back the
  room it grew above four times its share of a chunk, so the buckets of one template reserve at most
  24 MiB plus 544 bytes x threads^2 whatever the input. Measured (unit test, chunks of 524,288
  positions): 8 MiB per template on random letters at 4, 16 and 64 threads; 16 / 15 / 12 MiB on
  blocks of unmasked short repeats (each chunk's words on the cells of one or two owners, changing
  from chunk to chunk), where buckets that kept their room reserved 16 / 32 / 48 MiB (up to 2 x
  threads x 4 MiB); on (a) below 6.0 MiB at 4 threads and 7.8 MiB at 16. It acts only where a LOSAT search pool
  of two threads or more is current: the search's pool, or the private pool of a `LOSAT_X_QBATCHPAR`
  slot of two or more threads (never rayon's global pool); with one thread, a slot of one thread, a
  build without the `parallel` feature or WebAssembly without `wasm-threads` the original loop runs.
  Shadow: `LOSAT_X_DISCLUTPARSHADOW=1` (implies the switch) builds the tables again the other way
  (serially when the pool built them, otherwise in four parts one after another) and compares every
  element and the longest chain; a difference prints `[X_SHADOW] DISCLUTPAR ... DIFFERENT ...` and
  panics; at exit `main` prints `[X_SHADOW] DISCLUTPAR tables= pool= two_templates= cells= diffs=0
  identical` (tables compared, of them built on a pool, of them with two templates, elements
  compared). With `LOSAT_TIMING` (any switch state) `main` prints `[TIMING] disc_lookup_build: <s>s
  (builds=<n> parallel=<n>)` when a table was built: the summed wall time of the builds (allocation,
  fill and longest chain, not the shadow's second build; with `LOSAT_X_QBATCHPAR` it includes the
  slots' concurrent builds) and how many were filled on a pool. A second line, `[TIMING]
  disc_lookup_split: alloc= loop= pool_fill= chain= chunks= read_wall= read_max= read_sum= read_late=
  update_... link_... buckets_peak=`, splits the same builds: the allocation, the fill by the original
  loop, the fill on the pool and the longest chain (seconds, summed); for the pool fills (not the
  shadow's) their chunks and, per phase (read = A, update = B, link = C) summed over the chunks, the
  calling worker's wall time, the longest job, the jobs' summed time and the latest job start; and the
  most the buckets of one fill reserved. The jobs' times are taken only with `LOSAT_TIMING`. Measured on the workstation (2x Xeon
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
  the jobs' own work, whose sum is 3.8 times (4 threads) and 4.9 times (16 threads) the loop's whole
  fill; phase B (the owners' updates) alone sums to 2.4 and 3.9 times the loop, with balanced jobs
  (longest / mean 1.04 and 1.08); at 4 threads phase A is also costly (the longest reader 1.29 times
  the mean). The runs with the switch also have more system time (0.22 against 0.15 s at 4 threads,
  0.90 against 0.60 s at 16). What makes a word's work dearer in the phases is not measured;
  candidates are the first touch of the fresh tables by many threads (transparent huge pages are
  `madvise` there, so 4 KiB page faults) and cache lines of the buckets moving between cores (written
  in A, read and rewritten in B, read in C). Fewer joins would save at most the 0.01 s above. Fix round
  1's bucket bound costs nothing measurable: (a) at 16 threads, 4 rounds old (24c714fd) / new
  (b1e3df69) / off, all slow phase: walls 2.056 / 2.055 / 2.174 s (new/old 1.000, 0.999-1.001; new/off
  0.945), max RSS 91.3-95.3 / 91.3-94.3 / 88.9-92.9 MiB; (b) at 16 threads one round in a common phase
  4.419 / 4.419 s (the other mixed the phases). Every output equals the sha256 of the run without the
  switch. Tests: the unit tests of
  `x_disc_lookup_par.rs` (partitioned and pool fills equal the original loop for every template type,
  length and word size, one and two templates, 1-20 parts, chunks of 1 to 5,000 positions and random
  layouts of queries and masks, pools of 1, 2, 3, 5, 8 and 16 threads, a 2.6 M-position query in the
  default chunks on pools; the template indices read only the last `template_length` bases; the owners
  run on at least two workers of a 3-thread pool; the shadow names a change in any one array or the
  longest chain; the buckets' bound) and `LOSAT/tests/blastn_disc_lookup_par.rs`
  (off / on / shadow equal at 1 and 4 threads, `-subject` and `-db`, one and two templates, with
  `LOSAT_X_PAIRPARDISC`, `LOSAT_X_QBATCHPAR` slots and split query chunks, megablast with and without a
  template, the timing and split lines, and AP027131 against NZ_CP006932 with one and two templates,
  whose tables span two chunks).
- `LOSAT_X_PAIRPARDISC=1` (not in the strict set: off unless set to any value, whatever
  `LOSAT_X_DEFAULT` and `LOSAT_X_PAIRPAR` say) blastn `-task dc-megablast`, one query-subject pair
  on a pool of 2+ threads: the two-hit seed stage of a subject chunk with a discontiguous megablast
  table (one or two templates) and the diagonal hash (query block over 8000 letters, `scan_range`
  0, unmasked subject, one subject, no debug or trace). The pool threads scan strips of the subject
  with the discontiguous scan kernel (both templates' chains in NCBI's order, the query context of
  each pair), and the owner thread folds the pairs in scan order through a copy of the two-hit
  diagonal-hash code (window test, `s_TypeOfWord`, ungapped extension, every insert;
  `blastn/blast_engine/x_pair_seed_disc.rs`). The helpers scan at most two strips per thread
  ahead of the owner and reuse the owner's consumed buffers. The owner's fold is the slower side
  for dc-megablast (most pairs reach the diagonal hash), so the stage stops gaining at about 4
  threads. Anything else (contiguous tables, which `LOSAT_X_PAIRPAR` takes in one-hit mode) runs
  the original code. Shadow: `LOSAT_X_PAIRPARDISCSHADOW=1` (per chunk the hit list and the whole
  diagonal hash compared with the original code; one `[X_PAIRPARDISCSHADOW]` summary line on
  stderr at the end, with the number of chunks compared). With `LOSAT_TIMING` a
  `[TIMING] x_pairpardisc_seed` line gives the strips, seeds and the owner's fold, scan and wait
  time. Test: `LOSAT/tests/blastn_pairpar_disc.rs`.
- `LOSAT_X_QDUSTPAR=1` (BLASTN, all tasks with DUST on, which is the default of megablast, dc-megablast and
  blastn, and blastn-short with `-dust yes`; not in the strict set: off unless set to any value, whatever `LOSAT_X_DEFAULT` says) the DUST masks
  of a query batch are made on the search pool instead of one query after the other
  (`blastn/x_query_dust_par.rs`, called from `coordination.rs` `apply_dust_masking`, so by `query_masks` for
  every query batch and by `chunk_query_masks` for the parts of a split query; `-db` searches too). NCBI makes
  a fresh `CSymDustMasker` for each query; the mask of a query is a function of its letters and the three
  DUST parameters only (`DustMasker` keeps no state across calls, each call makes its own random generator for
  `N`, and the statics of `utils/dust.rs` are switches read once), so the per-query map is
  `par_iter().map().collect()` in query order (rayon splits the queries into ranges; the masks are per query,
  the jobs are not) with `SearchPool::install`: on the pool of the search that calls
  it (the whole search pool, or the private pool of a `LOSAT_X_QBATCHPAR` slot of two or more threads), never
  on rayon's global pool. It acts only with a pool of two or more threads and two or more queries in the
  batch; with one thread, a slot of one thread, a build without the `parallel` feature, WebAssembly without
  `wasm-threads`, or the switch off, the serial map runs (with `-dust no` there is no map at all). The verbose lines and totals stay as
  they were. Shadow: `LOSAT_X_QDUSTPARSHADOW=1` (implies the switch) also makes the masks serially and compares
  them element for element; a difference counts in `diffs`, prints
  `[X_SHADOW] QDUSTPAR batches= queries= intervals= diffs= DIFFERENT` and panics; at the end of a run `main`
  prints `[X_SHADOW] QDUSTPAR batches= queries= intervals= diffs=0 identical` (the batches, queries and masked
  intervals compared, all 0 when no batch had a pool and two queries; the line ends in `DIFFERENT` when
  `diffs` is not 0, which can happen at exit only when a `LOSAT_X_QBATCHPAR` slot's batch differed and the
  loop discarded that batch, in a build that unwinds). With `LOSAT_TIMING` (any switch state) `main` prints
  `[TIMING] query_dust: <s>s`, and with `LOSAT_X_STATS` `[X_STATS] QDUSTPAR dust_s=<s>`, only when at least one
  batch was masked (not for BLASTP, BLASTX, TBLASTN, TBLASTX or `-dust no`). `<s>` is the sum of the wall
  times of every DUST call of the process, not the wall time of the step: with `LOSAT_X_QBATCHPAR` it includes
  the slots' batches, which run concurrently, and those the loop discards; with `LOSAT_X_QDUSTPARSHADOW` it
  includes the serial reference pass. The ws numbers below have neither, so there it is the step's time.
  Measured on the workstation (2x Xeon Gold 6338, powersave, 16 idle physical cores of one node with the
  siblings idle, `numactl --membind`; `-outfmt 6`, subject 3 Mb of D. melanogaster 2L, `LOSAT_X_PAIRPARDISC=1`
  for the dc-megablast inputs, `LOSAT_X_STATS=1` in the timed pairs; off runs carry an unread
  `LOSAT_X_PADQDUSTPAR=1`), 16 threads, alternating off/on pairs, medians: (a) 1,775 cDNAs (5.0 Mb),
  dc-megablast, one query batch: DUST step (`dust_s`) 0.211 -> 0.019 s, wall 2.19 -> 2.01 s (on/off 0.92, 18
  pairs in two windows); (b) 3,492 cDNAs (10 Mb), two big batches: 0.286 -> 0.040 s, 2.95 -> 2.71 s (0.92, 4
  pairs); (c) blastn with (a)'s queries (6 batches): 0.211 -> 0.022 s, 1.83 -> 1.67 s (0.91, 8 pairs). Max RSS
  within the pairs: (a) 89-93 MiB off, 90-95 on (1-2 MiB higher with the switch on in every window, within the
  spread between runs); (b) 104-108 / 103-106; (c) 137-148 / 126-131. The speed of
  one thread on that machine changes between runs by up to 1.6x (the DUST step of (a) takes 0.131 s or 0.211
  s with the switch off); pairs run in different speed states give other ratios, so compare `dust_s` and
  pairs of the same state. Every output equals the sha256 of the run without the switch (and of NCBI
  `-subject`). Tests: the unit tests of `x_query_dust_par.rs` (pool == serial at 1, 2, 3 and 7 threads, with
  empty sequences and batches of one; and a probe called once per query shows that the masks are made by
  at least two workers of a 3-thread pool, and by none with the switch off, one thread or one query) and
  `LOSAT/tests/blastn_query_dust_par.rs` (off / on / shadow equal at 1 and 4 threads for the three tasks, with
  `-dust 10 32 2 -lcase_masking`, `-dust no`, `LOSAT_X_QBATCHPAR` slots of two threads, and split query chunks).
- `LOSAT_X_DIAGLV=1` (not in the strict set: off unless set to any value, whatever `LOSAT_X_DEFAULT` says)
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
  values read) and in the strips folds of `LOSAT_X_PAIRPAR` and `LOSAT_X_PAIRPARDISC`; the off-diagonal
  retrieves (`scan_range`) stay walks. With the switch off the table has no bounds, an insert tests the
  length of an empty vector, and the retrieves walk as before (callgrind, release: +0.005 to +0.10 %
  instructions against the base).
  Shadow: `LOSAT_X_DIAGLVSHADOW=1` (implies the switch) walks every skipped bucket and panics unless each
  of its cells lies below the bound and at or below the bucket's bound (the bound invariant), and unless
  the cell the retrieve's walk returns gives the seed's decisions (the `last_hit` test, then with two hits
  the window test) of no cell (the bound's sufficiency); one `[X_SHADOW] DIAGLV retrieves= skipped=
  verified= decisions= cells= identical` line on stderr at the end of the search (`cells`: skipped
  retrieves whose walk found the diagonal's cell). `LOSAT_X_PAIRPARSHADOW` and `LOSAT_X_PAIRPARDISCSHADOW`
  compare the bounds of the two tables too.
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
- `LOSAT_X_PAIRPARFOLD=1` (not in the strict set: off unless set to any value, whatever `LOSAT_X_DEFAULT`
  says) blastn, where the strips paths of `LOSAT_X_PAIRPAR` (one-hit fold: megablast, blastn) and
  `LOSAT_X_PAIRPARDISC` (two-hit fold: dc-megablast) run, on a pool of 2+ threads: the diagonal-hash fold
  of the chunk is split by bucket (`bucket % P`, one partition per thread) over the owner and `threads - 1`
  jobs of an `in_place_scope` (`blastn/blast_engine/x_pair_seed_part.rs`). Each participant scans the
  next strip when it lies within `2 x threads` strips of the slowest partition, cuts its seeds by
  partition (a stable counting sort) and folds partition p of the next strip of any partition it can
  lock, into p's own table, with the folds of the strips paths (`fold`, `fold_two_hit`, a const
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
  it stays a figure of the owner's wall time. Shadow: `LOSAT_X_PAIRPARSHADOW=1` /
  `LOSAT_X_PAIRPARDISCSHADOW=1` with the switch run the partitioned stage on the copies; with the switch
  their lines end with `fold=part partitions= merges=` (without the switch the lines are as before).
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
  (0.577, 4 pairs), with `LOSAT_X_QBATCHPAR=1` 1.620 -> 0.959 s; blastn with the first query set 2.832 ->
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
- `LOSAT_X_QBATCHPAR=1` (BLASTN, all tasks, `-subject` with two or more threads; not in the strict set:
  off unless set to any value, whatever `LOSAT_X_DEFAULT` says) NCBI's query batches are searched side by
  side ahead of the batch loop of `run_in_pool` (`blastn/blast_engine/x_qbatch_par.rs`). The size of each
  batch after the first comes from the successful initial extensions of the batch before (NCBI's
  `CBatchSizeMixer`), so a planner advances a copy of the loop's mixer with predicted extensions (the hit
  rate per letter of the last batch of the same size, times the planned batch's letters; 0 for a batch that
  is split into query chunks) and plans the keys (first record, end record, batch size) of the next `2k`
  batches. `k` slots search them ahead, each with `m` threads (`k x m <= n`): `m = 1` on a plain thread,
  `m > 1` on a pool of `m` threads whose first thread is the slot's. The loop computes every batch's key
  from its own mixer as before and takes a slot's result only under that key, so it alone decides the
  batch boundaries, the errors, the warnings, the appends and the order of the output; the first batch is
  always its own; its thread waits while slots search, and it searches on the whole pool only when no
  slot runs. After every batch the plan is made again from the true mixer: planned batches the new plan
  no longer holds are dropped (queued), cancelled (running) or discarded (finished). Before the loop writes
  anything (the reports, or what an error writes) every slot is cancelled and joined, and once a slot's
  batch returns an error no batch planned after it starts. Cancellation: a slot's search carries a flag
  that `search_query_batch` reads at its start, between subjects and 5 Mb subject chunks, between query
  chunks and before its traceback; a cancelled search ends with an error the loop never takes (the loop's
  own searches have no flag). Memory: a slot starts a batch only while the estimates of the running
  batches and the sizes of the finished results not yet taken, with it, stay within the budget plus the
  loop's next batch (always admitted). Waste: two mispredicted next keys within eight batches, or more
  discarded batches (finished results dropped and searches cancelled) than `2k` plus a quarter of the
  batches taken, stop the speculation for the search. A batch's `LOSAT_TIMING` block is printed when the
  loop takes the batch (discarded batches print none). Off (the loop as before) with `-db`, one thread,
  `BLEMIR_DEBUG`, `LOSAT_DEBUG_BLASTN`, `LOSAT_DEBUG_COORDS[_START]`, the blastn trace or
  `LOSAT_WASI_THREADS_DEBUG`, subjects of at least as many 5 Mb work units as threads, fewer than two
  batches after the first, or a budget of fewer than two slots; not in WebAssembly or without the
  `parallel` feature. Slots: one batch in flight is estimated at `4^w x 4` bytes of lookup table plus 21
  bytes per query letter plus 220 bytes per predicted extension (`w` the batch's lookup width; calibrated
  to 50 MiB per 1 Mb blastn batch and 116 MiB per 5 Mb dc-megablast batch), the budget of the batches
  beyond the loop's one is `LOSAT_X_QBATCHPARMEM` MiB, by default 96 MiB (blastn, blastn-short) or 28 MiB
  (megablast, dc-megablast) per thread beyond the first (0.8 of NCBI BLAST+ 2.17.0's max-RSS growth per
  thread on the audit input; provisional; dc-megablast's 5 Mb batches then need six threads or more),
  `k = min(1 + budget / estimate, n, units of work left)` lowered to the fewest slots that need as many
  rounds, and `m = n / k`; `LOSAT_X_QBATCHPARSLOTS=k` sets `k` (tuning). Without `LOSAT_X_STATS` the switch
  writes nothing of its own; with it, one line on stderr at the end of the search: `[X_STATS] QBATCHPAR
  batches= precomputed= loop= mispredicted= discarded= cancelled= panicked= waste_bound= k= m= idle=
  budget_mib= estimate_mib= peak_mib= off=` (`precomputed`: batches taken from a slot; `discarded`:
  finished or cancelled batches dropped (a search whose flag was set when it returned counts as cancelled,
  whatever it returned); `panicked`: slot batches that panicked, kept or dropped; `waste_bound`: the
  discards that stop the speculation, at the end; `idle`:
  threads no slot uses; `peak_mib`: the most memory admitted at once; `off`: why it did not, or stopped
  to, speculate). Shadow: `LOSAT_X_QBATCHPARSHADOW=1` (implies the switch) searches every batch the loop
  takes from a slot on the loop as well (on the whole pool while the slots run, so more than `n` threads
  are active) and compares every field (each query's hit list with every HSP field and edit script
  through its `Debug` text, the extensions, the masks, the Karlin-Altschul blocks, the search spaces, or
  the two errors); any slot panic, taken or not, is a difference too; the first difference prints the
  summary line and panics; at the end `[X_SHADOW] QBATCHPAR batches= compared= diffs= panicked=
  mispredicted= wasted=`. `LOSAT_X_QBATCHPARMISS=N` (test
  only) corrupts every N-th prediction. A slot's panic is kept and resumed only if the loop takes its
  batch; the release profile aborts on a panic, so there a panic in a batch a slot searches ends the
  process (no slot runs while output is written, which the writers check with debug assertions). Past an
  error at batch j, the slots may have searched up to the plan's window (2k - 1 batches) before it is
  known; once it is known, at most k - 1 still run and are cancelled when the loop meets j. With debug assertions, a check after the slots are
  joined fails if a slot started rayon's global pool. The counters and `LOSAT_TIMING` end-of-run lines of
  other switches (PAIRPAR, PAIRPARDISC, GREEDYPRELIM, the megablast lookup fill) and
  `LOSAT_X_DPCAPTURE` files count the work of the slots' threads too, discarded batches included, so
  they differ from a run without the switch. Tests: the unit tests of `x_qbatch_par.rs` and
  `LOSAT/tests/blastn_qbatchpar.rs`.
- `LOSAT_X_TBNPAR=1` tblastn: preliminary stage per (frame, chunk) and redo per match in
  parallel, one pool for the search. Shadow: `LOSAT_X_TBNPARSHADOW=1`.
- `LOSAT_X_TBNBATCH=1` (needs `TBNSSIDE` and its kept translations, so four or more query
  batches, and at least one batch per thread) tblastn query batches searched side by side,
  appended in input order.
- `LOSAT_X_TBNOIDCHUNK=1` (in the strict set since 2026-10-13) tblastn with 2+ threads and several
  subjects: one entry into the search pool takes four chunks of subjects per thread instead of one
  subject per thread, and each chunk is one task that the pool gives to whichever thread is free;
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
  `-db` 6.6-7.4 s; on the workstation 16.4 s against 62.1 s, NCBI 16.1 s). Nothing changes at one thread, nor where `TBNBATCH` searches the batches
  side by side (each batch on one thread). With `LOSAT_WASI_THREADS_DEBUG=1` the
  `[losat-thread-activity]` line of the subject stage ends with `pool_entries=` (the entries of
  one query batch); `LOSAT/tests/tblastn_oid_chunk.rs` reads it.
- `LOSAT_X_TBNOIDNEXT=1` (in the strict set since 2026-10-15) tblastn with 2+ threads and several
  subjects, where `TBNOIDCHUNK` makes the chunks: no wait between entries. The threads of the pool
  take the next chunk when they are free, as NCBI's search threads ask the database for the next
  chunk (`seqdbimpl.cpp:333-344`), instead of waiting at the end of each entry for its last chunk.
  The thread that owns the search takes the lists back chunk by chunk in OID order and enters them
  into the collector as before, and searches a chunk itself while the next one in order is not
  finished (`tblastn/x_oid_flow.rs`). At most `4 * threads` chunks are taken or finished beyond the
  oldest one not taken back; a thread waits before a new chunk when that window is full, so the
  memory of the waiting lists stays bounded. The chunks and the search of each subject are those of
  `TBNOIDCHUNK`, so the report is the same. `=0` turns it off; with `LOSAT_X_TBNOIDCHUNK=0` (or
  `LOSAT_X_DEFAULT=0` alone) there are no chunks and it changes nothing. With
  `LOSAT_WASI_THREADS_DEBUG=1` the search is one entry (`pool_entries=1`) and the line adds
  `flow_chunks=` and `flow_peak_waiting=` (the most finished chunks that waited at once);
  `LOSAT/tests/tblastn_oid_chunk.rs` reads them.
- `LOSAT_X_TBNOIDUNITS=1` (in the strict set since 2026-10-15) tblastn with 2+ threads and several
  subjects, where `TBNOIDCHUNK` makes the chunks (with or without `TBNOIDNEXT`): a subject longer
  than a chunk's budget (a chunk of its own, such as a chromosome) runs its (frame, chunk) units on
  the pool, as the single subject of a search does with `TBNPAR`, instead of one after another on
  the thread that took it. Threads that have no chunk to search (at the end of an entry, or when no
  chunk is left) take its units; its merge in NCBI order, linking, E-values and reaping run on its
  thread after its last unit, so the report is the same. NCBI's thread searches a subject whole, so
  at 16 and 32 threads the longest chromosome bounds both programs; this goes beyond NCBI. Other
  subjects keep their units on their thread. With `TBNOIDNEXT` the threads of the flow start
  together (each waits until all have entered the flow) so that a thread waiting for a unit can only
  take units, never a loop of the flow. `=0` turns it off. With `LOSAT_WASI_THREADS_DEBUG=1` each
  subject whose units ran on the pool prints a `prelim_units` stage line with
  `parallel_selected=true`; `LOSAT/tests/tblastn_oid_chunk.rs` counts them.
- `LOSAT_X_TBXNEXT=1` (TBLASTX, speed; in the strict set since 2026-10-17, on by default; acts only
  where the parallel subject loop runs: two or more threads and more than one subject, native builds) each thread of the pool takes the next subject from one
  counter when it is free (`tblastx/x_subject_next.rs`, one broadcast job per thread), with one worker
  state per thread made when it takes its first subject. Without it the loop gives rayon the subject indices
  as one indexed range (`for_each_init`), which rayon cuts into pieces before the search, so the long
  subjects a genome lists first (the chromosome arms) fall into one piece and one thread searches them in
  order. `process_subject` is unchanged and `report::final_hit_order` orders the HSP lists by query and
  subject index, so the reports are the same. Use it with `LOSAT_X_TBXFRAMES`: alone it holds the frames of
  several long subjects at once (3 fruit fly transcripts x the fruit fly genome database at 8 threads: 892
  MiB, NCBI 359 MiB). With both at 8 threads (node 0, CPUs 8-15, medians of 3 alternating pairs, reports
  identical): that search 11.48 -> 2.58 s (NCBI `-db` 2.73 s), max RSS 257 -> 316 MiB; with the genome FASTA
  as `-subject` 12.83 -> 3.14 s (410 MiB); 10 transcripts x the zebrafish cDNA database 8.51 -> 8.01 s.
- `LOSAT_X_TBNSTATE=1` (in the strict set since 2026-10-15) tblastn, the unit path (`TBNPAR`): the
  scan state of a subject's (frame, chunk) units (the diagonal table of the two-hit word finder, the
  offset pairs, the initial hits, the gapped alignment memory) is kept by the thread from one
  subject to the next, as NCBI's search thread keeps its word-finder and gapped alignment structures
  for every subject of its loop (`blast_engine.c:1002-1003,1412-1413`), and each unit ends with
  NCBI's step of the diagonal offset (`Blast_ExtendWordExit`, `blast_extend.c:166-171`: the offset
  moves on by the unit's length plus the window; the table is cleared only when the offset reaches
  INT4_MAX / 4) instead of the next unit clearing the whole table. After the step every entry lies
  more than a window behind any hit of the next unit, so its first hit on a diagonal only records
  it, as on a cleared table: the same initial hits by construction. Without the switch each subject
  makes its own state and each unit clears a table of next_pow2(query residues + window) entries of
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
  it over. `=0` turns it off; it acts at every thread count. With `LOSAT_X_DEFAULT=0` the reference
  loop runs instead of the unit path and the variable changes nothing. `LOSAT_X_TBNPARSHADOW=1`
  compares the unit path of every subject with the serial reference. Tests: `tblastn/search_init.rs`
  (the table kept across units, across the clear and taken over by another query set), the last test
  of `LOSAT/tests/tblastn_oid_chunk.rs`.
- `LOSAT_X_TBNXLATE=1` (in the strict set since 2026-10-15) tblastn: the six translated frames of a
  subject are read from tables, as NCBI makes them from the compressed subject
  (`blast_util.c:1071-1090`: one table of 64 residues per strand, `s_BlastGetTranslationTable`, and
  one table read per codon, `BLAST_TranslateCompressedSequence`). The subject's letters become
  ncbi2na in one pass; a gap or an ambiguous letter draws from the generator seeded with the
  subject's length, in subject order, as `s_Ncbi4naToNcbi2na` does
  (`blast_objmgr_tools.cpp:427-475`). The reference path reaches the same residues through letters
  (the subject as ncbi4na, as ncbi2na, as `ACGT`, its reverse complement, then `GeneticCode::get`
  per codon). The tables are filled from the reference path's own functions, so every residue is the
  reference's by construction. `=0` turns it off. It acts where the subject plans of the unit path
  are made (`TBNPAR`, on by default); with `LOSAT_X_DEFAULT=0` the reference loop translates each
  subject itself and the variable changes nothing. Shadow: `LOSAT_X_TBNXLATESHADOW=1` makes both
  translations for every plan and stops (panic) if a frame, a unit or an error text differs; one
  `[X_SHADOW] TBNXLATE` line on stderr at the end gives the plans and residues compared, or says
  that nothing was compared. With `LOSAT_X_STATS=1` a `[X_STATS] TBNXLATE table_plans=` line counts
  the plans made from the tables; `LOSAT/tests/tblastn_subject_translate.rs` reads both.
- `LOSAT_X_TBNSCANBF=1` tblastn (speed; not in the strict set: off unless set, to any value;
  `tblastn/x_scan_bf.rs`): the lookup scan of the unit path (`TBNPAR`, `search_seed.rs`
  `x_scan_unit`) takes the subject in blocks of 256 words. Each word's position and table index are
  written to a list whose length advances by the word's presence-vector bit, without a branch; the
  hits of the listed words are then copied in position order. When the offset array has no room for
  a word's hits, the scan stops at that word and the next call resumes there, as NCBI's
  `s_BlastAaScanSubject` (`blast_aascan.c:75-127`) does, so every call returns NCBI's pairs and scan
  range (test `x_scan_bf_equals_ncbi_scan_on_random_subjects` against
  `s_blast_aa_scan_subject_one_range` on random subjects, with offset arrays down to the longest
  chain). A backbone cell's hits are copied as three pairs, so the offset array has three entries
  more. With `LOSAT_X_DEFAULT=0` the reference loop scans and the variable changes nothing; tblastx,
  blastp and blastx keep their scans. perf stat on the workstation (one thread, a 411-residue
  protein x the human genome `-db`, the same output): branch misses 1.96 G -> 0.46 G, cycles
  136.0 G -> 96.7 G, instructions 285 G -> 348 G. Timed on the workstation (perf-si-20261012
  session 20, node 0, NCBI BLAST+ 2.17.0 `-db`, LOSAT without and with the switch in turn, one
  hash per case; time with the switch / without): human genome x one protein 0.73 at one thread
  (45.6 -> 33.2 s, NCBI 44.9 s) and 0.77 at 8 threads (6.79 -> 5.21 s, NCBI 7.35 s); zebrafish
  genome x one protein at one thread 0.79 (31.3 -> 24.6 s, NCBI 31.3 s). With many query proteins
  the bit is set at most words and the branch is predicted: 0.98-0.99 (fly genome x one batch at
  1 / 8 threads, E. coli MG1655 x 16 batches at 8, 16S x one batch at 8). For tblastn's table (word
  size 3, 5-bit letters) the candidate loop has its shifts as constants (`x_candidates::<3, 5>`;
  other tables take the general loop): a further 0.97 of the time on the one-protein cases
  (cycles 96.9 G -> 93.4 G), the same CPU time on the many-protein ones. An AVX2 form of the loop
  (`vpgatherdd` of the presence words) cut instructions by a third but was 2-7 % slower (IPC 3.6 ->
  2.3) and is not used.
- `LOSAT_X_TBNSTORE=1` tblastn `-subject` without `-subject_loc`: each subject record's letters are
  packed as soon as the reader has read it (`sequence/packed_subject.rs`, the shape of a BLAST
  database volume: ncbi2na 4 bases per byte with the remainder count in the last byte, the runs of
  letters other than `ACGT` as ncbi4na codes, the lower-case intervals), so the letters of the
  whole set are never held; the records keep their deflines. With two or more threads a second
  thread packs the records while the reader reads the next ones (`x_packed_while_read`: batches of
  up to 1,024 records or 4 MiB of letters, handed over one at a time), so the packing does not add
  to the read; the records and their stores are those of the one-thread read, in the same order. The unit path of the preliminary
  stage (`TBNPAR`, on by default) translates a subject's six frames straight from its packed bases
  (`x_packed_plan` in `tblastn/x_subject_translate.rs`: the codon tables of `TBNXLATE`, the 2-bit
  codes unpacked 49,152 bases at a time, each frame filled at its residues' places), as NCBI's
  `BLAST_GetAllTranslations` translates the compressed subject, so a thread holds the frames of
  the subject in hand and no letters, resolved bases or reverse complement; this is the
  translation whatever `TBNXLATE` says, and `LOSAT_X_TBNXLATESHADOW=1` also makes the reference
  plan from the resolved bases and stops if they differ. The reference loops (`LOSAT_X_DEFAULT=0`)
  read the packed bases as letters (`ACGT`, each ambiguous base as resolved), which resolve to
  themselves; with `-lcase_masking` the stage reads the letters. The redo (composition mode 2),
  its translations and the pairwise rows (outfmt 0) read windows of the subject (`SubjectNt` in
  `packed_subject.rs`: the HSP's codons, or the window around an HSP that NCBI translates), made
  in upper case from the store when read (every reader ignores the case); only the traceback of
  composition mode 0 makes a subject's whole letters again, exactly as the reader gave them (the
  comparison tests in `packed_subject.rs`). The packed bases are those NCBI's subject set-up makes
  (`encode_subject_ncbi2na_packed`, `s_Ncbi4naToNcbi2na` in `blast_objmgr_tools.cpp:427-475`; one
  pass over the letters). Memory: about a quarter of a byte per base held instead of one; per
  thread the frames of the subject in hand (two bytes per base) instead of four bytes per base
  (letters or resolved bases, reverse complement, frames). A subject range keeps the letters (the
  records are cut before the search), and so does a set with a byte that has no ncbi4na code
  (only an API caller can make one); the library entry `run_local` reads the caller's records.
  In the set since 2026-10-16 (on unless `0`, or `LOSAT_X_DEFAULT=0` and not set). With
  `LOSAT_X_STATS=1` a `[X_STATS] TBNSTORE
  subjects= stored_bytes= letters_made= resolved_made= packed_plans=` line gives the subjects held
  packed, the bytes of their stores, the times their whole letters and their resolved bases were
  made, and the plans translated from the packed bases; `LOSAT/tests/tblastn_subject_store.rs`
  reads it.
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
  the call to the reference function. In the strict set since 2026-10-14 (`=0` turns it off). Shadow:
  `LOSAT_X_BPGAPSHADOW=1` runs both for every restricted alignment and stops (panic) if
  `(a_offset, b_offset, best_score)` differ; one `[X_SHADOW] BPGAP compared= not_taken=` line on
  stderr at the end. With `LOSAT_X_STATS=1` a `[X_STATS] BPGAP fast= not_taken=` line.
- `LOSAT_X_PAIRPAR=1` TBLASTX seed stage (scan, two-hit test, ungapped extension) of each subject
  as units of (subject frame, chunk, query context), 36 per one-query pair, on the search pool (in
  order on the calling thread with one thread); each unit scans with a lookup table holding only its
  context's chains (order kept) on a fresh unmasked diagonal table, and the unit lists go through the
  reference's sort and merge in the reference order (`tblastx/x_pair_par.rs`, exactness argument in
  its header). Only when every query context, or every subject chunk, is at most `slack` = diagonal
  table size minus query length long, and `slack` plus twice the longest subject chunk stays under
  2^29 (it does unless `-window_size` is near 2^29); otherwise, and with the trace, debug,
  diagnostics and chunk switches, the reference loop runs (`LOSAT_X_STATS=1` counts `fallback_guard` / `fallback_gate`;
  `LOSAT_TIMING` prints one `x_pairpar` line). Strict-set candidate (S-E, round 5).
- `LOSAT_X_PAIRPARSHADOW=1` shadow of `LOSAT_X_PAIRPAR`: the unit path and the reference loop for
  every subject, the subject's combined ungapped HSP list compared field by field (the first
  difference stops the run with subject, frame, index and both values); one `[X_SHADOW] PAIRPAR`
  summary line at exit (subjects, units, HSPs compared, fallbacks).
- `LOSAT_X_PAIRLUT=1` (in the strict set; acts while `LOSAT_X_PAIRPAR` is on, not with its
  shadow) TBLASTX:
  the lookup tables of the unit path, one per context group, built directly from the query instead
  of being filtered out of the batch's full table, which is then not built: for each group the
  exact-match index of the group's lookup locations, then the count, layout and fill passes of
  `LOSAT_X_BXLUT` (`tblastx/lookup/x_pair_lut.rs`, exactness argument in its header), the groups
  in order on the calling thread (tuning only: `LOSAT_X_PAIRLUT_POOL=1` builds them on the search
  pool, a few ms faster per batch and more memory). The reference code holds a stub (the full
  table's header fields and `longest_chain`, no entry). Whether the unit path runs is decided before the lookup
  build by the function `prepare` decides with; the run stops if `prepare` disagrees or if the
  reference frame loop would scan the stub. The reference builder and `build_groups` run as before
  when the unit path does not run (gate, guard), in diagnostics mode, and when a lookup location
  does not lie inside one group. `LOSAT_X_STATS=1` prints
  `[X_STATS] PAIRLUT batches_direct= groups= fallback=`; `LOSAT_TIMING` prints one `x_pairlut` line
  (the direct build; it is part of `build_lookup`, and `sublookup` of the `x_pairpar` line is 0).
- `LOSAT_X_PAIRLUTSHADOW=1` shadow of `LOSAT_X_PAIRLUT` (while `LOSAT_X_PAIRPAR` is on, or with its
  shadow): the
  full table and `build_groups` as before, and the direct groups as well; every group table (cells,
  overflow array, presence vector, header, `longest_chain`), `q_min` / `q_max`, the stub's header
  and `longest_chain`, and the contexts are compared, and the first difference stops the run; the
  search uses the reference tables. One `[X_SHADOW] PAIRLUT batches= groups_compared= fallback=`
  line at exit.
- `LOSAT_X_BXPAR=1` blastx query-chunk, redo and per-context SEG parallelism;
  `LOSAT_X_BXPOOL=1` one pool for the search; `LOSAT_X_BXBATCH=1` (needs `BXPOOL`) small query
  batches searched together, output in input order; `LOSAT_X_BXCHUNK=64` subjects per worker;
  `LOSAT_X_BXWAVE` tuning only.
- `LOSAT_X_BXFLOW=1` (BLASTX, speed; not in the strict set: off unless set, to any value; acts with
  `BXBATCH` and `BXPOOL` at two or more threads): the queued query batches are searched as a flow
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
  threads 0.76; 5 cDNAs (one batch) 1.00.

Counters and tests: `LOSAT_X_STATS=1` (counters on stderr; DP and Newton counts need the
`xstats` cargo feature), `LOSAT_X_ADJMEMO_STATS=1` (distinct composition-adjustment inputs),
`LOSAT_FUZZ_CASES=<n>` (random groups of the linking fuzz test; random problems of the X-drop DP
kernel tests).

X-drop DP capture and replay (diagnostics, not part of the strict set):

- `LOSAT_X_DPCAPTURE=<dir>` records every Nth call of each thread that reaches the vector X-drop
  DP (`LOSAT_X_DPFAST`) as a self-contained problem (residues, matrix, gap costs, x_drop, result
  and edit script) in `<dir>/xdp-<pid>-t<k>.xdp`, one file per thread, appended
  (`utils/x_xdrop_capture.rs`, format in its header). `LOSAT_X_DPCAPTURE_EVERY=<n>` sets N
  (default 25). The search output does not change.
- `LOSAT_X_DPREPLAY=<dir>` input of the ignored test `replay_captured_problems`
  (`cargo test --release replay_captured_problems -- --ignored --nocapture`): every captured
  problem through the scalar loops, the 128-bit, the `LOSAT_X_DPW128` and the 256-bit kernels and
  the dispatch, compared with the recorded result, then the three kernels timed (median ns per cell
  and per call, traceback and score-only, 8- and 16-bit lanes). `LOSAT_X_DPREPLAY_REPS` (default 7),
  `LOSAT_X_DPREPLAY_OUT=<file>` (append the table to a file). The ignored test
  `w128_replay_captured_problems` does the same for the 128-bit and the `LOSAT_X_DPW128` kernels
  only and also runs in the wasm32-wasip1 test build (the libtest `.wasm` under a WASI host such
  as Node's `wasi` module with a root preopen), which times the simd128 kernels.
