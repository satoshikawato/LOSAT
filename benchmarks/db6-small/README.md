# Database search benchmark, small public databases (DB-6)

This directory holds the tables, scripts and raw row files of a compute-program
benchmark of LOSAT `-db` searches against NCBI BLAST+ 2.17.0 `-db` searches. It
characterizes program output and wall-clock execution; it is not a biological
interpretation of the input sequences. The summary is in the repository
[README](../../README.md#database-searches).

## What was measured and why

The question: on small public databases, how long does `LOSAT -db` (the default
configuration) take compared with NCBI BLAST+ 2.17.0 `-db` at 1, 8 and 16 threads, with
what peak memory, and is every LOSAT output the same bytes as NCBI BLAST+
`-db -num_threads 1`? It also times database construction (`makeblastdb` and
`LOSAT makedb`). Seven search cases cover megablast, blastn (`-task blastn`),
blastp and blastx.

Not covered: large databases (nr, nt), tblastn and tblastx database searches, and
LOSAT `-subject` searches; see [Not measured](#not-measured).

## Data

No FASTA file is included in this repository. Each input is a public data set that
can be fetched again. The table gives sequence and letter counts and the first 16
hex digits of the SHA-256 of the file used.

| Name | Source | Sequences | Letters | SHA-256 (first 16 hex) |
| --- | --- | --- | --- | --- |
| 16S | RefSeq 16S rRNA sequences of bacteria and archaea, collected 2026-10-12 | 28,036 | 40.6 Mb | 7c5eec160a6291c5 |
| Swiss-Prot | UniProtKB/Swiss-Prot release 2026_03 | 575,748 | 209.0 M aa | 1372f18e218d7610 |
| yeast proteome | UniProt proteome UP000002311 (taxon 559292, *Saccharomyces cerevisiae*), FASTA | 6,066 | 2.9 M aa | 8ca7d8ccfd337f03 |

Query sets:

- `q100` (megablast): 100 full-length 16S sequences (143,418 nt) taken from the 16S set.
- `amp100` (megablast and blastn task): letters 341-800 (1-based, 460 nt) of each
  `q100` record, 100 sequences (46,000 nt). They are an amplicon-like (V3-V4)
  stand-in made by slicing, not by primer matching.
- `sp_q1` (blastp): 1 Swiss-Prot protein of 379 aa. `sp_q100` (blastp): 100 Swiss-Prot
  proteins, 35,048 aa.
- the yeast proteome itself (blastp, all against all).
- `t10` (blastx): the first 10 records (10,954 nt) of a 100-transcript set used in
  an earlier LOSAT blastx benchmark; the 10 transcripts are not an extra public
  download.

The databases are built in the run directory from these files (timings below). The
searches of both tools use the `makeblastdb` default database (format version 5, as
a user gets it). One extra untimed LOSAT run per case at 8 threads searches the
`LOSAT makedb` database and must give the same bytes.

## Machine and build

- Workstation with two Intel Xeon Gold 6338 CPUs, governor `powersave`. Each timed
  run used 16 physical cores of one NUMA node (`numactl --physcpubind`) with memory
  bound to that node (`--membind`). One thread ran on the first picked core, eight
  on the first eight. Before each timed run the driver picked cores whose own and
  hyper-thread-sibling use was below 10 % and recorded how many other CPUs of the
  node were busy (above 50 %) at the start and end of every run; runs where that
  rose by 3 or more were kept as disturbed and repeated (14 runs in the small-tier
  run, whose times were within 1 % of the kept runs; `rows/db6_runs.tsv` keeps
  them with the busy counts). Other work ran on other CPUs of the machine during
  the runs.
- Small-tier run (2026-10-17): LOSAT `origin/dev` 0b7335ce (release build, SHA-256
  of the binary starts with 68fb2529df9ee280), run on node 0 (CPUs 16-31 with their
  hyper-thread siblings); default configuration with no `LOSAT_*` switch set
  except `LOSAT_BLASTX_WORKER_LOG`, which only `--all-features` builds need. Rust 1.92.0.
- NCBI BLAST+ 2.17.0 (`blastn`, `blastp`, `blastx`, `makeblastdb`).
- Later builds named in the section on the two former losses: 6d0c7b1e (PR #287's
  head, on dev 7a5a7f87) and 9560b7c1 (PR #289's head, the default set that
  includes `LOSAT_X_BPDBCUT`).

## Protocol

For each case and thread count (1, 8, 16): one untimed warmup of each tool, then
three timed rounds; the tool order alternates by round (NCBI first, LOSAT first,
NCBI first). The tables give the median wall time and the min-max range of the
three. Wall time is `time.perf_counter()` (monotonic) around the search process;
user and system CPU and maximum resident set size (`ru_maxrss`) come from `wait4`
(`scripts/s7_time.py`). Output is `-outfmt 6` written to a file, with the default
options of each program. Database builds use the same warmup and three timed runs
per tool and database.

## Exactness

In every case, every LOSAT output (1, 8 and 16 threads, warmups and timed runs, and
the run on the `LOSAT makedb` database) is NCBI BLAST+ 2.17.0 `-db -num_threads 1`'s
output byte for byte (SHA-256 compared). NCBI's own 8- and 16-thread outputs equal
its 1-thread output in all seven cases. Row counts and the first 16 hex digits of
the SHA-256 are in `rows/db6_summary.txt`.

## Search time

Small-tier run, LOSAT build 0b7335ce, median of 3 in seconds; the range is the
min-max of the three. Ratio is LOSAT / NCBI of the medians.

| Case | Threads | NCBI s (range) | LOSAT s (range) | Ratio | Max RSS MiB NCBI / LOSAT |
| --- | --- | --- | --- | --- | --- |
| MB-q100: megablast, 100 full-length 16S x 16S | 1 | 111.19 (110.71-111.77) | 65.90 (65.67-65.94) | 0.59 | 140 / 52 |
| | 8 | 17.22 (17.21-17.22) | 13.08 (13.01-13.09) | 0.76 | 147 / 146 |
| | 16 | 10.62 (10.61-10.90) | 6.59 (6.55-6.69) | 0.62 | 140 / 164 |
| MB-amp: megablast, 100 slices x 16S | 1 | 34.70 (34.66-34.74) | 21.64 (21.61-21.72) | 0.62 | 196 / 46 |
| | 8 | 6.74 (6.71-6.80) | 4.57 (4.48-4.58) | 0.68 | 206 / 118 |
| | 16 | 5.46 (5.28-5.65) | 2.32 (2.32-2.35) | 0.43 | 208 / 141 |
| BN-amp100: blastn task, 100 slices x 16S | 1 | 128.90 (128.60-128.96) | 54.69 (54.09-55.32) | 0.42 | 149 / 32 |
| | 8 | 18.73 (18.68-18.73) | 10.75 (10.68-11.32) | 0.57 | 158 / 105 |
| | 16 | 10.78 (10.78-10.88) | 5.46 (5.44-5.52) | 0.51 | 149 / 125 |
| BP-q1: blastp, 1 protein x Swiss-Prot | 1 | 3.70 (3.70-3.71) | 3.67 (3.66-3.68) | 0.99 | 273 / 13 |
| | 8 | 0.66 (0.66-0.69) | 0.53 (0.52-0.58) | 0.80 | 269 / 28 |
| | 16 | 0.51 (0.51-0.52) | 0.36 (0.36-0.36) | 0.71 | 263 / 44 |
| BP-q100: blastp, 100 proteins x Swiss-Prot | 1 | 249.70 (249.66-250.78) | 153.47 (153.33-153.75) | 0.61 | 405 / 119 |
| | 8 | 36.35 (36.27-36.35) | 20.53 (20.44-20.59) | 0.56 | 441 / 176 |
| | 16 | 20.23 (20.19-20.30) | 10.96 (10.75-11.00) | 0.54 | 440 / 183 |
| BP-yeast: blastp, yeast proteome all against all | 1 | 885.72 (885.57-886.05) | 342.91 (342.57-343.81) | 0.39 | 80 / 32 |
| | 8 | 114.34 (114.31-114.37) | 103.79 (102.91-103.94) | 0.91 | 172 / 54 |
| | 16 | 61.36 (59.93-61.39) | 97.91 (97.88-99.06) | **1.60** | 212 / 62 |
| BX-t10: blastx, 10 transcripts x Swiss-Prot | 1 | 30.19 (30.18-30.20) | 28.25 (28.19-28.27) | 0.94 | 294 / 219 |
| | 8 | 4.12 (4.12-4.12) | 7.11 (7.09-7.14) | **1.73** | 290 / 221 |
| | 16 | 2.23 (2.23-2.23) | 4.91 (4.88-4.93) | **2.20** | 282 / 222 |

The range of 32 of the 42 sets of three runs is at most 2 % of the median. The
widest are BP-q1 LOSAT at 8 threads (0.52-0.58 s, 11 %), MB-amp NCBI at 16 threads
(5.28-5.65 s, 7 %) and BN-amp100 LOSAT at 8 threads (10.68-11.32 s, 6 %); none of
them changes which side of 1 a ratio is on. LOSAT's maximum RSS is at most NCBI's
in every case and thread count except MB-q100 at 16 threads (164 vs 140 MiB).

The three ratios in bold are the two losses of this run (BP-yeast at 16 threads;
BX-t10 at 8 and 16 threads). Their causes and the numbers after the fixes are in
the next section. NCBI blastx at 16 threads took 2.23 s here and 2.5-3.2 s in the
follow-up window of the same day.

Per-run values (wall, user and system CPU, RSS, row count, output hash, load,
busy-CPU counts) are in `rows/db6_runs.tsv`; the per-case summary is
`rows/db6_summary.txt`.

## Database builds

Small-tier run, median of 3 (min-max in `rows/db6_summary.txt`), one thread; neither
tool has a thread option for building.

| Database | `makeblastdb` (default, v5) | `makeblastdb -blastdb_version 4` | `LOSAT makedb` |
| --- | --- | --- | --- |
| RefSeq 16S | 0.69 s | 0.70 s | 0.17 s |
| Swiss-Prot | 8.47 s | 8.26 s | 2.02 s |
| yeast proteome | 0.16 s | 0.15 s | 0.03 s |

The sequence, header and index files of `LOSAT makedb` are byte for byte those of
`makeblastdb -blastdb_version 4` (0 differing bytes in `.nsq`/`.nhr`/`.nin` and
`.psq`/`.phr`/`.pin`). The sequence-file hashes of the three builds of each
database agree (`sq_sha16` in `rows/db6_summary.txt`). The timed searches use the
`makeblastdb` default (version 5) databases.

## Two former losses and their fixes

The small-tier run found two cases where LOSAT was slower than NCBI. Both were
traced and fixed afterwards; the README table shows the numbers after the fixes
(marked there). Timings are medians of the number of rounds stated.

### BP-yeast: blastp, yeast proteome all against all, 16 threads

Cause: the 2.9 MB database made 3 blocks of 1 MiB, so at most 3 threads had work
(LOSAT used 361 s of CPU in 99 s of wall time). Fix: `LOSAT_X_BPDBCUT` (#231; blocks
cut by database size and thread count), in the default set since #289 (merged into
dev 2026-10-18).

Evidence:

1. Same-window check at 16 threads (2026-10-17, build 0b7335ce with
   `LOSAT_X_BPDBCUT=1`, 3 alternating rounds with NCBI): 32.7 / 32.9 / 35.1 s, against
   NCBI 60.5-61.6 s and 99.2-99.3 s for the default of that build. Rows:
   `rows/db6_followup_runs.tsv`.
2. Post-fix numbers: paired run on the same workstation, 2026-10-18, build 9560b7c1
   (#289's head, the new default; node 1, 16 cores), with and without
   `LOSAT_X_BPDBCUT=0` in the same rounds (rows: `rows/bp_yeast_bpdbcut_pairs.tsv`;
   outputs identical):

| Threads | Default (build 9560b7c1) | `LOSAT_X_BPDBCUT=0` | On / off | Rounds | Max RSS MiB (median) |
| --- | --- | --- | --- | --- | --- |
| 16 | 22.889 s | 71.442 s | 0.320 | 3 | 52 |
| 8 | 38.623 s | 77.455 s | 0.499 | 2 | 45 |
| 1 | 279.109 s | 278.272 s | 1.003 | 1 | 32 |

That run had no NCBI arm. The README ratios (0.32, 0.34 and 0.37 at 1, 8 and 16
threads) divide these LOSAT medians by NCBI's small-tier medians of 2026-10-17
(885.7, 114.3 and 61.4 s); they are cross-run ratios. NCBI's maximum RSS in the
small-tier run was 80, 172 and 212 MiB at 1, 8 and 16 threads.

### BX-t10: blastx, 10 transcripts x Swiss-Prot, 8 and 16 threads

Cause: one query batch; the subjects ran in waves of threads x 64 that each waited
for their slowest chunk, so about 7 of 16 threads were busy at 16 threads
(7.1 vs 4.1 s at 8 threads and 4.9 vs 2.2 s at 16 in the small-tier run). The
opt-in `LOSAT_X_BXFLOW` changed nothing here (7.14-7.28 vs 7.16-7.68 s at 8 threads).
Fix: `LOSAT_X_BXSUBJFLOW` (#287, default ON, merged into dev 2026-10-18): the subjects
of the batch run as a flow of 64-subject tasks.

Post-fix numbers: paired run with NCBI in the same rounds, 2026-10-17 23:1x, build
6d0c7b1e (#287's head, on dev 7a5a7f87), node 1, 16 cores, on a machine with a
1-minute load of about 28 (other work ran on other CPUs). Medians of 3; outputs
identical to NCBI's. Rows: `rows/bx_t10_flow_pairs.tsv` (arms: `def` = default with the
flow, `off` = `LOSAT_X_BXSUBJFLOW=0`, `ncbi`; threads 2, 4, 8, 16).

| Threads | LOSAT default (flow) | LOSAT, flow off | NCBI | LOSAT / NCBI |
| --- | --- | --- | --- | --- |
| 8 | 3.899 s | 7.216 s | 4.116 s | 0.95 |
| 16 | 2.237 s | 5.002 s | 2.271 s | 0.99 |

At 1 thread the README uses the small-tier value (LOSAT 28.3 s vs NCBI 30.2 s,
ratio 0.94); the flow does not act at one thread. The 2- and 4-thread rows of the
same run (LOSAT 14.7 and 7.39 s vs NCBI 15.9 and 7.99 s) are in the row file.

## Files

The scripts are copied as they were run on the workstation; their paths, state
directories and markers are that machine's, so they are a record of the method
rather than something to run unchanged.

`scripts/`
- `ws_db6_bench.sh`: the DB-6 driver: builds the databases, runs the warmup and
  the three timed rounds of each case at 1, 8 and 16 threads, records every run,
  and hashes every output.
- `ws_db6_build.sh`: builds LOSAT at a given commit on idle cores (untimed).
- `ws_db6_queue.sh`: picks 16 idle physical cores of one node and starts the bench.
- `ws_db6_followup.sh`: the follow-up run of the two losses (`rows/db6_followup_runs.tsv`).
- `ws_s16_pairs.sh`: the paired run of 2026-10-17 for `LOSAT_X_BXSUBJFLOW`
  (`rows/bx_t10_flow_pairs.tsv` holds its BX-t10 rows; its other cases are not part
  of this benchmark).
- `cpu_busy.py`: per-CPU use from `/proc/stat`, used to pick idle cores and to count
  busy CPUs.
- `db6_summary.py`: makes the per-case summary (`rows/db6_summary.txt`) from
  `rows/db6_runs.tsv`.
- `s7_time.py`: the `wait4` timer (monotonic wall time, CPU time, max RSS) that the
  drivers call.

`rows/`
- `db6_runs.tsv`: every run of the small-tier run (warmups, timed, disturbed repeats
  and builds); columns `case tool threads run rc wall_s user_s sys_s maxrss_kib rows
  sha16 load busy0 busy1`. Run label `w` is the warmup, `1`-`3` the timed rounds
  (`1d1` and so on a disturbed run that was repeated), `p` the pilot run that chose
  the size of the blastn task query set, and `x` the run on the `LOSAT makedb`
  database (tool `losatdb`).
- `db6_summary.txt`: the output of `db6_summary.py` for it (tables, exactness line
  per case, build medians).
- `db6_followup_runs.tsv`: the follow-up window of 2026-10-17 on build 0b7335ce
  (BP-yeast at 16 threads: `def`, `cut` = `LOSAT_X_BPDBCUT=1`, `ncbi`; BX-t10 at 8 and
  16 threads: `def`, `flow` = `LOSAT_X_BXFLOW=1`, `ncbi`; `timing` = one untimed
  `LOSAT_TIMING` run).
- `bx_t10_flow_pairs.tsv`: the BX-t10 rows of the 2026-10-17 paired run (build 6d0c7b1e).
- `bp_yeast_bpdbcut_pairs.tsv`: the yeast proteome rows (`input` `Yb`) of the
  2026-10-18 paired run (build 9560b7c1; arms in the `tag` column: `new` = default,
  `-BPDBCUT` = `LOSAT_X_BPDBCUT=0`; the `env` column lists the environment, with a
  padding variable `XOSAT_PAD` that only equalizes the size of the environment).

## Not measured

- Human proteome all against all and about 4,000 proteins x Swiss-Prot: hours of
  booked cores for the same question the shorter blastp cases answer.
- tblastn and tblastx: the small tier has no nucleotide genome database. Middle-size
  databases (genomes, gene annotations) are not part of this snapshot.
- LOSAT `-subject` (the FASTA path), and the configuration with the default set
  turned off (`LOSAT_X_DEFAULT=0`): this snapshot covers `-db` with the default.
- WebAssembly: the engine reads a database only natively.
- Database open and read time apart from the search: the database driver has no
  separate timing line; BP-q1 (one short query) bounds the fixed cost.
