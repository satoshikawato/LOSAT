# BLASTP unified P-values: the last query of a query batch (PD-LOSAT-NCBI-DEFECTS exception 16)

Line RG item X5, 2026-11-20. NCBI BLAST+ 2.17.0 (the official Linux x86-64 build) and LOSAT, `blastp
-comp_based_stats 1u`, `2u` and `3u`.

## What the runs show

NCBI scales each e-value of a query with the length, length adjustment and search space of the
**next** query of its batch (`blast_kappa.c:3543`, `:3590-3592`, `:3675-3683`, `:134-210`). For the
last query of a batch, and so for every search of one query, it reads those values past the end of
its array of query contexts (heap bytes), and its report of that query is not fixed. A copy of the
query after it in the same batch gives NCBI the query's own values, and its report of the query is
then the same in every run. LOSAT writes that report.

`probe.sh` runs each case and records, per run, the SHA-256 (12 digits) and number of the outfmt 6
rows of the named query (`results_*.tsv`: case, run, threads, SHA-256, rows).

| Case | Query and search | NCBI | LOSAT |
|---|---|---|---|
| A | `BDT63074.1` (897 residues) alone, `2u -evalue 1000` | 7 different reports in 40 runs (2,127, 2,129 and 2,132 rows): run 1 gave one report in its 10 runs at 1 thread and 6 in 10 at 4 threads; run 2 gave one report in all 20 (one of run 1's at 4 threads) | - |
| B | the same query with a copy after it | one report (1,944 rows) in every run: 1 thread (3 runs), 4 and 8 threads, and with a third query after the copy | - |
| C | `BDV03099.1` (104 residues) before `BDT63074.1` | one report (783 rows) at 1 and 8 threads, the same with another protein of the next query's length; with its own copy as the next query another report (782 rows) | - |
| D | case A and B without the `u` (mode 2) | the copy changes nothing (382 rows, the same SHA-256) | - |
| E | LOSAT on cases A and C | - | case A = NCBI's report of case B; case C = NCBI's at 1 and 8 threads |
| F | `BDV03099.1` before a query of 50 `X` (no search space) | every e-value 0.0; with `-searchsp 100000000` a defined report | refused (exit 1); with `-searchsp` = NCBI's report |
| G | `BDT63114.1` (86 residues) alone against `tests/db_search/paja`, `-max_target_seqs 1`: a full hit list | 3 different reports in 10 runs (8, 1, 1); with a copy after it one report in 10 runs | = NCBI's report with the copy at 1 and 8 threads |

Case G is the input of decision X5-4: NCBI's redo ends matches early only when every hit list of the
batch is full (`redo_alignment.c:1560-1580`, `compo_heap.c:405-409`), and it redoes the matches in
query order (`blast_hspstream.c:288-319`), so with a copy after the query, whose list is still empty,
no match of the query is ended early. LOSAT with unified P-values ends no match early.

## The runs

- `results_1.tsv`: 2026-11-20 16:59 JST, ws (2 x Xeon, Linux), `probe.sh` of commit d29a4a6ca (cases A
  to G), LOSAT of commit d98e4389d.
- `results_2.tsv`: 2026-11-20 17:07 JST, the same machine, `probe.sh` of commit d98e4389d (cases A to F;
  case G was added after it), LOSAT of commit d98e4389d.

The two runs differ only in case A, whose reports change from run to run.

```sh
bash docs/evidence/unified_p_last_query/probe.sh <NCBI BLAST+ 2.17.0 bin> <this repository> <work dir> [losat]
```

Every run is capped (`prlimit --core=1 --as=4000000000 timeout 300`). The fixture rows `db.x5u.*`
(`LOSAT/tests/db_search_manifest/blastp.x5u.tsv`, `psiblast.tsv`) and `x5u.*`
(`LOSAT/tests/outfmt0_manifest.tsv`) hold NCBI's frozen outputs of searches of the same kinds (case
G's search is `db.x5u.cp.c.mts1.2u.6`; outfmt 0, 5, 6, 7 and 11; `-max_target_seqs`; `BATCH_SIZE`;
`-searchsp`; `-dbsize`; `-use_sw_tback`; the refusals).
