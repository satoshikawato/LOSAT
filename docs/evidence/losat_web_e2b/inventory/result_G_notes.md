# Range G result notes (S08 final code, commit 0d533ba76)

Result TSV: `result_G.tsv` (32 rows, one per row of `G.tsv`). Scratch: `rscratch_G/` (inputs, `cmp.sh`, one
`<tag>.{ncbi,losat}.{out,err,rc}` triple per run, `stress_all.log`, `mkstress.py`, `rows_data.py` = the row texts).
LOSAT binary: `s08/native/release/LOSAT` (built from 0d533ba76; `git status` of `LOSAT/src` is clean, so the
source I read is the final code). NCBI oracle: `/home/kawato/micromamba/bin/tblastx` 2.17.0.

## 0. Counts per `s08_result`
ported 9 (rows 2, 3, 4, 6, 9, 15, 16, 17, 28), reused 4 (18, 19, 20, 21), faithful 7 (5, 7, 11, 13, 14, 22, 30),
n/a 8 (10, 23, 24, 25, 26, 27, 29, 31), exception 1 (12), rejected 1 (8), deferred 1 (1), GAP 1 (32), UNSURE 0.

## 1. How the final code implements the hit list (what I read and checked)
* `algorithm/tblastx/report.rs:765 final_hit_order(hits, hitlist_size)` is the whole of range G. Called once from
  `run_impl.rs:1058 write_tblastx_outputs` with `args.max_target_seqs.unwrap_or(500)` on the complete, e-value-filtered HSP
  vector (global `q_idx`, `s_idx` = oid). It groups by (q_idx, s_idx) in a BTreeMap (so oid ascending, queries in input order),
  puts each group in `score_compare_hsps` order (BLAST_LinkHsps end state, link_hsps.c:1802-1810), builds one
  `HitList<TblastxHspList>::new(N)` per query, feeds the subjects with `update()` (= `Blast_HitListUpdate`, `blastn/hsp.rs:1023`),
  then `hit_list.sort_by_evalue()` (= `Blast_HSPResultsSortByEvalue`, `hsp.rs:1096`) and finally sorts the HSPs inside every list
  by e-value (`blast_seqalign.cpp:1574-1577`) while flattening.
* `TblastxHspList` (`report.rs:671-722`) is the `HitListEntry` impl: `best_evalue` = minimum e-value over the HSPs
  (`s_BlastGetBestEvalue`), `first_score` = raw score of `hsps[0]` **in the list's current order**, `sort_by_evalue` =
  `common.rs:380 evalue_compare_hsps`. Because `update()` sorts the lists by e-value only when the query overflows, the score
  vs e-value state of `hsp_array[0]` is reproduced exactly (the inventory's central finding), and the final subject sort
  reads the same state.
* Second hit list of the traceback stage, `BlastHSPStreamClose` order, `s_BlastPruneExtraHits` are not modelled; they are identity
  for ungapped tblastx (prelim size = final size = N, second list never full), see rows 24-26, 29.
* Hit list size: `-max_target_seqs` (>= 1) or 500; outfmt 0 shows `num_descriptions`/`num_alignments` = (N,N) or (500,250)
  (`run_impl.rs:1154-1157`, `report/pairwise.rs:2800-2813`).
* The `-max_target_seqs < 5` warning is `process_options` (`run_impl.rs:759`), called by `search_cli` (:781) and `run_local` (:828).

## 2. Oracle runs (all: stdout + stderr + exit status compared, `rscratch_G/cmp.sh Q S tag fmt [args]`)
Every run below was byte-identical for stdout, stderr and exit status unless stated.
* Exact-tie inputs from the inventory (`tie3.fa`, `tie3b.fa`, `tie3c.fa`, `tie2AB.fa`, `tie2BA.fa`, `zero.fa`, query `Q.fa`), default and
  `-max_target_seqs 1,2,3,5`, outfmt 6 (all), outfmt 7 and 0 (`1,2,3`). The inventory found LOSAT-base wrong on tie2AB/tie3 (default
  `SB SA`, N=1/N=2 cut); the final code gives NCBI's `SA SB`, `SB`, `SB SA`, ... in every case.
* `multi.fa`, `multi_shuf.fa`, `Q2.fa` (2 queries), `Q3.fa` (2 batches): default and N=1,2,3,5 in outfmt 6 and 0.
* `copies600.fa` default (500 of 600 identical copies kept: `C599..C100`), N=7,50,333,499,500,600 in outfmt 6/7/0; `copies300.fa` default 0/6/7;
  `mix600.fa` (600 copies with 0 to 2 percent mutations) default and N=5,40,100.
* Stress: `st1..st3`, `sx4..sx11` (random seeds; 2 to 5 queries cut from LC738874, 30 to 96 subjects that are exact or mutated copies,
  sub-windows, LC738875 windows or random DNA; 16 to 68 subjects with hits per query): default, N=1,2,3,4,5,6,8,10,15,20,25,33 in outfmt 6 and 7,
  N=2 and 7 in outfmt 0. 160 comparisons in `stress_all.log` plus 27 for st1-st3: no difference.
* `-evalue 1e-170,1e-100,1e-60,1e-20,10` with N=1,2,3 on `multi.fa` and `tie3.fa`; outfmt 7 default of `multi.fa`.
* Interleaved too-short subjects (`tie3_short.fa`: 1-, 2-, 3-nt and all-N records between SA/SB/SC), N=1,2,3, default, outfmt 6/7.
* `Qpair.fa` (3 queries) x tie3/tie3b/tie3c/tie2AB x N=1,2,3,default.
* Unsearched/invalid query files of the outfmt0 fixtures x N=1,2,3 x outfmt 6/7/0; `BATCH_SIZE=700` with the batch query, N=2.
* Threads: `-num_threads 4` on tie3/tie3b/tie3c/tie2AB/zero N=1,2,3; `-num_threads 8` on Q2 (outfmt 0, N=3) and copies600:
  stdout identical; stderr differs only by NCBI's `'num_threads' is currently ignored when 'subject' is specified.` (approved exception, row 12).
* Duplicate ids: `dup300.fa` (copies300 with adjacent records sharing an id) default and N=11 in outfmt 0/6/7.
* `-max_target_seqs` values: see row 1.
* Warning text/order: N=1,2,4 stderr; N=5 none; empty query N=2 and N=4 (warning then `Query is Empty!`), empty subject N=2 (exit 3).

## 3. GAP row 32: `-culling_limit` with the hit list (and the pre-existing HSP-culling difference)
NCBI: with `-culling_limit N>0` the preliminary writer is the culling writer (`setup_factory.cpp:330-342, 365-400`): all surviving subjects go into
the query hit lists without a cut, the cut at `-max_target_seqs` happens in the traceback stage, the HSP lists are score-sorted by the culling pipe, so the final subject sort
always reads the score state. LOSAT: `run_impl.rs:3315-3330` applies `hsp_culling::apply_culling` per subject, then the ordinary collector hit list (`final_hit_order`) runs, so
an overflowing query is cut and ordered in the e-value state.

(1) Hit-list difference (isolated: HSP culling removes nothing here):
```
cd rscratch_G
tblastx -query Q.fa -subject tie3.fa -culling_limit 100 -max_target_seqs 2 -outfmt 6   # NCBI
LOSAT tblastx -query Q.fa -subject tie3.fa -culling_limit 100 -max_target_seqs 2 -outfmt 6
```
NCBI: 197 rows, subject runs `SA` (100 rows) then `SB` (97 rows); first row
`Q1	SA	50.407	123	61	0	718	350	752	384	2.67e-167	137`.
LOSAT: 197 rows, `SB` (97) then `SA` (100); first row the same fields with `SB`. md5 66e635176f841bc3c47e9711e59f20f1 (NCBI) vs
c24c5a9f132c8fc2b5fbfd37c45def41 (LOSAT). Files: `tie3_cull_m2.{ncbi,losat}.out`. The same inputs with default N, N=1 and N=3 are identical (tie3), and tie3c/tie2AB are identical for all N
(the tie is decided the same way there); it needs a tie of best e-values where the highest-score HSP is not the first HSP after the e-value sort, and more subjects with hits than N.
(2) HSP-level culling itself differs from NCBI, independent of the hit list and of N (pre-existing; LOSAT-base gives the same numbers):
`tblastx -query Q.fa -subject multi.fa -culling_limit 100 -outfmt 6`: NCBI 489 rows, LOSAT 646 (first difference: LOSAT keeps three W15 rows `...1297 1353 824 880 1.8 23.5`, `...826 882 1.8 23.5`, `...818 874 2.5 23.1` that NCBI culls);
`-culling_limit 2`: 130 vs 254 rows; `tie3.fa -culling_limit 1`: NCBI 60 rows, LOSAT 0 rows; `-culling_limit 2`: 117 vs 126; `-culling_limit 5`: 203 vs 210 (`cl_*`, `multi_cull_*` files).
No inventory range lists the culling algorithm itself (E row 1 only mentions it); I report it here because the culling row of range G is the only one that touches it.

## 4. Things found in the inventory that are wrong or outdated
* Notes section 4 (working tree edits `prune_to_hitlist_size`, `best_score = max raw_score`) is superseded: the final code has no `prune_to_hitlist_size`; the cut is the generic
  `HitList` port as the inventory itself recommended in section 5. Row 28's `divergent` and row 2's description of LOSAT-base are historical.
* Row 6 (`GetPrelimHitlistSize`, `reusable`): the final code does not call `get_prelim_hitlist_size`; it uses N directly (identity for tblastx). Classified `ported`.
* Row 31 remark "two subject records with the same ID next to each other count as one subject": after the hit list at most N subjects remain, so PruneSeqalign cannot cut anything for outfmt 6/7; for outfmt 0
  `dup300.fa` (adjacent equal ids) prints 250 blocks in both programs, i.e. the run counting did not make NCBI print more than 250.
* Row 4: the inventory could not check the order against the subject FASTA diagnostics. The order itself is right in the code (subject read, then `process_options`, then `Query is Empty!`), but the
  FASTA diagnostics are missing in LOSAT TBLASTX altogether (next section), so the combined stderr differs for inputs that trigger them (not a range G defect).
* Stale comment in `run_impl.rs:3941-3945` (quotes the qsort of `Blast_HitListSortByEvalue` above `let hits = final_hits.unwrap_or_default();`, where no sort happens).
* Row 1 domain: besides the clap-vs-USAGE exit code, NCBI accepts hex `0x10`; LOSAT's `positive_usize` rejects it (`ncbi_integer` in `value_parsers.rs:190` exists but is not wired to TBLASTX `-max_target_seqs`).
* Not frozen as fixtures (only the unit test `hit_list_keeps_ncbis_subjects_on_ties` and my runs cover them): the exact-tie inputs `tie3.fa`/`tie2AB.fa`/`zero.fa` and `copies600.fa` default (500 of 600).
  The existing `hitlist.*` fixtures use the 260-subject `many` set.

## 5. Incidental findings outside range G (not counted as G rows; for the owner of the input-reading range)
* FASTA diagnostics: `FASTA-Reader: Ignoring invalid residues at position(s): On line 2: 9-13` is printed by NCBI for a record with digits (`badres.fa`, `badresq.fa`), never by LOSAT TBLASTX
  (stderr differs). For a query with digits NCBI drops them, LOSAT keeps their positions: `tblastx -query badresq.fa -subject multi.fa -outfmt 6 -max_target_seqs 2`: NCBI query coordinates `37 11` / `40 11` / `28 8` / `32 12` (4 rows), LOSAT `42 16` / `45 16` / `37 17` (3 rows); for the same digits in a subject (`-query Q.fa -subject badres.fa`) NCBI prints `Q1 s1 55.556 9 4 0 1733 1707 34 8 0.68 14.8`, LOSAT `Q1 s1 62.500 8 3 0 1733 1710 39 16 2.0 13.4`.
  Range F row "CTrans_table::x_InitFsaTable" / the FASTA-reader rows are the likely owners.
* A subject file without any `>` line (`bad.fa` = `not fasta`): NCBI reads it with `FASTA-Reader: Ignoring invalid residues at position(s): On line 1: 2, 5`, exit 0; LOSAT exits 1 with
  `Error: failed to read subject FASTA bad.fa / Caused by: Expected > at record start.` (stdout identical, empty). Probably part of the deferred input-file error/"no data" items.

## 6. "What the port must get right" (inventory section 6) against the final code
1. N per query, prelim size N, only subjects with an HSP after the final `-evalue` reap: done (rows 2, 6, 14, 15, 17).
2. Kept set = N best by (e-value class with 1e-180, first-HSP score in the e-value-sorted list, later oid first): done by the `HitList` port (rows 18-21); 462 NCBI-vs-LOSAT runs in `rscratch_G` (`*.ncbi.rc` files), stdout identical in all except the `-culling_limit` runs (section 3), `-max_target_seqs 0x10` (row 1) and the two digit-FASTA runs (section 5).
3. `hsp_array[0]` state: score order unless the query overflowed, then e-value order, also for the kept lists' subject order: done (row 28, tie runs).
4. stderr warning for N < 5 in every format: done (row 4).
5. outfmt 0: descriptions up to N (500), alignment blocks only up to 250 unless `-max_target_seqs`: done (row 3).
Not done: `-culling_limit` (row 32).
