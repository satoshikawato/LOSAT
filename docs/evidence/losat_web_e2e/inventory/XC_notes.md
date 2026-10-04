# Range XC notes: TBLASTX `-culling_limit` and the hit-list trimming that moves with it

Table: `XC.tsv` (35 rows). Scratch: `scratch_XC/`. NCBI paths are relative to `c++/` (pinned commit 598d8ae6), LOSAT paths to `LOSAT/src`
(base copy 78c06fe61). NCBI oracle: `/home/kawato/micromamba/bin/tblastx` 2.17.0.

## 0. Result in one paragraph

LOSAT rejects `-culling_limit N>0` for TBLASTX (`algorithm/tblastx/args.rs:873-880`). The rejection can be lifted: the whole NCBI path was
traced and re-implemented as an executable specification (`scratch_XC/cullsim.py`, Python, a line-by-line transliteration of
`hspfilter_culling.c` plus the stage glue). It reproduces NCBI's `tblastx` output EXACTLY (row set and row order, outfmt 5 compared) in every one of
about 200 `(culling_limit, -max_target_seqs, other option)` combinations on 5 input sets, including exact-tie subjects (identical copies), reverse
complements, several queries, `-evalue 1000/1e-10`, `-seg no`, `-sum_stats false`, `-query_gencode 4`, `-window_size 0`, `-num_threads 4`, `-max_hsps 1..3`
and the Int4 wraparound values. So a faithful port is "transliterate cullsim.py to Rust at the place described in section 4".

NCBI's culling is NOT one filter at one place. It is (a) a preliminary-stage writer that culls across ALL subjects with merit N+3 and does not
cut the hit list, (b) the generic hit-list cut at `-max_target_seqs` in the traceback stage, and (c) a second culling pass (pipe) with merit N over
the hit list that survived the cut. LOSAT's `hsp_culling.rs` is a good transliteration of the interval tree but is wired as one pass per subject,
with a self-domination bug, the wrong query length, the wrong HSP order and the e-value filter after the culling.

## 1. Call path followed (tblastx bl2seq, `-query Q -subject S -culling_limit N`)

1. `CTblastxApp::Run` (`app/blast/tblastx_app.cpp`) -> `CLocalBlast lcl_blast(queries, opts_hndl, db_adapter); lcl_blast.Run()`.
2. Arguments: `CHspFilteringArgs::SetArgumentDescriptions` (`blast_args.cpp:3291-3331`: `-culling_limit` `CArgInteger`, constraint `>= 0`, `eExcludes`
   `-best_hit_overhang`/`-best_hit_score_edge`; also in tblastx via `tblastx_args.cpp:87`) -> `ExtractAlgorithmOptions` (`:3335-3340`) ->
   `CBlastOptions::SetCullingLimit` (`blast_options_cxx.cpp:1837`) -> `CBlastOptionsLocal::SetCullingLimit` (`blast_options_local_priv.hpp:1319-1341`):
   `N<=0` returns; else `hsp_filt_opt->culling_opts->max_hits = N`, `culling_stage = eBoth`, `culling_limit = N`.
3. `SetupInternalData` (`blast_aux_priv.cpp:225-240`): `CSetupFactory::CreateHspWriter` (`setup_factory.cpp:296-362`) picks the culling writer
   (`BlastHSPCullingParamsNew`; `if (culling_max > 1) culling_max += 3`), and `CreateHspPipe` (`:365-403`) registers a culling PIPE for the
   traceback stage (`culling_max = N`, no +3).
4. Preliminary stage per subject (`BLAST_PreliminarySearchEngine`, `blast_engine.c:1348-1554`; ascending oid, single thread for `-subject`):
   `s_BlastSearchEngineCore` (`:702-915`) -> `BLAST_LinkHsps` (ends with `Blast_HSPListSortByScore`, `link_hsps.c:1761-1810`),
   `Blast_HSPListReevaluateUngapped`, relink, `s_Blast_HSPListReapByPrelimEvalue` (prelim_evalue = `-evalue`, `blast_parameters.c:868,933`),
   `Blast_HSPListGetBitScores`, then `BlastHSPStreamWrite` (`blast_hspstream.c:339`) -> first call `s_BlastHSPCullingInit`, then
   `s_BlastHSPCullingRun` per subject list (all queries' contexts in one list, score order).
5. `Blast_RunTracebackSearchWithInterrupt` (`blast_traceback.c:1809`) -> `BlastHSPStreamClose` (`blast_hspstream.c:133`) -> `s_FinalizeWriter` ->
   `s_BlastHSPCullingFinal` (builds the per-query hit lists from the trees, score-sorts every list, NO size cut).
6. `BLAST_ComputeTraceback_MT` ungapped branch (`blast_traceback.c:1495-1735`): `BlastHSPStreamToHSPStreamResultsBatch` (ascending oid), per list
   `Blast_HSPListGetBitScores` + `Blast_HSPResultsInsertHSPList(results, list, hitlist_size)` = `Blast_HitListUpdate` (`blast_hits.c:3243`): THE cut at `-max_target_seqs`.
   Then `BlastHSPStreamTBackClose` (`:1721`) = `s_BlastHSPCullingPipeRun` (`hspfilter_culling.c:699`), then `s_FilterBlastResults`
   (`-max_hsps`, `-qcov_hsp_perc`, `-subject_besthit`, `:1763-1766`), `Blast_HSPResultsSortByEvalue` (`:1768-1772`), `s_BlastPruneExtraHits`.
7. Seq-align build sorts the HSPs of each subject by e-value (`blast_seqalign.cpp:1574-1577`); formatting is unchanged (outfmt 0 prints
   `Expect(n)` with the `num` of the sum-statistics linking, also for HSPs whose partners were culled).

## 2. The algorithm that must be ported (what cullsim.py does)

Per HSP: `cid` = context (query index x 6 + frame index, tblastx; the `isBlastn` branch is not taken), `sid` = oid, `begin/end` = `hsp->query.offset/end`
= frame-relative amino-acid offsets (LOSAT `sort_query_offset/sort_query_end`, `q_aa_start/q_aa_end`), `score` = the HSP's own raw score,
subject offset = frame-relative (`sort_subject_offset`). One tree per context, root `[0, qlen]` with `qlen` = frame length in amino acids
(`query_info->contexts[ctx].query_length` = `(L - (|frame|-1)) / 3`, LOSAT `QueryContext::aa_len`). `s_DominateTest`, `s_FullPass`, `s_SaveHSP`,
`s_ProcessHSPList`, `s_ProcessCTree`, `s_MarkDownCTree`, `s_ForkChildren` (threshold 20) as in `hspfilter_culling.c`; LOSAT's `hsp_culling.rs`
already has them (faithful) except for the points in the table.

Stage sequence (this is what `cullsim.py:simulate_full` does and what was validated):

1. Input: the complete HSP vector after all subject jobs and after the `-evalue` filter (`run_impl.rs:4236-4262`), grouped per `(query, subject)`,
   each group in `Blast_HSPListSortByScore` order (`common::score_compare_hsps`; `final_hit_order` already does this).
2. Prelim pass: `culling_max = if N > 1 { N.wrapping_add(3) } else { N }` (i32). Subjects in ascending oid; HSPs in that order into the tree of their context.
   `merit` is i32 and decremented with `wrapping_sub`.
3. Final: per context ascending: rip (node list, left, right), regroup per `(query, oid)`, score-sort each group. No cut.
4. Hit list: per query `HitList::new(max_target_seqs or 500)`, `update()` with the groups in ascending oid (existing code; the same `HitList` as the default path).
5. Pipe pass: per query: HSPs of every kept list sorted by e-value (`s_EvalueCompareHSPs`: e-value with 1e-180 epsilon, then `ScoreCompareHSPs`), list
   `best_evalue = hsp[0].evalue`, lists sorted by `s_EvalueCompareHSPLists` (best e-value with epsilon, then `hsp[0].score` higher first, then LARGER oid first);
   then every HSP in that order into NEW trees with `merit = N` (no +3); Final as in 3.
6. Existing tail: `sort_by_evalue` of the hit list (reads `hsps[0]` = best-SCORE HSP because the lists are score-sorted after step 5), prune to N (no-op),
   HSPs of each subject by e-value. (With `-max_hsps K` from range XA: trim each score-sorted list to K before the subject sort, `best_evalue` stays stale; verified.)

Everything is plain integer arithmetic (Int8 in `s_DominateTest`) plus the existing comparators. No NCBI data tables are needed.

## 3. Evidence and fixture candidates (`scratch_XC/`)

* Inputs (cut from LC738874/LC738875 with `mk.py`): `q_0.fa`/`s_0.fa` (3 kb windows, 58 HSPs), `q_200000.fa`/`s_200000.fa` (33), `Qm.fa` (3 queries) +
  `Sm.fa` (7 subjects: window of LC738875, an identical copy of the query window, 8 % and 5 % mutated copies, random DNA), `Q6.fa`/`S6.fa` (6 kb, 432 HSPs at `-evalue 1000`),
  `Qt.fa` (2 queries) + `St.fa` (10 subjects incl. 4 identical copies c1..c4, reverse complements, mutated copies; exact ties for the hit list).
* `runfix.sh` writes 98 NCBI outputs to `fixtures/` (`<tag>.out/.err/.rc`; md5 in `fixtures.md5`, line counts in `fixtures_summary.tsv`). Key counts, outfmt 6:
  `q_0/s_0`: cl 0/1/2/3/5/10 = 58/20/31/35/39/44. `Qm/Sm` N=500: cl 0/1/2/5/10/20 = 427/101/178/274/367/422; N=1: 217/98/136/182/214/217.
  `Qt/St` N=3: cl 1/2/3/7 = 85/172/231/270. Also outfmt 0 and 7 (`*_f0`, `*_f7`), `-evalue 1000/1e-10`, `-max_hsps 1,2` (`*_mh*`), the big values
  (`*_cl2147483644..47`), argument spellings (`t1_clraw*`, `t1_clbad*`: error text in `.err`), best_hit exclusion (`t1_cl*_bh*`), `-subject_besthit` (`t1_cl2_sbh`: 27 rows).
  The `.err` files of runs with `-max_target_seqs` < 5 carry NCBI's `Examining 5 or more matches is recommended` warning.
* `cullsim.py` (library) + `fullcheck.py Q S cl_list N_list` (env `EXTRA="-evalue 1000"`, `MAXH=2`) compares the simulation to NCBI outfmt 5 (row set AND row order).
  `ablate.py`, `ablate2.py` switch single NCBI mechanisms off and show what each one changes (results in `ablate*_*.txt`):

| switched off / changed | what it models | example effect (rows vs NCBI) |
|---|---|---|
| no second (pipe) pass | LOSAT has one pass | q_0/s_0 cl=2 +8, cl=5 +3; Qm/Sm cl=2 +95 |
| per-subject trees | `apply_culling` per subject | Qm/Sm cl=1 +94, cl=2 +90, cl=5 +98 |
| self-domination (clone as `x`) | `x_for_dom` | q_0/s_0 cl=1 0 rows (NCBI 20); Qm/Sm cl=1 0 (101); cl=5 -23 |
| no +3 in the prelim pass | N used as is | Qm/Sm cl=5 -1; `-evalue 1000` cl=2 +7/-9, cl=5 +9/-17 |
| qlen = nucleotide length | `ctx.orig_len` | Q6/S6 `-evalue 1000` cl=1 +1; Qm/Sm cl=5 +2/-4 |
| HSP input order reversed / by e-value | not score order | Qt/St `-seg no` cl=20 +40/-83 (+28/-30) |
| fork keeps parent order | LOSAT `fork_children` | 0 differences (equivalent) |
| all LOSAT mechanisms together (per-subject + nt qlen + self-dom, 1 pass, no +3) | the base code if enabled | Qm/Sm cl=2 +28/-11, cl=5 +67/-1; Q6/S6 cl=2 -87 |

* LOSAT: the base binary only reports `Error: -culling_limit N is not supported by LOSAT's TBLASTX ...` (exit 1) for N>0; no LOSAT run of the culling path is possible without a build. S08 measured the
  old code (`result_G_notes.md` section 3: `tie3.fa -culling_limit 1`, NCBI 60 rows, LOSAT 0; `multi.fa -culling_limit 100`: 489 vs 646) which the table above explains.

## 4. Where the port goes in LOSAT

* Remove the check at `args.rs:873-880`; parse `-culling_limit` with `value_parsers::ncbi_integer` (i32, `>= 0`). Keep `culling_limit` as i32 in `TblastxArgs`
  (`web_api.rs:364` builds `TblastxArgs { culling_limit: 0 }` and its tblastx extra-argument parser (`web_api.rs:366-470`, no `-culling_limit` branch) must learn `-culling_limit`; `web/adapter/src/run.rs:396` has a test that
  expects `-culling_limit 2 is not supported`, and `docs/web/abi_v2.md` (validate row) lists it as rejected: update both).
* Delete the call at `run_impl.rs:3629-3640` (per-subject culling before the e-value filter). New function next to `report::final_hit_order` (called from
  `write_tblastx_outputs`, `run_impl.rs:1214`) that implements steps 1-6 of section 2 on the collected, e-value-filtered `Vec<TblastxHsp>`. It needs per
  HSP: `hit.q_idx`, `hit.s_idx`, `hit.query_frame` (context), `hit.sort_query_offset/_end`, `hit.sort_subject_offset`, `hit.raw_score`, `hit.e_value`, and per query the
  frame aa lengths (`(L - (|f|-1)) / 3`; `run.queries` has the lengths). `hsp_culling.rs` keeps its tree code; change `process_hsp_list` to skip the entry by identity
  (index), take `merit: i32` with wrapping ops, and replace `apply_culling` by a state object (`CullingWriter { trees, culling_max }` with `run(list, oid)` and `finalize()`).
* The default path (culling_limit = 0) must stay byte-identical; `final_hit_order` keeps its current behavior when `culling_limit == 0`.
* Order dependence: the survivors depend on the feed order (feeding the reverse order changes 6-83 rows). Always walk subjects in ascending oid and each
  `(query, subject)` group in `score_compare_hsps` order, even when the subject jobs ran in parallel (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES 2 covers only the missing warning).

## 5. Surprises

* `eBoth`: `-culling_limit` culls twice. The prelim pass uses N+3 (only if N > 1), the pipe uses N. Skipping either one changes the output.
* The hit list is NOT cut in the preliminary stage with culling. With `-max_target_seqs 1` NCBI still keeps all subjects that have a surviving HSP until the traceback
  stage; the cut is made among the CULLED lists (so the best subject is judged after the N+3 culling), and the second pass then culls only inside the kept subjects.
* `-culling_limit 2147483645` is a no-op, `2147483646` and `2147483647` behave like `-culling_limit 1` (signed overflow of `culling_max + 3` and of the merit; the shipped
  binary wraps). 2147483648 and above is an argument error. Per the NCBI defect policy this is deterministic and should be reproduced (port with wrapping i32).
* `-culling_limit 0 -best_hit_overhang X` is an NCBI argument error (the exclusion tests presence). `-subject_besthit` is not excluded.
* The mark-down shortcut (`x` covers a whole node range -> every HSP below loses merit without a score test) makes the result depend on the TREE SHAPE, i.e. on `qlen`
  and the midpoints; that is why passing the nucleotide length changes rows.
* The `-evalue` reap happens before the writer; sum-statistics cutoffs depend on `-evalue`, so the pre-cull HSP set differs per `-evalue` (the NCBI-vs-sim comparisons always
  feed the sim with the same `-evalue` run).
* `Expect(n)` in outfmt 0 keeps the linking `num` after culling (e.g. `Expect(13)` with fewer HSPs displayed); LOSAT's `TblastxHsp.num` already does this.

## 6. Decisions for the session

1. Port (recommended; deterministic, validated, M-sized: tree code exists) or keep rejecting. Nothing found in NCBI is a crash or undefined read except the benign Int4 wraparound.
2. Reproduce the Int4 wraparound for `-culling_limit` 2147483645..2147483647 (row 6) or reject those three values explicitly (they are meaningless inputs; reproducing costs nothing).
3. `-best_hit_overhang/-best_hit_score_edge`, `-max_hsps`, `-subject_besthit`, `-qcov_hsp_perc` are range XA; only their order relative to the culling is fixed here (rows 2, 30, 32).
4. The rows' `impact` is `medium` for the algorithmic rows because the option itself is a non-default value; the defects change the output for every N.
