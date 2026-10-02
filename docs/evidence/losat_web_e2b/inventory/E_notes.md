# Range E notes: what the tblastx report receives from the search

Files: TSV `E.tsv` (60 rows). Scratch inputs/outputs: `scratch_E/` (all oracle files named below live there).
NCBI line numbers are CR-stripped, pinned commit 598d8ae6. Oracle = `/home/kawato/micromamba/bin/tblastx` 2.17.0.

LOSAT caveat: the LOSAT working tree was being edited while I read it (uncommitted `TblastxHsp`, `TblastxQueryStats`,
`TblastxBatch`, `num`/`hsp_link_num` on `UngappedHit`, query batching in `run_in_pool`, generic `SubjectGroup<T>` in
`common.rs`). Every "LOSAT" statement below is for the tree as it was when I last read each file; line numbers drift,
function names are stable.

## 0. Useful oracle trick (not in the allowed option list, so only for inspection)

`tblastx -outfmt 8` prints the Seq-align set as ASN.1 text and `-outfmt 11` prints the BLAST archive, which contains
the masks (`masks {... frame plus1 ...}`), the ka-blocks, and `search-stats` ("Effective search space", "Length
adjustment") exactly as the formatter receives them. They were used to confirm every row below. Parsers:
`scratch_E/parse8.py` (Seq-align -> python), `scratch_E/parse0.py` (outfmt 0 HSP blocks), `scratch_E/kab.py`
(independent model of `Blast_ScoreBlkKbpUngappedCalc`; it reproduces NCBI's lambda to 5 digits).

## 1. Call path of the data (tblastx -query Q -subject S)

1. `tblastx_app.cpp:176-207`: loop over query batches (`CBlastInput`, batch = cumulative query length >= 10002 nt,
   `blast_input_aux.cpp:130-134`, `blast_input.cpp:134-176`). One `CLocalBlast::Run()` per batch.
2. `blast_app_util.cpp:206-210`: `-subject` => `CLocalDbAdapter(subjects, opts, true)` = DbScanMode (unless env
   BL2SEQ_LEGACY). Therefore (a) the subject set is treated as a database: total length and sequence count feed the
   statistics; (b) `local_blast.cpp:288-291` does NOT set `eSequenceComparison`; the result set is a normal database
   result set: **one `CSearchResults` per query**, all subjects of that query inside one `CSeq_align_set`.
3. Prelim stage: `SetupInternalData` -> `CreateScoreBlock` -> `BLAST_MainSetUp` (`blast_setup.c:563-667`):
   SEG per translated frame, mask applied (X=21) to the query buffer, masks converted to nt (`BlastMaskLocProteinToDNA`),
   THEN `BlastSetup_ScoreBlkInit` computes the per-frame Karlin blocks from the masked frames. Masks are turned into
   `TMaskedQueryRegions` by `Blast_GetSeqLocInfoVector` inside `CreateScoreBlock` (`setup_factory.cpp:187-194`).
4. Search + sum statistics (engine range): per subject `BLAST_LinkHsps` ends with `Blast_HSPListSortByScore` and sets
   `best_evalue` (`link_hsps.c:1802-1810`); bit scores are set at `blast_engine.c:1541`; the list is written to the HSP
   stream (`Blast_HitListUpdate`, `blast_hits.c:3243-3299`).
5. Traceback stage for ungapped tblastx (`blast_traceback.c:1556-1696`): no subject fetch, no
   `BLAST_OneSubjectUpdateParameters`; only `Blast_HSPListGetBitScores`. Then (`blast_traceback.c:1768-1777`)
   `Blast_HSPResultsSortByEvalue` (TotLen>0 in DbScanMode) and `s_BlastPruneExtraHits(hitlist_size)`.
6. `LocalBlastResults2SeqAlign` -> `s_BlastResults2SeqAlignDatabaseSearch_OMF` -> `BlastHitList2SeqAlign_OMF` per query
   -> per subject: `Blast_HSPListSortByEvalue`, then `BLASTUngappedHspListToSeqAlign` (gapped=false).
7. `BlastBuildSearchResultSet` -> `BuildBlastAncillaryData` (one `CBlastAncillaryData` per query, built from the
   PRELIM-stage `sbp`/`query_info`, after the whole search) -> `CSearchResultSet`.
   `CLocalBlast::Run` then calls `SetFilteredQueryRegions` (`local_blast.cpp:294`).
8. Formatter (`blast_format.cpp:1411-1590`, `x_PrintTabularReport:759-836`): `PrepareBlastUngappedSeqalign`, then
   outfmt 0 / 6 / 7 printing, footer from the ancillary data.

## 2. The Seq-align of one subject (oracle: `qS8_vs_sS8` = `scratch_E/qS8.o8`)

* One `CSeq_align` per subject: `type diags`, `segs std { ... }`, one Std-seg per HSP, `dim 2`, `ids {Query_1, Subject_1}`,
  two `int` Seq-locs (query first) with `from`, `to`, `strand`, `id`. No gaps, no merged HSPs, nothing removed.
* nt coordinates (0-based inclusive): plus frame `from = 3*off + frame - 1`, `to = 3*end + frame - 2`; minus frame
  `from = L - 3*end + frame + 1`, `to = L - 3*off + frame`; L = nucleotide length of that sequence
  (`blast_seqalign.cpp:1345-1377`). Tabular = from+1 / to+1, reversed for minus (LOSAT `convert_coords` equals this).
* Score list inside each Std-seg, in this order (`blast_seqalign.cpp:1162-1228`):
  `score` (int), `blast_score` (int), [`sum_n` (int) only if `hsp->num > 1`], `e_value` (real; `sum_e` when num>1;
  stored as 0.0 if < 1e-180), `bit_score` (real), `num_ident` (int, always), [`comp_adjustment_method` never for
  tblastx], `num_positives` (int, if >0), `hsp_percent_coverage` (real, uses the NUCLEOTIDE query length).
  The 0/6/7 formatters read only `score`, `bit_score`, `e_value|sum_e`, `sum_n`, `num_ident`
  (`align_format_util.cpp:230-270`).
* `hsp->num`: `BLAST_LinkHsps` sets every num to 1 (`link_hsps.c:1773-1776`) then each chain member gets the chain
  size (`link_hsps.c:1049,1061`). All members of a chain share `sum_e`.
* The formatter then calls `PrepareBlastUngappedSeqalign` (`showalign.cpp:3162-3213`): subjects with >1 Std-seg are
  split into one Seq-align per HSP (scores moved to Seq-align level), a subject with exactly one HSP is kept as is.
  "# N hits found" (outfmt 7) is the number of those Seq-aligns before any -num_alignments/hitlist prune.

## 3. What decides the order

HSPs inside a subject (Seq-align and therefore outfmt 0, 6 and 7): `Blast_HSPListSortByEvalue`
(`blast_seqalign.cpp:1577`): `s_EvalueCompareHSPs` = e-value ascending (`s_EvalueComp`: two values below 1e-180 are
equal), then `ScoreCompareHSPs`: raw score DESC, subject.offset ASC, subject.end DESC, query.offset ASC, query.end DESC.
All four offsets/ends are the frame-relative protein BlastSeg values, the frame/context is not compared. The list is
sorted only if the pre-check finds an inversion (then glibc `qsort`).

Equal raw score in two frames: the frame's own Karlin block decides the e-value (`BLAST_LinkHsps` uses
`kbp[hsp->context]->Lambda/logK` for the sum, `eff_searchsp`/lengths of the first HSP's context of the strand group,
`link_hsps.c:560-571,913`). With ideal blocks in every frame the e-values of single HSPs are practically identical
(they cancel the search space), so the order falls to the offsets. With a biased frame the block differs.
Oracle `scratch_E/qE.fa` vs `scratch_E/sC.fa`, `-seg no`:

| Seq-align row | raw | e-value | bits | Frame | note |
|---|---|---|---|---|---|
| 0 | 519 | 9.86e-68 | 240.7 | +3/+1 | normal frame (ideal block) |
| 1,2 | 518 | 1.38e-67 | 240.3 | -1/-1, -2/-1 | |
| 8 | 519 | 8.85e-61 | 217.6 | +1/+1 | frame +1 has lambda 0.286, K 0.0978 |

Row 8 has the same raw score as row 0 but is listed after rows with raw 518, 488, 479 because its e-value is larger.
outfmt 6 (`qE_sC_segno.o6`) and outfmt 0 (`qE_sC_segno.o0`) have exactly this order. The footer of that run is
`0.286 0.0978 0.222` (first valid context = frame +1, which keeps its own block) and search space 142506; the same
query with default SEG prints `0.318 0.134 0.401` and 145180 because the masked frame is replaced by the ideal block.

Subjects inside a query: hit-list order = `Blast_HSPResultsSortByEvalue` (`blast_traceback.c:1771`, because DbScanMode has
TotLen>0): `s_EvalueCompareHSPLists`: best_evalue (min over HSPs) ASC, then `hsp_array[0]->score` DESC, then oid DESC
(`blast_hits.c:3071-3107`). The Seq-align builder iterates the hit list, so input order is irrelevant. Oracle:
`sDup.fa` (dupA, dupB, dupC identical) prints dupC, dupB, dupA in outfmt 0, 6, 7; `sMulti.fa` (444, 7000, 438 nt)
prints 444, 438, 7000.

`hsp_array[0]->score` at that moment is the MAXIMUM raw score of the subject: `BLAST_LinkHsps` finishes with
`Blast_HSPListSortByScore` (`link_hsps.c:1803`) and nothing re-sorts the list before step 5 (unless the hit list had
overflowed and `Blast_HitListUpdate` sorted by e-value at `blast_hits.c:3273,3284`). LOSAT (`common.rs`
`ncbi_query_subject_groups`, EvalueCompare mode) takes the score of the first HSP AFTER the e-value sort. These agree
except when two subjects tie on best_evalue (usually both 0.0) and the best-e-value HSP is not the max-score HSP.
Row `link_hsps.c 1802-1810` (divergent, low).

`-max_target_seqs N` (hitlist_size): the heap in `Blast_HitListUpdate` keeps the N best subjects under the same
comparator; "newer hits equal to the worst replace it" so with identical subjects the highest oids survive
(oracle: N=1 -> dupC only). stderr: `Warning: [tblastx] Examining 5 or more matches is recommended` for N < 5.
LOSAT TBLASTX parses `max_target_seqs` (`algorithm/tblastx/args.rs:57`) but nothing applies it (row `missing`).

Comparator-equal HSPs (same e-value, score, frame-relative offsets/ends, different frames) occur for repeats:
`scratch_E/qPA.fa` (poly-A with a CGTACGT spacer) vs `scratch_E/sPA.fa`, `-seg no` gives 1440 HSPs; rows 4-6 of
`qPA.o8` are query 100-186, 99-185, 98-184 (frames -1,-2,-3) with identical keys. Their relative order is the incoming
list order; LOSAT replays with a stable sort over `hsp_list_order`. Good regression input (UNSURE row).

## 4. Ancillary data (footer) per query

`CBlastAncillaryData(program, query_number, sbp, query_info)` (`blast_results.cpp:72-115`): first context
(+1,+2,+3,-1,-2,-3) with `is_valid`; takes `eff_searchsp` (and `length_adjustment`) of that context, then
`kbp_std[ctx]` (ungapped block), `kbp_gap` (NULL for ungapped tblastx), `kbp_psi*` (filled but only printed by
psiblast/deltablast), `gbp` (freed for ungapped, NULL). `s_InitializeKarlinBlk` copies only when `Lambda >= 0`.

Block value: `Blast_ScoreBlkKbpUngappedCalc` (`blast_stat.c:2736-2829`): composition of the masked frame (stop codons
counted, only X ignored) against Robinson background, `Blast_KarlinBlkUngappedCalc`; context invalid if that fails;
`check_ideal` (tblastx, blastx, rpstblastn): if computed `Lambda >= ideal Lambda` the ideal block replaces it.
Ideal = 0.3176/0.1340/0.4012 (printed `0.318 0.134 0.401`). Stop codons raise lambda a lot, so almost every real frame
ends at the ideal block; only stop-free biased frames keep their own block. My first model that ignored stops wrongly
predicted non-ideal blocks.

Search space: `BLAST_CalcEffLengths` (`blast_setup.c:699-850`), computed once per batch in `BLAST_GapAlignSetUp`
with the TOTAL subject length/3 and the NUMBER of subjects (DbScanMode); per-subject update
(`blast_engine.c:1434-1444`, `blast_traceback.c:1598-1630`) does not run. Oracle: qS8 (324 nt) vs three subjects
(444+7000+438 nt): `Length adjustment: 22`, `Effective search space used: 220246` = (108-22)*(2627-3*22); each
subject alone gives 12369 / 198746 / 12183. The footer is printed once per query after all subjects.
LOSAT: `compute_eff_lengths_tblastx` + `eff_searchsp_per_context` (run_impl.rs) follow this and the working tree already
collects the first valid context into `TblastxQueryStats` (karlin, eff_searchsp, seg_masks).

Three distinct "invalid" situations (bytes from `q2.o0`, `q4.o0`, `multi.o0`):

1. Whole batch invalid (`BlastScoreBlkCheck` != 0; e.g. a lone `AC`, all-N, only-stop query): unsearched path
   (`local_blast.cpp:177-225`). Per query: NULL align set, ancillary (-1,-1,-1, 0) with BOTH blocks non-NULL, no masks,
   stderr warning (a WARNING, so exit status 0). Report footer after "***** No hits found *****\n\n\n":
   `\nLambda      K        H\n   -1.00    -1.00    -1.00 \n\nGapped\nLambda      K        H\n   -1.00    -1.00    -1.00 \n\nEffective search space used: 0\n`.
   outfmt 7: header only, no "# 0 hits found" (NULL set). outfmt 6: nothing. stderr, per query, when its result is
   formatted: `Warning: [tblastx] Query_N <title>: Could not calculate ungapped Karlin-Altschul parameters due to an
   invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options ` (trailing
   space). Batch = cumulative length >= 10002 nt, so `qBig(10100 nt) + q2 + qN` gives warnings only for Query_2 and Query_3.
2. Query without any valid context inside a batch that has valid contexts (`qSmallplus2.fa`): searched; its set is the
   empty non-NULL set; ancillary blocks NULL and search space 0; footer = `\n` `\n` `\n` + `Effective search space used: 0`
   (no Lambda table, no Gapped); outfmt 7 prints `# 0 hits found`; NO stderr warning (the per-context warning is
   written only for non-translated programs, `blast_stat.c:2787-2790`).
3. Valid query with some invalid frames (`q4.fa` = `ACGT`: frames +1,+2 have one residue, frame +3 none): footer from the
   first valid frame: `0.278 0.0810 0.188`, `Effective search space used: 148` (= 1 residue * 148 aa). With stops-only
   frame +1 (`qStop.fa`, `-seg no`) the footer comes from frame +2: `0.278 0.108 0.280`, space 3645, length
   adjustment 13. With default SEG `qStop.fa` is fully masked -> case 1.

## 5. Query masks

* SEG (`BlastSetUp_Filter`, window/locut/hicut 12/2.2/2.5, `overlaps=TRUE`) runs per valid frame on the protein
  buffer; intervals inclusive; `BlastSeqLocCombine(.,0)` merges only `previous.right > next.left`
  (`blast_filter.c:971-1011`; LOSAT `combine_masks` is identical); no reversal for protein frames.
* Mask applied as X (21) to the working buffer; `sequence_nomask` kept; Karlin blocks computed afterwards.
* `BlastMaskLocProteinToDNA` (`blast_filter.c:892-953`) then rewrites each frame's intervals in nt coordinates with the
  quirky end rounding: plus frame `[3l+f-1, 3r+f-1]`, minus frame `[L-3r+f+1, L-3l+f]` (so a minus mask covers 3(r-l) bases;
  a single-residue minus mask is empty). `L` = sum of the three frame lengths + 2 (`blast_query_info.c:140-161`).
  LOSAT `algorithm/blastx/query_setup.rs:protein_to_dna_masks` is the same formula.
* `Blast_GetSeqLocInfoVector` (`blast_aux.cpp:903-959`): per query, contexts 0..5, one `CSeqLocInfo` per interval with
  the frame (+1,+2,+3,-1,-2,-3), dropped if empty or equal to the whole query [0,L-1] (cannot happen for tblastx frames),
  no merging afterwards. Oracle (`qS8.o11`, L=324): plus1 `141-210`; plus2 `106-139,151-169`; plus3 `152-170`;
  minus1 `129-194,78-116` (descending DNA = ascending protein order), minus2 `149-172`, minus3 `130-213,40-72`.
* The masks reach only outfmt 0: `results.GetMaskedQueryRegions` -> `CDisplaySeqalign ... eLowerCase`
  (`blast_format.cpp:1547-1581`). Lowercase appears only on the Query row, only in the HSP's own frame: `qS8.o0` has
  `KGLSER...DFDTTaaaaaaaDYT`, `...FFSIVccccccccCI`, `*YNsssssssRLY`; with `-seg no` all upper case and 44 instead of 10 HSPs.
  Identities still count masked residues as matches (nomask).
* Subject masks: none without `-lcase_masking` (`seqinfosrc_seqvec.cpp:139-163`, `blast_input_aux.cpp:236`).
* LOSAT: engine masks per frame in `run_impl.rs` (`SegMasker`, `aa_seq` X, `aa_seq_nomask`, `seg_masks` end-exclusive);
  `TblastxQueryStats.seg_masks` carries them raw (not combined, not converted). Missing for the report:
  `combine_masks`, `protein_to_dna_masks` (both `pub` in `algorithm/blastx/query_setup.rs`, reusable), context order,
  and the minus-frame / empty-interval drops.

## 6. Oracle commands run (all under `scratch_E/`, run from that directory)

```
tblastx -query qA.fa -subject sA.fa -outfmt 0|6|7|8     # 8 kb / 7 kb windows (LC738874:233000-241000, LC738875:264000-271000), 206 HSPs
tblastx -query qT.fa -subject sT.fa -outfmt 0|6|8       # duplicated coding block, equal-score HSPs in several frames
tblastx -query qE.fa -subject sC.fa -seg no -outfmt 0|6|7|8|11   # biased frame +1 (SATQEK-rich); also default SEG
tblastx -query qS8.fa -subject sS8.fa [-seg no] -outfmt 0|8|11   # poly-Q insertion -> SEG masks, lowercase
tblastx -query qS8.fa -subject sDup.fa [-max_target_seqs 1|2] -outfmt 0|6|7   # tie order, hitlist
tblastx -query qS8.fa -subject sMulti.fa -outfmt 0|11            # 3 subjects, search space 220246
tblastx -query q2.fa|q4.fa|qN.fa|qStop.fa -subject sS8.fa [-seg no] -outfmt 0|6|7|11   # invalid frames / unsearched
tblastx -query qMulti.fa|qBigplus2.fa|qSmallplus2.fa -subject ... -outfmt 0|7        # batch behaviour, per-query footers
tblastx -query qPA.fa -subject sPA.fa -seg no -outfmt 8          # comparator-equal ties
```

## 7. What the port must do for range E (each item with its NCBI source)

1. Subject and HSP order, in all three formats: queries by index; subjects by `s_EvalueCompareHSPLists`
   (best_evalue with the 1e-180 rule, then raw score of `hsp_array[0]` = max raw score of the subject, then oid DESC;
   `blast_hits.c:3071-3107`, run at `blast_traceback.c:1771`); HSPs by e-value then `ScoreCompareHSPs` on the frame-relative
   offsets (`blast_hits.c:1329-1457`, `blast_seqalign.cpp:1577`). Keep comparator-equal HSPs in incoming order.
2. Apply hitlist_size / `-max_target_seqs` (`blast_hits.c:3243-3299`, `blast_traceback.c:1777`) and its stderr warning
   (`query_warnings.rs:few_matches_warning`); currently ignored by the TBLASTX engine.
3. Footer / ancillary data per query = first valid of the 6 frames (+1,+2,+3,-1,-2,-3): that frame's (masked-composition,
   check_ideal) block and its `eff_searchsp`, no Gapped table, no Gumbel columns (`blast_results.cpp:72-115`).
   Distinguish the three invalid situations of section 4, using the 10002-nt query batch as the unit
   (`blast_input.cpp:134-176`, `local_blast.cpp:177-225`, `blast_stat.c:2815-2823`), including NULL vs empty align set
   ("# 0 hits found" present or absent, `tabular.cpp:1277-1282`) and the stderr warning only for a fully invalid batch.
4. Carry per HSP: query frame, subject frame, `hsp->num` (sum_n shown when >1, e-value name sum_e), score, bit score of
   the HSP's own context, e-value (clamp to 0 below 1e-180 only at print time), num_ident (nomask), nt from/to/strand of
   both rows (formulas above). `Hit.num_positives` currently equals identities and `Hit.query_length` is 0: neither is
   printed by 0/6/7, but do not rely on them.
5. Query masks for outfmt 0: per frame SEG intervals, `combine_masks`, `protein_to_dna_masks` with `L = sum of frame
   lengths + 2`, frame order +1..-3, drop empty/whole-query, lowercase only on the Query row in the HSP's frame; no
   subject masks; `-seg no` -> no masks and Karlin blocks from the unmasked frames.

## 8. Open questions / UNSURE

* Comparator-equal HSP order (`qPA` case) depends on the engine's incoming order; confirm LOSAT reproduces rows 4-6 of
  `qPA.o8` (and the rest) before relying on the stable replay.
* `best_score` divergence (max raw score vs score of best-e-value HSP) was derived from the source, no oracle case built.
* `Blast_HSPListSortByEvalue` uses glibc `qsort`; ties are only reproducible if the incoming order is reproduced.
* LOSAT `max_target_seqs` is parsed but not applied; check whether a later stage outside `algorithm/tblastx` applies it
  (grep found nothing).
* The second SEG/ScoreBlk pass of the 6-argument `CBlastTracebackSearch` constructor is not used for `-subject`; not
  verified for any other front end (out of scope).
