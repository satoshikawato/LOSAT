# TBLASTN Stage D continuation — pinned NCBI call and input map

Branch: `feature/tlosan-tblastn-v0.2.0`; Stage C fixed-input gate: [STAGE_C_GATE_20260924.md](../tlosan_stage_c/STAGE_C_GATE_20260924.md). Pinned NCBI source: `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`. This is a **Stage D diagnostic checkpoint, not a Stage D pass**. The Rust TBLASTN path still accepts NCBI-recorded C-stage parameters as comparison inputs and the public CLI still reports unimplemented. No Stage E work is authorized by this evidence.

**2026-09-24 correction:** [CHECKPOINT_PARAMETERS_LINKING_20260924.md](CHECKPOINT_PARAMETERS_LINKING_20260924.md) supersedes the call-state hypothesis in this note. A new comparison-only probe measured positive `BlastSeqSrcGetTotLen` in all six saved local `-subject` runs, so the conditional `BLAST_OneSubjectUpdateParameters` call at `blast_engine.c:1434-1443` did **not** run. The cutoffs supplied to Stage C were initialized before subject iteration. The row and remaining-work section below have been corrected to use the observed initial setup state.

## Reproduction and authority

Run from the repository root:

```bash
python3 docs/evidence/tlosan_stage_d/run_ncbi_d_call_trace.py \
  /tmp/tlosan-stage-d-replay \
  docs/evidence/tlosan_stage_c/seg_hard_query_20260924 \
  docs/evidence/tlosan_stage_c/multi_query_20260924 \
  docs/evidence/tlosan_stage_c/run_20260923
```

A fresh replay to `/tmp/tlosan-d-replay-20260924-2` matched every saved file byte for byte (`diff -rq` exit 0). The focused Rust effective-length oracle test passed (1/1).

Two `multi_query` raw trace lines preserve NCBI's warning with its original trailing space. `git diff --cached --check` flags those evidence bytes; check source, tests, runner and prose with the raw `.trace` path excluded, and retain the trace byte-for-byte.

The runner pins `/home/kawato/micromamba/bin/tblastn` SHA-256 `e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0`, records each query/subject checksum and exact command, and checks comparison-only `LD_PRELOAD` results against the same unprobed NCBI command for stdout **and ordinary stderr bytes**. Full output, call/HSP traces, manifests and SHA-256 checksums are in [run_20260924](run_20260924). The two profiles are NCBI defaults (composition mode 2, sum statistics true, default SEG) and a control with `-comp_based_stats 0 -sum_stats false`; both use `-subject`, code 1, and one thread. Their differences are option effects, not LOSAT discrepancies. Outfmt 6 here is a numeric D diagnostic, not Stage E byte-parity evidence.

## NCBI function, input state and execution order

| Order | Pinned NCBI source | Required state and observed call |
| --- | --- | --- |
| 1 | `blast_setup.c:964-985`; conditional update at `blast_engine.c:1407,1434-1443` | Initial setup calculates effective lengths and hit/word/link cutoffs before subject iteration. The saved local sources have positive total length (362, 6377, 4694), so `BLAST_OneSubjectUpdateParameters` is skipped; [seqsrc_callstate_20260924](seqsrc_callstate_20260924) proves this for all six runs. `blast_setup.c:729-735,770-847` divides translated subject length by 3 inside effective-length calculation. `blast_engine.c:1446-1467` separately divides `stat_length` by 3. |
| 2 | `blast_parameters.c:774-815`, `blast_engine.c:870-905` | `do_sum_stats=true` creates link parameters. For translated gapped searches `-max_intron_length 0` selects `(DEFAULT_LONGEST_INTRON-2)/3`, observed `longest_intron=40`, so this is **uneven-gap linking**, not disabled linking. Preliminary `BLAST_LinkHsps` gets raw nucleotide `subject->length` (362 or 6377 here), before the preliminary E-value reap. With sum statistics disabled it calls `Blast_HSPListGetEvalues` with translated `stat_length` (120 or 2125). |
| 3 | `link_hsps.c:1765-1810`, `blast_hits.c:1811-1926` | Uneven-gap link resets HSP `num=1`, computes individual E-values with raw subject length `/3`, then links, score-sorts and recomputes best E-value. The probe records HSP raw score, context, frame, internal offsets, exact `double` E-value, `num`, and list order before/after. Observed gap decay is `0.10000000000000001`, and subject length 6377 is passed to link versus 2125 to its E-value call. `Blast_HSPListGetBitScores` uses `(score*Lambda-logK)/ln(2)` but the default Kappa path normalizes scores separately. |
| 4 | `blast_traceback.c:1481-1499`, `blast_kappa.c:2418-2478,2981-3125` | Composition mode 2 enters `Blast_RedoAlignmentCore_MT`, even on local `-subject`; the probe records one redo call per fixture. `blast_kappa.c:3290-3450,3520-3695` consumes ordered preliminary lists, converts incoming HSPs, redoes alignments, builds HSP lists and removes contained HSPs. |
| 5 | `blast_kappa.c:390-445,3695-3782` | After redo, `s_HitlistEvaluateAndPurge` calls link again with raw subject length (or direct E-value without sum statistics), then `Blast_HSPListReapByEvalue`, normalize score, identity, and local result writeback. The composition P-value E-value adjustment branch explicitly covers BLASTP/BLASTX, not TBLASTN. On the retained multiple-query fixture, query context 0 has `query_length=120`, `eff_searchsp=182004`, `length_adjustment=33`; context 1 has `70`, `88074`, `28`. Kappa's first link has `cutoff_small_gap=0` after preliminary `17` for that fixture; this reflects call timing. |
| 6 | `composition_adjustment/compo_heap.c:89-101,252-275,330-391`, `blast_traceback.c:1690-1777` | Kappa heap comparator orders E-value, score and subject index. The CLI local `-subject` path passes a non-NULL `seqSrc` (observed `redo` call state), so `blast_kappa.c:3719-3755` can use `BlastCompo_HeapWouldInsert` and heap insertion; the `!seqSrc` writeback branch is a different C++ API call state. Common post-pipes apply optional filters, sort results by E-value **only if** database total length is positive, then prune to hitlist size. Local `-subject` order must use its actual `seqSrc`/total-length state, not database output ordering. |

## Retained scores, rank and deletions

- Hard-SEG one-hit fixture: default NCBI output is raw **640**, bits **251**, E **1.73e-92**, while the control is raw **656**, bits **257**, E **6.52e-95**. The default preliminary link input is raw 656, and Kappa's redo/relink changes the final score; the full precision E-value before formatting is in the call trace. The first default link E-value is `7.2489135614429429e-95`, whereas direct control E-value is `6.5240222052986489e-95`. The parameter/branch difference is intentional.
- Multi-query fixture: default preliminary link takes **21** ordered HSPs. Kappa redo relinks **10** HSPs for query 0 and **9** for query 1; E-value reap leaves **7** and **6**. Query 0 deletion order at the reap list is `(-3,78,118,49,92,raw 771)`, `(-3,12,19,1430,1437,745)`, `(-1,33,46,962,975,677)`. Query 1 deletion order is `(-1,8,21,962,975,696)`, `(-3,53,69,49,65,644)`, `(3,6,21,1150,1165,592)`. Tuples are `(subject frame, query start/end, internal subject start/end, post-redo scaled raw score)`. The final local output has 13 rows in query/result order, recorded with raw, bit, E-value and coordinates. The control has 15 rows; its preliminary list/input cutoffs differ because its options differ.
- The 2026-09-23 retained fixture gives 11 default output rows and 11 control rows. Default trace records 24 link calls around one Kappa redo; the control has no link/redo calls. These calls and all HSP values are in the saved traces; a single hit-count equality is not accepted as parity.

## First unimplemented boundary and remaining work

The initial local parameter state is now computed in `stage_d_stats.rs` and compared with NCBI hit/word/link probes; a test-only Rust C→D path feeds those values into preliminary search, uneven-gap linking and preliminary reap for the multi-query and positive-link fixtures. The first unported D function boundary is translated-subject Kappa redo: [CHECKPOINT_KAPPA_WINDOWS_20260924.md](CHECKPOINT_KAPPA_WINDOWS_20260924.md) verifies its translated window ranges and amino-acid bytes, but `redo_alignment.rs` still explicitly rejects `subject_is_translated`, and there is no public integrated TBLASTN D pipeline. The saved post-Kappa HSP inputs also pass function-level Rust relink/reap comparisons, but those inputs come from NCBI. Implement the source-defined composition-mode-2 translation/range callbacks and redo first, then integrate containment, score normalization, identity, ranking and deletion order. The existing BLASTP path's `do_link_hsps=false` is not a TBLASTN rule. Every Rust port requires the NCBI file, line and snippet immediately above it.

Before closing D, compare Rust versus NCBI exact raw score, full-precision bit/E-value, rank and deletion order on code 1 fixtures. For non-default subject codes, isolate only the selected translation effect under `PD-TLOSAN-LOCAL-GENCODE-32`; use the comparison-only NCBI C++ API with `FindGeneticCode(32)` for code 32. Do not use NCBI `-db` statistics or headers as local `-subject` expected output. Keep the public CLI's explicit unimplemented error. Stage E outfmt 0/6/7 byte comparisons remain blocked until all of D matches.
