# PS notes: BLASTP scoring, statistics and composition-based statistics

Rows: `PS.tsv` (65 rows). Scratch and evidence: `scratch_PS/` (inputs `q1.faa s5.faa qa.faa sa.faa qh.faa sh.faa`, oracle outputs `n_*`, LOSAT outputs `l_*`, scripts `cmp_matrices.py`, `cmp_stats.py`, `fl.py`, `addrows.py`). Oracle: NCBI 2.17.0 `blastp`; LOSAT: `LOSAT-before` (base commit). LOSAT rejects everything except BLOSUM62 11/1, comp_based_stats 2, so no non-default path of LOSAT could be run: "faithful" rows for matrix-generic code rest on reading the Rust and the C side by side, not on execution.

## Call path followed (NCBI, in execution order)

1. Arguments (`blast_args.cpp`): `CGenericSearchArgs::ExtractAlgorithmOptions` (255-300: `-evalue`; `BLAST_GetProteinGapExistenceExtendParams` when `-matrix` is given, so `-matrix X` alone resets both gap costs to the BEST pair of X; explicit `-gapopen`/`-gapextend` win individually), `CMatrixNameArg` (`SetMatrixName` keeps the string as typed), `CCompositionBasedStatsArgs` + `s_SetCompositionBasedStats` (795-886; first character selects the mode, 2nd character `u` the unified P; the arg is an unconstrained string), `CGappedArgs` (`-ungapped`), `-use_sw_tback` (flag, only sets `eTbackExt = eSmithWatermanTbck`; `ePrelimGapExt` stays score-only DP).
2. Validation: `BLAST_ValidateOptions` -> `BlastScoringOptionsValidate` (gapped: `Blast_KarlinBlkGappedLoadFromTables(NULL,...)`, `std_matrix_only = FALSE` for blastp so IDENTITY is allowed; ungapped: NO matrix lookup), `BlastHitSavingOptionsValidate` (evalue <= 0), `s_BlastExtensionScoringOptionsValidate` (CBS needs gapped), IDENTITY word size > 5.
3. `CSetupFactory::CreateScoreBlock` (setup_factory.cpp:143-156: IDENTITY resets CBS to 0 with a per-query warning) -> `BLAST_MainSetUp` -> `BlastSetup_ScoreBlkInit` -> `Blast_ScoreBlkMatrixInit` (name upper-cased into `sbp->name`; `options->matrix` keeps the typed case) -> `Blast_ScoreBlkMatrixFill` -> `BlastScoreBlkProteinMatrixLoad` (`NCBISM_GetStandardMatrix`) -> `BlastScoreBlkMaxScoreSet` -> `Blast_ScoreBlkKbpUngappedCalc` (+ `Blast_ScoreBlkKbpIdealCalc`; per-query kbp_std from the query composition; invalid contexts get a warning) -> `Blast_ScoreBlkKbpGappedCalc` (gapped only: `Blast_GumbelBlkCalc`, `Blast_KarlinBlkGappedCalc` per context; for `-ungapped` `sbp->gbp` is freed and `kbp_gap` stays NULL).
4. `BLAST_CalcEffLengths` (`BLAST_GetAlphaBeta`, `BLAST_ComputeLengthAdjustment`; gapped uses kbp_gap_std, ungapped uses sbp->kbp with row-0 alpha/beta).
5. `LookupTableWrapInit_MT` (plain `BlastAaLookup*` or compressed alphabet: `SCompressedAlphabetNew`, `RPSfindUngappedLambda`, `_PSIMatrixFrequencyRatiosNew`).
6. `BLAST_GapAlignSetUp`: `BlastScoringParametersNew`, `BlastExtensionParametersNew`, `BlastHitSavingParametersNew` (Spouge cutoff with `cbs_stretch` 5 for modes 2/3, or `BLAST_Cutoffs` when ungapped), `BlastInitialWordParametersNew/Update` (x_dropoff_init from the per-query ungapped lambda, gap trigger; ungapped branch with CUTOFF_E_BLASTP 1e-300).
7. Preliminary search (`aa_ungapped.c` word finders; `BLAST_GetGappedScore` with optional `s_ChainingAlignment`, restricted alignment for evalue <= 10; `BLAST_GetUngappedHSPList` + `Blast_HSPListReevaluateUngapped` + ungapped E-values/bit scores for `-ungapped`).
8. Traceback dispatch `BLAST_ComputeTraceback_MT` (blast_traceback.c:1481-1508): `compositionBasedStats > 0 || eTbackExt == eSmithWatermanTbck` -> `Blast_RedoAlignmentCore_MT` (`Blast_RedoOneMatch` or `Blast_RedoOneMatchSmithWaterman`, `Blast_AdjustScores` with `Blast_ChooseMatrixAdjustRule`, `Blast_CompositionMatrixAdj` or `Blast_CompositionBasedStats`), else `Blast_TracebackFromHSPList`; then `s_FilterBlastResults` (`-max_hsps`), `Blast_HSPResultsSortByEvalue`, `s_BlastPruneExtraHits`.
9. Report: `CBlastFormat` footer (Lambda/K/H of kbp_std per query, gapped block and a/alpha/sigma from the Gumbel block, `Matrix: <typed name>`), `showalign.cpp:3595-3604` Method label.

## Data the port needs (where, size)

| Data | NCBI location | Size | LOSAT today |
|---|---|---|---|
| Packed matrices | `src/util/tables/sm_*.c` | 9 x 625 | 8 of 9 present and identical (IDENTITY missing) |
| Karlin tables (Lambda,K,H,alpha,beta) | `blast_stat.c:183-587` | 98 rows | present for 6 matrices; PAM30 (4 rows missing, one alpha wrong), PAM70 (2 rows missing, ungapped beta wrong), IDENTITY missing |
| Gumbel columns C, alpha_v, sigma | same arrays, columns 9-11 | 98 x 3 numbers | 2 rows only (BLOSUM62 11/1, BLOSUM45 14/2) |
| Frequency ratios 28x28 | `matrix_freq_ratios.c` | 7 more x 784 doubles (~1400 lines) | BLOSUM62 only |
| Joint probabilities + background | `composition_adjustment/matrix_frequency_data.c` | 7 more x (400+20) doubles | BLOSUM62 only |
| Unified P table (565 bins) and BLOS62 scores | `unified_pvalues.c` | already ported | ported |
| Smith-Waterman | `smith_waterman.c` (615), `redo_alignment.c:1309-1555`, `blast_kappa.c:843-884,1812-1877` | ~600 C lines | not ported (explicit bail in redo_alignment.rs:2486) |

Gap pairs accepted per matrix (same arrays): BLOSUM45 13, BLOSUM50 15, BLOSUM62 11, BLOSUM80 9, BLOSUM90 7, PAM250 15, PAM30 10, PAM70 8, IDENTITY 1 (plus the ungapped row, which is not a valid user pair). LOSAT's TBLASTN module (`algorithm/tblastn/scoring.rs`) already has these pair lists for validation.

## What already exists in LOSAT (shared code)

Generic and verified in this pass: matrix tables for 8 matrices, `protein_score`, score-frequency/lambda/H/K code, length adjustment, Spouge StoE/EtoS, gapped DP (matrix as parameter, BLOSUM62 monomorph), `BlastGetStartForGappedAlignment`, positives/identities, display matrix (pinned for 8 matrices), kappa machinery (`Blast_RedoOneMatch`, heap, windows), RE-based adjustment, mode-1/mode-3 arms, unified-P formulas, generic ungapped extension, generic lookup builder `build_ncbi_lookup_for_profile` (dead code), generic ideal-lambda `ideal_karlin_params_for_matrix`, `matrix_score_bounds`.
Hard-wired to BLOSUM62 (the places to touch): `build_ncbi_lookup` / `prepare_blosum62_lookup_query_for_word_size` / `build_blosum62_compressed_lookup`, `extend_*_blosum62` calls (blast_engine.rs:4689, 4826), `compute_blosum62_ideal_karlin_params`, `BlastCompositionWorkspace::new_blosum62`, `build_matrix_info`, `blast_choose_matrix_adjust_rule` (BLOSUM62_BG), `lookup_protein_gumbel_params`, the chaining test, `validate_requested_blastp_support`.

## Surprises (things the session must know)

1. Default-option defect (row 47): outfmt 0 prints `Method: Compositional matrix adjust.` for every HSP; NCBI prints `Composition-based stats.` for alignments where `Blast_ChooseMatrixAdjustRule` returned `eCompoScaleOldMatrix` (37 of 125 HSPs on qa/sa). Known since the TBLASTX/segmask investigation; cheap to fix (carry `matrix_adjust_rule`; TBLASTN already has the formatter branch).
2. Default-option stderr/footer defect (row 9): a query with no valid ungapped Karlin block (all X) gets an NCBI warning and `-1.00` footers; LOSAT prints neither.
3. `-evalue` >= about 9e4 (row 36): LOSAT prints score-0 HSPs that NCBI does not; cause not found (UNSURE). `-task blastp-short` (default evalue 20000, 5x stretch) lives close to this regime.
4. `-comp_based_stats Xu` (unified P, row 41): NCBI 2.17.0 gives different output on repeated identical runs for single-query inputs (8 different results in 12 runs). Needs a policy decision: reject, or exception.
5. `-matrix IDENTITYX -ungapped` crashes NCBI (`free(): invalid size`, exit 134); other unknown names give `BLAST engine error: Error: Unknown error code -1` (exit 3) when `-ungapped` skips the table lookup (row 63).
6. The matrix name is kept as typed in the footer (`Matrix: blosum62`) and in the chaining test (`strcmp(..., "BLOSUM62")`, case-sensitive), but upper-cased everywhere else. LOSAT's enum loses the typed string (rows 54, 58).
7. NCBI uses `kFixedReBlosum62 = 0.44` as the target relative entropy for EVERY matrix in `Blast_CompositionMatrixAdj`; LOSAT reproduces that (row 46).
8. `-comp_based_stats` takes any string; only the first two characters matter (row 37). Two spellings (`2U`, `Du`) that look odd are valid and select unified P.
9. `-use_sw_tback` also redirects `-comp_based_stats 0` to the kappa code (dispatch at blast_traceback.c:1486); `eSmithWatermanTbckFull` / `eSmithWatermanScoreOnly` are not reachable from the CLI.
10. For ungapped blastp `do_sum_stats` stays FALSE (no `Expect(n)`), so only the plain `BLAST_Cutoffs` / `KarlinStoE_simple` path is needed; the footer has no Gapped block and no `Gap Penalties` line.
11. `-max_hsps` is faithful (outfmt 6/7 identical for N = 1,2,3,5,100); outfmt 0 differs only by the Method label.
12. LOSAT's `lookup_protein_params` silently falls back to the first gapped row for an absent pair, `lookup_protein_params_gapped` hard-codes 11/1 (used by TBLASTN-era code): do not reuse them for a general pair.

## Suggested order for the porting session

1. Fix tables (PAM30/PAM70), add IDENTITY, add the three Gumbel columns (data, M).
2. Wire the selected matrix through lookup, ideal lambda, per-query Karlin, ungapped extension (all S); lift the gate for BLOSUM62 gap pairs, then the other matrices with CBS 0 (needs Blast_TracebackFromHSPList, L) or CBS 2 (needs matrix_frequency_data + matrix_freq_ratios, M + L).
3. Method label (S), invalid-context warning (S), typed matrix string (S).
4. Modes 1 and 3: enable and compare with the oracle (S each).
5. `-ungapped` (M+M+S), blastp-short/blastp-fast (depend on the above), Smith-Waterman (L), compressed alphabets for other matrices (M).

## Decisions requested

- Unified P (`Xu`): reject or approved exception (non-reproducible NCBI output).
- `IDENTITYX` with `-ungapped`: NCBI crash; proposal reject.
- Whether to keep the typed matrix string through the whole pipeline (needed for footer and the chaining strcmp).
- Root-cause of the score-0 HSPs (row 36) before porting new -evalue defaults (blastp-short).
