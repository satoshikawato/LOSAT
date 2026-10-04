# RP notes: BLASTP report (outfmt 0, 6, 7) with non-default options and inputs

58 rows in `RP.tsv` (22 faithful, 18 divergent, 14 unported, 2 n/a, 2 exception, 0 rejected). Scratch in `scratch_RP/`
(`<name>.nc.*` = NCBI blastp 2.17.0, `<name>.lo.*` = LOSAT-before; `run2.sh NAME args...` runs both).

## Call path followed (NCBI, bl2seq `blastp -query Q -subject S`)

1. `blast_format.cpp`: `PrintProlog` 348-441 (bl2seq with a user subject is "DbScan", so the full prolog is printed:
   version, Altschul 1997 reference, composition-based-statistics reference, `PrintDbReport(top)` "Database: User
   specified sequence set (Input: <path as typed>)."), then per query `PrintOneResultSet` 1411-1590:
   `AcknowledgeBlastQuery` (Query=, Length=), no-hit message or `x_DisplayDeflines` (1-line table) + `CDisplaySeqalign`
   (alignments) + `x_PrintOneQueryFooter` 445-479 (`PrintKAParameters`, search space); after the last query
   `PrintEpilog` 2204-2297 (database block, `Matrix:`, `Gap Penalties` if gapped, `Neighboring words threshold`
   as a double, `Window for multiple hits`).
2. Table: `CShowBlastDefline::x_InitDeflineTable` 1051-1160 groups the per-HSP Seq-align list by subject and calls
   `GetSeqAlignSetCalcParams` (`align_format_util.cpp` 4247-4310: highest bit score HSP, strict `>`), then
   `x_DisplayDefline` 753-1002 (title cut at 68 with `...`). Note the NCBI bug at 1131-1139 (last row: the Score width
   is assigned from the length of the total-score string).
3. Alignments: `x_ShowAlnvecInfo` 3613 -> `x_PrintDefLine` (heading from `CDeflineGenerator::GenerateDefline`,
   `create_defline.cpp` 3928-4106) -> `x_DisplaySingleAlignParams` 3985 -> `x_DisplayAlignInfo` 3570 (Method text from
   `comp_adjustment_method`) -> `s_DisplayIdentityInfo` 304 -> `x_DisplayRowData` 1435, `x_OutputSeq` 2485 (lowercase of
   masked query residues). The Seq-align set comes from `blast_seqalign.cpp` (`s_BlastHSP2SeqAlign` drops score-0 HSPs,
   line 672; one flat Seq-align per HSP, sorted by `Blast_HSPListSortByEvalue`).
4. outfmt 6/7: `x_PrintTabularReport` 759-815, `CBlastTabularInfo::PrintHeader` 1266-1285 (`# N hits found` = number of
   Seq-aligns = HSPs; omitted together with `# Fields` when the search was not run).
5. Option inputs of the report: `CFormattingArgs::ExtractAlgorithmOptions` (`blast_args.cpp` 2885-2980): outfmt <= 4
   shows 500 descriptions and 250 alignments unless `-max_target_seqs`/`-num_*` are given.

## What I ran

`blastp` (NCBI) and `LOSAT-before` on: 7 pairs of the `tests/fasta` proteomes x `-evalue 10/1000` and `-seg yes` x outfmt 0/7
(`sweep/`, 42 runs: all outfmt 7 identical; outfmt 0 differs only by the Method text and, with `-seg yes`, the lowercase
Query residues); a 3-query x 1,176-subject run at `-evalue 1e8` (`manyA/B/C`); 78 protein subject titles
(`titles/title_sweep_results.txt`); query-title variants; fully masked / all-X queries (`lc*`, `tq*`); 60 kb and 90 kb
self/tandem hits (`huge0`, `tot`); the NCBI option matrix (`nc/*.out`: 8 matrices, `-comp_based_stats 0/1/3`, `-ungapped`,
tasks, gaps, word sizes 2-7, thresholds) to read what the report prints for each value.

## Findings that matter most

1. **Method text (high)**: LOSAT always prints "Compositional matrix adjust."; NCBI prints "Composition-based stats." for
   alignments whose `matrix_adjust_rule` is `eCompoScaleOldMatrix` (4% to 32% of HSPs in the 7 default-option pairs). LOSAT
   computes the rule but `Hit` has no field (`redone_hit_from_alignment_owned` in `blastp/kappa.rs` drops it; `comp_adjust_method`
   is `None` in `build_pairwise_hits`). Carry `comp_adjustment_method` (0/1/2 as `blast_kappa.c:331-342`) through `Hit` and
   the HSP copies. Needed for cbs 0 (no text), 1 and 3 as well.
2. **Protein subject titles (medium)**: BLASTP/TBLASTN/BLASTX still print the defline raw. 29 of 78 titles differ: trailing
   `.,;~ `, double spaces, ` ,`, `,,`, `( `, ` )`, `. [`/`, [` (protein only), HTML entities. `report/defline.rs` already
   ports the nucleotide pipeline; needs an `is_prot` flag for `x_CleanAndCompress` (`. [` -> ` [`, `, [` -> ` [`), a protein
   entry point, and either an `HtmlDecode` port (~2,100 entity names, `ncbistr.cpp:4223-4508`) or the BLASTN-style rejection
   (`check_shown_subject_title`). `x_AdjustProteinTitleSuffix` and `x_SetPrefix/x_SetSuffix` cannot change a FASTA subject
   (`m_Source` empty). Tab in a defline: NCBI cuts at the TAB (owned by IN).
3. **Zero-score HSPs (medium)**: the S08 "extra one-residue HSP" is an HSP with raw score 0 (`Score = 4.6 bits (0)`): NCBI
   drops it in `s_BlastHSP2SeqAlign` after the hit list was pruned; LOSAT prints it in outfmt 0/6/7 and counts it in
   `# N hits found`. Appears at `-evalue` 1000 and above (22 of 41,873 HSPs at `-max_target_seqs 1000 -evalue 1e8`).
4. **-seg yes lowercase (medium)**: `QueryFrame.seg_masks` already holds the intervals; the writer renders `query_nomask_sequence`
   in uppercase.
5. **Alignment count (medium)**: default `-num_alignments` is 250 (descriptions 500) when `-max_target_seqs` is not given;
   LOSAT writes an alignment block for every subject in the hit list (up to 500). `BlastnPairwiseReport` has the fields.
6. **Epilog threshold (medium)**: printed as a double (`11.5`, `19.3`, `21`, `20.25`); LOSAT casts to i32. Use `cpp_default_double`.
7. **Invalid queries (low)**: all-`X`, all-`*` or fully SEG-masked queries: NCBI warns on stderr, prints `-1.00` footers (all
   queries invalid) or empty footers (mixed batch) and, for outfmt 7 with all queries invalid, omits `# Fields` and
   `# N hits found`. LOSAT prints computed values; `query_warnings::invalid_query_warning` and
   `write_tblastn_unsearched_query_footer` already exist.
8. Smaller: Score-column width quirk of the last table row (row 15), first-HSP vs highest-bit-score rule (row 14, UNSURE,
   equal in 27,260 subjects), `O` residue replaced by X with a warning and different ungapped Karlin values (row 12).

## Data the port needs

- `src/algo/blast/core/blast_stat.c` 183-430: `blosum45/50/62/80/90`, `pam250/30/70_values` (14, 16, 12, 10, 8, 16, 11, 9 =
  96 rows) with 11 columns each: gap open, extend, INT2_MAX, Lambda, K, H, a, b, C, Alpha, Sigma. LOSAT `stats/tables.rs:214-314`
  has the first 7 columns and matches NCBI for BLOSUM45/50/62/80/90 and PAM250 (script comparison), but PAM30 lacks 4 rows
  (15/3, 14/2, 14/1, 13/3) and has `alpha` 0.34 instead of 0.35 at 10/1; PAM70 lacks 11/2 and 12/3 and has ungapped `beta` -0.9
  instead of -0.7. The Gumbel columns (C, Alpha, Sigma; ungapped `a_un`, `Alpha_un` = row 0 columns 6 and 9) exist only for
  BLOSUM62 11/1 and BLOSUM45 14/2 (`stats/spouge.rs:64`).
- Protein score matrices: `utils/matrix.rs` equals `util/tables/sm_*.c` for all 8 matrices (625 values each, checked) so the
  positives rule is already right for every matrix.
- `blast_stat.c:3696-3742` `Blast_GumbelBlkLoadFromTables` (formulas for b, Beta, Tau are in `spouge.rs`).
- `ncbistr.cpp:4223-4508` HTML entity table (if HtmlDecode is ported).
- Option defaults printed by the epilog come from the option range (PA): BLOSUM45 14/2 window 60 threshold 14, PAM30 9/1 window 15
  threshold 16 (blastp-short), word_size 5/6/7 thresholds 19.3/21/20.25.

## Surprises

- bl2seq protein subjects: the title shown is the whole defline including the ID (no `-parse_deflines`), so TPA/MAG prefix removal
  only matters when the ID itself starts with the prefix (`MULTISPECIES:` after an ID is untouched).
- `# N hits found` and the `Method` text are per HSP (flat Seq-align list), not per subject.
- NCBI invalid-query footers differ between "all queries invalid" (`-1.00` blocks, search space 0) and "mixed" (empty blocks, search
  space 0); the first comes from `CLocalBlast::Run` building ancillary data with `(-1,-1)` when the preliminary search was not run.
- `-max_target_seqs N` makes descriptions = alignments = N; only the default (no option) shows 250 alignments.
- NCBI `Warning: [blastp] Examining 5 or more matches is recommended` (hitlist < 5) is missing in LOSAT (owned by PV/PA).

## Decisions for the session

- Method text: port the per-HSP field (row 26) before any other BLASTP report work; it also gates cbs 0/1/3 reports (rows 27-29).
- Subject titles: port or reject entity decoding (row 21); I recommend reusing `defline.rs` with `is_prot` and the BLASTN rejection for
  entities, as the sweeps show the NCBI crash exception (PD-LOSAT-NCBI-DEFECTS 2) applies to the same `x_CleanAndCompress` code.
- Rows 4, 6, 24, 43 are owned by other ranges (IN/PV/PA); listed because the report bytes or stderr differ.
- Rows 32 (matrix positives), 39 and 40 (sort stability) need no report work now; row 40 matters only when `-task blastp-fast` is ported.
