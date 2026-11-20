# BLASTN Bit-Parity Fix Plan

Date: 2026-05-13

Goal: make LOSATN output bit-perfect with NCBI BLASTN for the supported
pure-Rust scope, without adding any runtime dependency on NCBI BLAST.

NCBI BLAST remains the only source of truth. NCBI binaries may be used only as
comparison oracles during diagnostics and validation. LOSAT runtime code must
not call NCBI BLAST, link to NCBI BLAST, or use NCBI BLAST as a fallback.

## Scope

Primary scope for this plan:

- `blastn --task blastn`
- `blastn` default megablast mode, only after task-blastn deltas are under
  control
- outfmt 6/7 parity for hit count, order, coordinates, bit score, e-value,
  identity, mismatches, gap opens, and alignment length
- single-query/single-subject comparison first, then multi-query/multi-subject

Out of scope:

- Delegating unsupported behavior to NCBI BLAST
- Refactoring for style without a confirmed parity bug
- Accepting hit-count thresholds as success. Hit-count deltas are diagnostics
  only; acceptance is bit-perfect output parity.

## Current Read

Do not compare LOSAT `--task blastn` against NCBI default `blastn`. NCBI default
is megablast unless `-task blastn` is specified. The repository currently has
both NCBI outputs:

- NCBI default/megablast:
  `LOSAT/tests/blast_out/NZ_CP006932.NZ_CP006932.blastn.out`
- NCBI task blastn:
  `LOSAT/tests/blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out`
- LOSAT task blastn:
  `LOSAT/tests/losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.out`

Current artifact comparison, using coordinate key
`query, subject, qstart, qend, sstart, send`:

| Case | LOSAT hits | NCBI hits | Common coords | LOSAT-only | NCBI-only | Common score diffs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `MelaMJNV.PemoMJNVA` task blastn | 2708 | 2729 | 2675 | 33 | 54 | 0 |
| `PemoMJNVA.PeseMJNV` task blastn | 2931 | 2940 | 2859 | 72 | 81 | 0 |
| `NZ_CP006932` self task blastn | 11479 | 12340 | 11144 | 335 | 1196 | 0 |

Interpretation:

- For the current task-blastn artifacts, same-coordinate HSPs have matching bit
  scores.
- The main remaining BLASTN problem is not raw scoring for shared HSPs. It is
  which HSPs survive seed, gapped extension, traceback, endpoint purge,
  re-evaluation, containment purge, and final hit-list pruning.
- The self comparison amplifies low-score and repetitive-sequence edge cases,
  so it should be used as the stress case after smaller virus-pair deltas are
  localized.

## Primary NCBI References

Read these before changing Rust code. Every parity code change must include the
relevant NCBI snippet comment with source path and line numbers.

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:154-193`
  - nucleotide initial word and gapped extension defaults
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:231-270`
  - `SetHitSavingOptionsDefaults` and `SetMBHitSavingOptionsDefaults`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_stat.c:674-684`
  - reward 2, penalty -3 Karlin tables and even-score rule
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:218-259`
  - ungapped x-drop and packed nucleotide score table setup
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:454-468`
  - gapped x-drop conversion from bit units to raw score units
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/na_ungapped.c`
  - nucleotide word scan, one-hit/two-hit behavior, ungapped extension
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3810-4105`
  - preliminary gapped alignment, interval tree containment, seed selection,
    HSP initialization
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3323-3390`
  - `BlastGetStartForGappedAlignmentNucl`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:350-613`
  - traceback loop, start offset selection, subject-range adjustment, final
    traceback, identity/length test, interval-tree insertion
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:633-692`
  - post-traceback endpoint purge, ambiguity re-evaluation, score re-sort, final
    containment purge
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:479-642`
  - `Blast_HSPReevaluateWithAmbiguitiesGapped`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1109-1132`
  - `Blast_HSPGetAdjustedOffsets`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1455`
  - score/e-value HSP comparators
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2268-2535`
  - common-start/common-endpoint sort keys and purge behavior
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_itree.c:146-166`
  - interval tree initialization
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_itree.c:797-847`
  - `s_HSPIsContained`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_itree.c:930-995`
  - `BlastIntervalTreeContainsHSP`

Primary Rust files:

- `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- `LOSAT/src/algorithm/blastn/blast_engine/mod.rs`
- `LOSAT/src/algorithm/blastn/extension.rs`
- `LOSAT/src/algorithm/blastn/alignment/gapped.rs`
- `LOSAT/src/algorithm/blastn/alignment/greedy.rs`
- `LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs`
- `LOSAT/src/algorithm/blastn/interval_tree.rs`
- `LOSAT/src/algorithm/blastn/hsp.rs`
- `LOSAT/src/algorithm/blastn/ncbi_cutoffs.rs`
- `LOSAT/src/stats/tables.rs`

## Phase 0: Freeze Valid Comparison Conditions

Objective: prevent false deltas caused by mismatched NCBI/LOSAT options.

Tasks:

1. Add a comparison manifest for BLASTN cases. Each row must record:
   `task`, `reward`, `penalty`, `gap_open`, `gap_extend`, `word_size`,
   `window_size`, `scan_range`, `dust`, `soft_masking`, `strand`, `evalue`,
   `max_target_seqs`, `max_hsps_per_subject`, `num_threads`, `outfmt`, query,
   and subject.
2. Split output directories or file names by NCBI task:
   - `*.megablast.*` for default megablast
   - `*.task_blastn.*` for `-task blastn`
3. Stop using `NZ_CP006932.NZ_CP006932.blastn.out` as the task-blastn oracle.
   Use `NZ_CP006932.NZ_CP006932.task_blastn.out`.
4. Add a small parser that reports:
   - coordinate-key matches
   - LOSAT-only coordinate keys
   - NCBI-only coordinate keys
   - same-coordinate score/e-value/identity differences
   - score-bin and alignment-length distributions
5. Select one first divergent HSP per failing case and store it in a manifest:
   `case_id`, `query_id`, `subject_id`, `qstart`, `qend`, `sstart`, `send`,
   NCBI score, LOSAT score if present, and expected first debug stage.

Acceptance:

- Every BLASTN parity discussion names the exact NCBI task and option set.
- Current task-blastn baseline is reproduced with no same-coordinate score
  differences.
- Each case has a single target HSP for trace-driven debugging.

## Phase 1: Add Stage-Level Differential Tracing

Objective: find the first pipeline stage where a specific HSP diverges.

Add tracing behind environment variables. The default runtime path must remain
silent and behaviorally unchanged.

Proposed variables:

- `LOSAT_TRACE_BLASTN_HSP="qstart,qend,sstart,send"`
- `LOSAT_TRACE_BLASTN_CONTEXT=<context_idx>`
- `LOSAT_TRACE_BLASTN_SUBJECT=<subject_id_or_index>`
- `LOSAT_TRACE_BLASTN_STAGE=<seed|ungapped|prelim|traceback|purge|all>`

Trace points:

1. Lookup hit retrieval:
   - subject offset
   - lookup index
   - query offsets returned from lookup
   - context index
   - query context offset
   - strand
2. Ungapped extension:
   - seed query/subject offsets
   - diagonal table decision
   - x-drop used
   - ungapped start/end
   - raw ungapped score
   - cutoff comparison
3. Preliminary gapped extension:
   - seed selected for DP or greedy extension
   - `x_drop_score_only = min(x_drop_gapped, ungapped_score)` for blastn
   - preliminary query/subject start/end
   - preliminary raw score
   - cutoff comparison
   - interval tree containment result before gapped extension
4. Preliminary common-endpoint purge:
   - sort key
   - common start/end group
   - kept HSP
   - removed or trimmed HSP
5. Traceback:
   - start selected by `BlastGetOffsetsForGappedAlignment` or
     `BlastGetStartForGappedAlignmentNucl`
   - subject `start_shift`
   - adjusted subject length
   - traceback x-drop
   - final query/subject start/end
   - raw score
   - edit script
   - identity/mismatch/gap counts
6. Post-traceback re-evaluation:
   - `extra_start`
   - re-evaluated score
   - e-value
   - deletion decision
   - identity/length test decision
7. Final endpoint purge and interval tree purge:
   - common endpoint groups
   - containment query/subject ranges
   - `min_diag_separation`
   - contained/not-contained result
   - HSP added to tree
8. Hit-list update and pruning:
   - HSP sort key
   - HSP-list sort key
   - `hitlist_size`
   - `max_hsps_per_subject`
   - final retained output order

Acceptance:

- For a selected NCBI-only HSP, LOSAT can report whether it disappeared at
  lookup, ungapped extension, preliminary gapped extension, traceback,
  re-evaluation, endpoint purge, interval containment, or final pruning.
- For a selected LOSAT-only HSP, LOSAT can report the first stage where NCBI
  would have rejected or altered it.
- Trace blocks include NCBI source comments immediately above the corresponding
  Rust logic when code is modified.

## Phase 2: Verify NCBI Option and Statistics Parity

Objective: remove hidden option/stat differences before algorithm changes.

Tasks:

1. Confirm task defaults from NCBI:
   - `blast_nucl_options.cpp:154-193` for word/gapped defaults
   - `blast_nucl_options.cpp:231-270` for hit-saving defaults
2. Verify LOSAT task config:
   - `task blastn`: reward 2, penalty -3, gap open 5, gap extend 2,
     `eDynProgScoreOnly`, `eDynProgTbck`, `min_diag_separation=50`
   - megablast: reward 1, penalty -2, gap open 0, gap extend 0,
     greedy score-only and greedy traceback, `min_diag_separation=6`
3. Verify Karlin parameter lookup:
   - reward 2, penalty -3, gap 5/2 must use `lambda=0.625, K=0.41`
   - reward 2, penalty -3, ungapped 0/0 must use `lambda=0.55, K=0.21`
   - even-score rule must be applied where NCBI applies it
4. Verify x-drop conversion:
   - ungapped x-drop uses ungapped lambda and `ceil`
   - gapped and final x-drop use smallest gapped lambda and C cast semantics
5. Verify effective search space:
   - context-specific `eff_searchsp`
   - plus and minus query contexts
   - length adjustment uses gapped parameters for final e-values

Acceptance:

- Debug output can print a compact NCBI/LOSAT option/stat block for each test
  case.
- All option/stat values are either proven identical to NCBI or have a tracked
  discrepancy with NCBI line references.

## Phase 3: Localize Seed and Ungapped Differences

Objective: prove whether missing/extra HSPs originate before gapped extension.

Tasks:

1. For one NCBI-only HSP from `MelaMJNV.PemoMJNVA`, trace all lookup hits within
   a small window around its query/subject coordinates.
2. Compare lookup behavior with NCBI `na_ungapped.c` and `blast_nalookup.c`:
   - lookup table type
   - direct/PV/two-stage path
   - query offset encoding
   - subject packed sequence scan step
   - mask-at-hash behavior
   - ambiguous-base handling
3. Verify diagonal table behavior:
   - one-hit mode for nucleotide default window size
   - diag-table offset advancement by subject length and window
   - diag reset threshold
   - hash vs array mode threshold
4. Verify ungapped extension:
   - exact vs inexact extension path
   - x-drop sign and threshold
   - boundaries at query context ends and subject ends
   - use of masked query for extension and unmasked query for identity where
     NCBI does so
5. Record whether the selected divergent HSP produces a matching ungapped HSP
   before gapped extension.

Acceptance:

- If selected HSPs already diverge in seed/ungapped stages, patch this phase
  before touching gapped extension.
- If selected HSPs match through ungapped stage, freeze seed/ungapped code and
  move to Phase 4.

## Phase 4: Fix Preliminary Gapped Extension Parity

Objective: make the preliminary gapped HSP list match NCBI before traceback.

Primary Rust target:

- `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- `LOSAT/src/algorithm/blastn/alignment/gapped.rs`

NCBI behavior to mirror:

- `blast_gapalign.c:3810-3832`: interval tree initialized with query length and
  subject length plus one.
- `blast_gapalign.c:3900-3919`: ungapped HSP is checked against prior gapped
  HSPs before new gapped extension.
- `blast_gapalign.c:4033-4046`: for non-greedy blastn, if the ungapped HSP is
  long enough, seed offsets advance by 3 before DP gapped alignment.
- `blast_gapalign.c:4058-4077`: only `gap_align->score >= cutoff` creates a
  new HSP.
- `blast_gapalign.c:4087-4088`: accepted preliminary gapped HSPs are added to
  the interval tree immediately.

Tasks:

1. Instrument preliminary gapped HSP list before endpoint purge:
   `(context, q_off, q_end, s_off, s_end, seed_q, seed_s, score)`.
2. Verify `x_drop_score_only` for blastn:
   - LOSAT currently uses `min(x_drop_gapped, ungapped_score)`.
   - Confirm against NCBI `blast_gapalign.c:2959-2963` and the actual caller.
3. Verify DP extension boundary behavior:
   - subject packed-vs-blastna buffers
   - 4-base boundary rounding
   - inclusive/exclusive stops
   - gap open/extend signs
4. Verify interval tree insertion timing:
   - add only accepted preliminary gapped HSPs
   - keep tree state across the subject's preliminary gapped loop
   - reset only where NCBI resets
5. Verify subject chunk merging:
   - overlap size
   - offsets added after per-chunk purge
   - duplicate removal across chunk overlap
   - order of merge and score sort

Acceptance:

- For selected divergent HSPs, preliminary gapped HSP presence/absence matches
  NCBI.
- Preliminary gapped coordinates and raw scores match NCBI for matching HSPs.

## Phase 5: Fix Traceback Start and Boundary Parity

Objective: make final traceback coordinates and raw scores match NCBI.

Primary Rust targets:

- `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- `LOSAT/src/algorithm/blastn/alignment/gapped.rs`
- `LOSAT/src/algorithm/blastn/alignment/greedy.rs`

NCBI behavior to mirror:

- `blast_traceback.c:436-460`: choose start with
  `BlastGetOffsetsForGappedAlignment` if both gapped starts are zero, otherwise
  use `BlastGetStartForGappedAlignmentNucl` for blastn.
- `blast_traceback.c:466-472`: apply `AdjustSubjectRange`.
- `blast_traceback.c:503-512`: run greedy or DP traceback over the adjusted
  subject range.
- `blast_traceback.c:583-600`: update HSP with traceback, run identity/length
  test for non-greedy, then adjust subject offset by `start_shift` and add to
  interval tree.

Tasks:

1. Add trace output for traceback start selection and `start_shift`.
2. Verify `BlastGetStartForGappedAlignmentNucl` against NCBI line-by-line:
   - how it scores candidate start points
   - how it handles gaps near the start
   - exact tie-break behavior
3. Verify `AdjustSubjectRange`:
   - shift amount
   - adjusted subject length
   - compensation when writing final subject offsets
4. Verify traceback DP:
   - recurrence equations
   - x-drop stop conditions
   - band/window limits
   - tie-break order among match, insertion, deletion states
   - edit script op semantics
5. Verify identity statistics after traceback:
   - use unmasked query where NCBI uses `query_nomask`
   - subject sequence must be canonical plus subject sequence for blastn
   - gap-open counting follows NCBI edit-script walk

Acceptance:

- Same-coordinate score differences remain zero.
- High-score coordinate differences in `NZ_CP006932` self are explained and
  reduced before focusing on low-score tail HSPs.

## Phase 6: Fix Common-Endpoint Purge and Re-evaluation

Objective: match NCBI decisions for HSPs sharing starts or ends.

Primary Rust target:

- `LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs`

NCBI behavior to mirror:

- `blast_hits.c:2268-2321`: `s_QueryOffsetCompareHSPs`
  - context ascending
  - query start ascending
  - subject start ascending
  - score descending
  - query end descending
  - subject end descending
- `blast_hits.c:2332-2387`: `s_QueryEndCompareHSPs`
  - context ascending
  - query end ascending
  - subject end ascending
  - score descending
  - query start descending
  - subject start descending
- `blast_hits.c:2392-2452`: `s_CutOffGapEditScript`
- `blast_hits.c:2455-2535`: `Blast_HSPListPurgeHSPsWithCommonEndpoints`
- `blast_traceback.c:637-668`: first purge trims, re-evaluation follows, second
  purge deletes for blastn.

Tasks:

1. Add tests for comparator tie-breaks from NCBI comments.
2. Add tests for `s_CutOffGapEditScript` on:
   - cut inside a substitution run
   - cut at a run boundary
   - cut through insertion/deletion runs
   - failed cut where NCBI leaves behavior unchanged
3. Verify first purge returns `extra_start` equivalent to NCBI, not just a new
   vector length.
4. Verify greedy traceback special case:
   - NCBI sets `extra_start = 0` for greedy traceback.
   - This affects megablast more than task blastn, but keep it in the same
     purge test suite.
5. Verify re-evaluation:
   - `Blast_HSPReevaluateWithAmbiguitiesGapped` score and deletion decision
   - e-value recomputation using context `eff_searchsp`
   - identity/length test immediately after re-evaluation
6. Remove any extra top-level duplicate purge that is not present in the NCBI
   traceback flow. If such code is retained for legacy output paths, prove it
   does not affect blastn task output.

Acceptance:

- Common-endpoint groups produce the same kept/trimmed/deleted HSP as NCBI.
- The low-score LOSAT-only and NCBI-only tail shrinks without changing shared
  HSP scores.

## Phase 7: Fix Interval Tree Containment Parity

Objective: make containment-based deletions match NCBI.

Primary Rust target:

- `LOSAT/src/algorithm/blastn/interval_tree.rs`

NCBI behavior to mirror:

- `blast_itree.c:146-166`: tree initialization.
- `blast_itree.c:797-847`: score, frame sign, start containment, end
  containment, and `MB_HSP_CLOSE` logic.
- `blast_itree.c:930-995`: query-indexed traversal.
- `blast_traceback.c:674-688`: final score-sort followed by tree reset and
  containment purge.

Tasks:

1. Add unit tests for `s_HSPIsContained` exact macro semantics:
   - score equal
   - score lower
   - score higher
   - plus and minus query contexts
   - subject frame sign
   - `min_diag_separation=0`
   - `min_diag_separation=50`
2. Verify context offset handling:
   - plus context offset
   - minus context offset
   - query-concatenated coordinate length equals NCBI query block length
3. Verify tree traversal:
   - middle comparison uses the same inclusive/exclusive behavior as NCBI
   - midpoint subject tree stores and traverses HSPs in the same way
   - leaf `leftptr` use for query context offset is equivalent to NCBI state
4. Add trace output showing:
   - candidate HSP
   - tree HSP that contained it
   - endpoint containment booleans
   - diagonal close booleans
5. Compare final containment removals for the self case:
   - LOSAT current debug shows only 187 removals from the final tree pass in
     one old timing artifact.
   - Confirm whether NCBI removes more, or whether earlier stages already
     differ.

Acceptance:

- For selected LOSAT-only HSPs, containment trace explains whether NCBI would
  delete them and why.
- Final containment pass removes the same HSPs as NCBI for the focused cases.

## Phase 8: Fix Query/Subject Coordinate Transforms

Objective: keep internal canonical coordinates separate from output
coordinates, matching NCBI.

Primary Rust targets:

- `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- `LOSAT/src/algorithm/blastn/hsp.rs`
- `LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs`

NCBI behavior to mirror:

- `blast_hits.c:1109-1132`: `Blast_HSPGetAdjustedOffsets`
- `blast_query_info.c:248-287`: context offset and total query length logic
- `blast_traceback.c:651-655`: re-evaluation query pointer uses the HSP context

Tasks:

1. Audit all uses of output coordinates in internal pruning paths.
2. Ensure score sorting, endpoint purge, containment purge, and re-evaluation
   use internal contexted offsets.
3. Apply output coordinate adjustment only when formatting or constructing the
   final report row.
4. For minus query context:
   - query output coordinates are flipped
   - subject output coordinates follow NCBI adjusted-offset behavior
   - internal subject offsets remain canonical plus-subject offsets
5. Add regression tests with:
   - plus query HSP
   - minus query HSP
   - same internal HSP rendered to output
   - endpoint purge on minus context

Acceptance:

- No pruning/filtering code depends on already-adjusted output coordinates.
- Minus-strand output matches NCBI while internal pruning remains stable.

## Phase 9: Final Hit-List Sorting and Pruning

Objective: match NCBI output order and final hit-list limits.

Primary Rust targets:

- `LOSAT/src/algorithm/blastn/hsp.rs`
- `LOSAT/src/algorithm/blastn/blast_engine/run.rs`

NCBI behavior to mirror:

- `blast_hits.c:1330-1455`: HSP score/e-value comparators.
- `blast_hits.c:3077-3106`: HSP-list e-value comparator.
- `blast_hits.c:3383-3400`: result sorting and purge.
- `blast_traceback.c:836-892`: max-HSP-per-subject and extra-hit pruning.

Tasks:

1. Verify HSP sort key:
   - score descending
   - subject offset ascending
   - subject end descending
   - query offset ascending
   - query end descending
2. Verify HSP-list sort key:
   - best e-value
   - best score
   - subject oid tie-break
3. Verify hit-list pruning:
   - `max_hsps_per_subject`
   - `hitlist_size`
   - `max_target_seqs`
   - order of sorting before pruning
4. Verify outfmt 7 comments:
   - hit count line
   - field order
   - NCBI formatting of e-value and bit score

Acceptance:

- Same retained HSP set as NCBI before formatting.
- Output row order matches NCBI.

## Phase 10: Regression and Acceptance Suite

Objective: prevent regressions while closing BLASTN parity.

Focused cases:

1. `MelaMJNV.PemoMJNVA` task blastn
   - small residual delta
   - good first case for endpoint/traceback issues
2. `PemoMJNVA.PeseMJNV` task blastn
   - small residual delta
   - includes minus-strand examples
3. `NZ_CP006932` self task blastn
   - stress case for repetitive low-score HSPs and containment
4. `EDL933.Sakai` megablast
   - run after task blastn is stable
5. `Sakai.MG1655` megablast
   - run after task blastn is stable

Validation sequence:

1. Run focused unit tests for the changed function only.
2. Run one selected comparison case with tracing disabled.
3. Diff coordinate-key sets and same-coordinate values.
4. Run all BLASTN comparison cases.
5. Run broader `cargo test` only after the active parity sweep has no known
   unresolved NCBI divergence in the touched stage.

Acceptance:

- For each focused case:
  - hit count matches
  - coordinate-key set matches
  - raw score/bit score/e-value matches
  - identity/mismatch/gap-open/alignment length matches
  - output row order matches
- Trace code remains gated and inactive by default.
- Every code change includes NCBI reference comments with file path and line
  numbers.

## Recommended Work Order

1. Phase 0: fix comparison discipline and choose one divergent HSP per case.
2. Phase 1: add trace gates for selected HSPs.
3. Phase 4: compare preliminary gapped HSP lists first, because current shared
   HSP scores already match and the hit-set delta likely starts at survival
   decisions.
4. Phase 5: fix traceback start/boundary if preliminary gapped HSPs match but
   final coordinates differ.
5. Phase 6: fix endpoint purge/re-evaluation for short low-score tail HSPs.
6. Phase 7: fix interval containment only after the input HSP list to the tree
   is proven to match NCBI.
7. Phase 9: fix final sorting/pruning after the retained HSP set matches.

Do not start by tuning thresholds or adding broad filters. If a LOSAT-only HSP
is wrong, identify the exact NCBI stage that removes it and port that behavior.

## Immediate Next Diagnostic

Use `NZ_CP006932` task blastn only with the correct oracle:

- LOSAT:
  `LOSAT/tests/losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.out`
- NCBI:
  `LOSAT/tests/blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out`

Select:

- one NCBI-only high-ish score HSP from the first 50 NCBI-only keys
- one LOSAT-only high-ish score HSP from the first 50 LOSAT-only keys
- one low-score perfect-match LOSAT-only HSP

For each, trace through:

1. lookup retrieval
2. ungapped extension
3. preliminary gapped extension
4. traceback
5. endpoint purge pass 1
6. re-evaluation
7. endpoint purge pass 2
8. final interval tree purge

The first stage that disagrees determines the first code patch.
