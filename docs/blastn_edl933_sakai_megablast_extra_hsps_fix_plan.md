# EDL933 vs Sakai Megablast Extra-HSP Fix Plan

Date: 2026-05-17

Goal: eliminate the current LOSATN default-megablast survivor-set delta for
`EDL933_vs_Sakai`, where LOSATN reports 5743 hits and NCBI BLAST+ reports 5718
hits. The target is bit-perfect NCBI BLAST+ parity, not a hit-count threshold.

NCBI BLAST is the sole source of truth. NCBI binaries may be used only as
diagnostic and validation oracles. LOSAT runtime code must not call, link,
embed, wrap, or delegate to NCBI BLAST.

## Current Diagnosis

Comparison by final outfmt coordinate key
`query, subject, qstart, qend, sstart, send`:

| Bucket | Count |
| --- | ---: |
| NCBI hits | 5718 |
| LOSAT hits | 5743 |
| Common coordinate keys | 5716 |
| NCBI-only coordinate keys | 2 |
| LOSAT-only coordinate keys | 27 |
| Net LOSAT surplus | +25 |
| Same-coordinate bit-score diffs | 0 |
| Same-coordinate E-value diffs | 0 |
| Same-coordinate pident/length/mismatch/gapopen diffs | 0 |

Interpretation:

- This is not a final statistic-formatting bug for shared HSPs.
- The defect is survivor-set divergence: LOSAT keeps different HSPs after
  preliminary greedy gapped extension, traceback, common-endpoint purging,
  ambiguity re-evaluation, final containment, and hit-list output.
- The most important event is a LOSAT-only merged HSP that replaces two NCBI
  HSPs.

NCBI-only:

```text
NC_002655.2 NC_002695.2 q=3551561-3559225 s=3484372-3492029 bit=14037
NC_002655.2 NC_002695.2 q=3559120-4038382 s=3491966-3971195 bit=8.839e+05
```

LOSAT-only:

```text
NC_002655.2 NC_002695.2 q=3551561-4038382 s=3484372-3971195 bit=8.975e+05
```

Trace from current LOSATN shows the merged HSP is produced during greedy
traceback from a low-score overlap seed:

```text
prelim=q3559119..3559225 s3491947..3492053 raw_score=106
traceback final=q3551560..4038382 s3484371..3971195 raw_score=485785
```

After that, common-endpoint purging naturally keeps the merged high-score HSP
and deletes the two NCBI-equivalent component HSPs. The later LOSAT-only short
perfect matches around `q=3559120..3559225` are also survivor-set fallout near
the same boundary.

## Primary NCBI References

Read these exact NCBI sections before changing Rust code. Every code patch must
include the relevant NCBI snippet comment with path and line numbers.

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3860-3922`
  - builds the temporary ungapped HSP and performs
    `BlastIntervalTreeContainsHSP` before preliminary gapped extension.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4012-4031`
  - megablast greedy seed selection and replacement with
    `greedy_query_seed_start` / `greedy_subject_seed_start`.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4058-4085`
  - accepted preliminary HSP initialization, save, and immediate interval-tree
    insertion.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2539-2569`
  - `s_BlastGreedyGapAlignStructFill`.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2762-2936`
  - `BLAST_GreedyGappedAlignment`, including traceback vs score-only behavior
    and best-start estimation.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:520-548`
  - preliminary `GetGappedScore` call followed by
    `Blast_HSPListPurgeHSPsWithCommonEndpoints(..., TRUE)` and score sort.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:633-692`
  - post-traceback common-endpoint trimming, ambiguity re-evaluation, second
    purge, score resort, and final containment purge.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2268-2535`
  - common-start/common-end comparators and
    `Blast_HSPListPurgeHSPsWithCommonEndpoints`.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_itree.c:809-847`
  - `s_HSPIsContained` and `MB_HSP_CLOSE`.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_itree.c:930-995`
  - `BlastIntervalTreeContainsHSP` traversal.

Primary LOSAT files:

- `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- `LOSAT/src/algorithm/blastn/alignment/greedy.rs`
- `LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs`
- `LOSAT/src/algorithm/blastn/interval_tree.rs`
- `LOSAT/src/algorithm/blastn/hsp.rs`
- `LOSAT/tests/compare_blastn_parity.py`

## Implementation Plan

### Phase 0: Freeze the Repro Case

Objective: make the failing case reproducible and impossible to confuse with
task-blastn or stale artifacts.

Tasks:

1. Add or update a manifest entry for `EDL933.Sakai.megablast` with exact
   options:
   `task=megablast`, reward `1`, penalty `-2`, gap open `0`, gap extend `0`,
   word size `28`, window size `0`, scan range `0`, DUST on, outfmt `6`,
   one thread.
2. Regenerate only this LOSAT output from the current binary.
3. Use the existing NCBI BLAST+ output as oracle unless a fresh NCBI run is
   explicitly requested for validation.
4. Save a small comparison artifact listing all 2 NCBI-only and 27 LOSAT-only
   coordinate keys.

Acceptance:

- The comparison still reports `common=5716`, `NCBI-only=2`, `LOSAT-only=27`,
  and no shared-coordinate statistic diffs.
- The manifest and command line make clear that this is default megablast, not
  `--task blastn`.

### Phase 1: Prove the First Divergent Stage

Objective: confirm whether LOSAT first diverges before preliminary greedy
extension, during score-only greedy start estimation, during final greedy
traceback, or during post-traceback survivor filtering.

Trace targets:

```text
Merged LOSAT-only target:
q=3551561-4038382 s=3484372-3971195

NCBI component targets:
q=3551561-3559225 s=3484372-3492029
q=3559120-4038382 s=3491966-3971195

Boundary short-HSP target:
q=3559120-3559201 s=3491990-3492071
```

Tasks:

1. Trace the overlap seed that produces the merged HSP:
   `q3559119..3559225 s3491947..3492053 raw_score=106`.
2. Trace the two NCBI component HSPs through:
   lookup, ungapped extension, preliminary interval-tree containment,
   preliminary greedy score-only extension, preliminary purge, final traceback,
   pass-1 endpoint purge, re-evaluation, pass-2 endpoint purge, final interval
   tree containment, and hit-list output.
3. Add focused trace output only if the current trace cannot show:
   preliminary HSP sort order, preliminary common-endpoint purge groups, and
   the exact tree HSP that should contain the overlap seed.
4. Do not patch algorithmic behavior until the first divergent stage is known.

Acceptance:

- The first stage where LOSAT differs from NCBI semantics is identified with
  NCBI source references.
- The explanation accounts for both the merged long HSP and the short perfect
  LOSAT-only boundary HSPs.

### Phase 2: Audit Preliminary Greedy HSP Generation

Objective: make the preliminary score-only greedy HSP list match NCBI before
traceback, because post-traceback purge is currently behaving consistently with
the HSP list LOSAT gives it.

Tasks:

1. Verify `greedy_gapped_alignment_score_only` in
   `LOSAT/src/algorithm/blastn/alignment/greedy.rs` against
   `BLAST_GreedyGappedAlignment` in `blast_gapalign.c:2762-2936`.
2. Confirm all of these NCBI details:
   - compressed subject pointer and remainder handling
   - grow-only greedy memory `max_dist`
   - right and left extension order
   - score conversion when `gap_open == 0 && gap_extend == 0`
   - `fwd_start_point` and `rev_start_point` interpretation
   - score-only best-start estimation choosing the middle of the longer exact
     match run inside the optimal alignment box
3. Verify the returned `greedy_query_seed_start` and
   `greedy_subject_seed_start` are the values later used for final traceback,
   not the original ungapped midpoint.
4. If LOSAT's score-only greedy start differs from NCBI, patch only that logic.

Acceptance:

- For the overlap seed, LOSAT produces the same preliminary range, score, and
  greedy seed start as NCBI.
- If NCBI would produce the overlap preliminary HSP too, the next phase must
  explain why NCBI does not let it merge across the component boundary.

### Phase 3: Audit Preliminary Containment and Purge Timing

Objective: ensure low-score overlap seeds are skipped or purged at the same
time as NCBI.

Tasks:

1. Verify the interval-tree precheck in `run.rs` uses the ungapped HSP range,
   not a preliminary gapped range, matching `blast_gapalign.c:3908-3919`.
2. Verify the interval tree is updated immediately after every accepted
   preliminary gapped HSP, matching `blast_gapalign.c:4058-4085`.
3. Verify LOSAT preliminary purge before traceback matches
   `blast_engine.c:520-548`:
   `Blast_HSPListPurgeHSPsWithCommonEndpoints(..., TRUE)` followed by score
   sorting.
4. Re-check whether `purge_prelim_hits_with_common_endpoints` exactly matches
   NCBI's common-start/common-end comparator order from `blast_hits.c:2268-2535`.
5. Confirm no LOSAT optimization changes qsort tie behavior where equal-score
   or near-equal endpoint groups affect the keeper.

Acceptance:

- The preliminary HSP list entering traceback contains the same component HSPs
  and does not contain an NCBI-invalid merged precursor.
- If the list still differs, the retained extra preliminary HSP has a documented
  NCBI rule that should have removed it.

### Phase 4: Audit Final Greedy Traceback Start Reuse

Objective: prevent score-only overlap seeds from causing an NCBI-invalid final
traceback span.

Tasks:

1. Verify final traceback uses the HSP start offsets saved by
   `Blast_HSPInit`, matching `blast_gapalign.c:4058-4076` and
   `blast_traceback.c:350-613`.
2. Verify LOSAT does not recompute or shift the start from a different seed
   after preliminary purge/sort.
3. Compare the final traceback for the overlap seed against NCBI's
   `do_traceback=TRUE` greedy path in `blast_gapalign.c:2762-2936`.
4. If the score-only HSP should survive but final traceback boundaries differ,
   patch the `greedy_gapped_alignment_with_traceback` path.

Acceptance:

- The overlap seed either no longer reaches traceback, or its final traceback
  matches NCBI and no longer forms `q3551561..4038382 / s3484372..3971195`.

### Phase 5: Verify Post-Traceback Filtering Is Not Masking the Real Bug

Objective: avoid papering over a preliminary/traceback defect with a late
filter that NCBI does not have.

Tasks:

1. Keep `purge_hsps_with_common_endpoints_ex` behavior tied to
   `blast_hits.c:2455-2535`.
2. Keep final interval-tree containment tied to `blast_itree.c:809-847` and
   `blast_itree.c:930-995`.
3. Do not add a custom "split merged HSP" or "remove short boundary HSP"
   heuristic. If such behavior is not in NCBI, it is not allowed.
4. After fixing the first divergent stage, re-run the full EDL933/Sakai
   comparison and inspect any remaining LOSAT-only rows.

Acceptance:

- The +25 net surplus is reduced by implementing NCBI timing/order, not by a
  LOSAT-specific final filter.
- Same-coordinate statistic diffs remain zero.

## Validation Plan

Run validation only after the patch for the identified NCBI divergence is in
place.

Minimum commands:

```bash
cd LOSAT
cargo build --release
python tests/compare_blastn_parity.py \
  --case EDL933_vs_Sakai \
  --ncbi tests/blast_out/EDL933.Sakai.blastn.megablast.out \
  --losat tests/losat_out/EDL933.Sakai.losatn.megablast.out \
  --limit 30
```

Expected final acceptance for this case:

- `ncbi_hits == losat_hits == 5718`
- `ncbi_only == 0`
- `losat_only == 0`
- no same-coordinate statistic diffs
- output row order matches NCBI BLAST+ outfmt 6

Regression checks after EDL933/Sakai is fixed:

- `Sakai_vs_MG1655` default megablast
- task-blastn manifest cases from
  `docs/losatn_blastn_bitwise_parity_remediation_plan.md`
- any existing unit tests touching:
  `greedy.rs`, `purge_endpoints.rs`, `interval_tree.rs`, and BLASTN output
  ordering

## Patch Rules

- Every Rust behavior change must include NCBI C/C++ reference comments with
  file path and line numbers immediately above the ported logic.
- Do not introduce behavior that has no NCBI equivalent.
- Do not call NCBI BLAST from LOSAT runtime, build scripts, unsupported
  feature paths, or fallback paths.
- Keep changes scoped to the first proven divergent stage.
- Do not treat hit-count improvement as success unless coordinate rows,
  values, and order are bit-perfect.
