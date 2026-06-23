# BLASTN v0.1.0 Supported-Parity Completion Plan

Status: planned.

Goal: make BLASTN eligible for v0.1.0 supported scope by eliminating every
current NCBI BLAST+ difference in `tests/blastn_parity_manifest.tsv`, then
recording fresh release evidence from a clean release branch.

This plan is only for promoting BLASTN. If v0.1.0 ships with BLASTN marked
experimental, this plan is not a release blocker.

## Mandatory Rules

- NCBI BLAST C/C++ source is the only behavior authority.
- Do not guess a fix from LOSAT output alone.
- Do not add filters, thresholds, tie-breakers, or special cases unless the same
  behavior is present in the corresponding NCBI source path.
- Do not call NCBI BLAST executables from LOSAT runtime, build code, fallback
  paths, or unsupported-feature handling.
- Every Rust code change must include an NCBI source comment immediately above
  the changed logic, including file path, exact line numbers, and the relevant
  C/C++ snippet.
- If the corresponding NCBI source cannot be found, do not implement the
  behavior. Leave BLASTN experimental instead.

## Current Baseline

Baseline command, from `LOSAT/`:

```bash
python3 tests/compare_blastn_parity.py \
  --manifest tests/blastn_parity_manifest.tsv \
  --limit 3
```

Current artifact summary:

| Case | Mode | NCBI hits | LOSAT hits | Current difference |
| --- | --- | ---: | ---: | --- |
| `PesePMNV.MjPMNV.task_blastn` | task-blastn | 241 | 241 | exact |
| `PmeNMV.MjPMNV.task_blastn` | task-blastn | 208 | 208 | exact |
| `PmeNMV.PesePMNV.task_blastn` | task-blastn | 1431 | 1431 | exact |
| `PeseMJNV.PemoMJNVB.task_blastn` | task-blastn | 11668 | 11660 | 8 NCBI-only HSPs |
| `MelaMJNV.PemoMJNVA.task_blastn` | task-blastn | 2729 | 2724 | 5 NCBI-only HSPs |
| `PemoMJNVA.PeseMJNV.task_blastn` | task-blastn | 2940 | 2939 | 1 NCBI-only HSP |
| `MjeNMV.MelaMJNV.task_blastn` | task-blastn | 2668 | 2664 | 5 NCBI-only HSPs, 1 LOSAT-only HSP |
| `MjPMNV.MlPMNV.task_blastn` | task-blastn | 54402 | 54401 | 1 NCBI-only HSP |
| `NZ_CP006932.NZ_CP006932.task_blastn` | task-blastn | 12340 | 12340 | exact |
| `EDL933.Sakai.megablast` | megablast | 5718 | 5718 | exact |
| `Sakai.MG1655.megablast` | megablast | 6476 | 6476 | 2 pident diffs, 5 length/mismatch/gapopen diffs |

Observed shape:

- task-blastn shared-coordinate HSPs currently have zero bit-score, E-value,
  pident, length, mismatch, and gapopen differences.
- task-blastn remaining defects are final survivor-set differences: LOSAT loses
  a small number of NCBI HSPs and keeps one LOSAT-only HSP.
- megablast remaining defects keep the same hit count, coordinates, bit score,
  and E-value, but disagree in edit-script-derived statistics for a few shared
  HSPs.

Before any patch, refresh this baseline and save:

- LOSAT commit.
- Rust version.
- NCBI `blastn -version`.
- Exact commands.
- The full `compare_blastn_parity.py` output.
- Input and output file checksums.

Suggested artifact directory:

```bash
mkdir -p .tmp/blastn_supported_parity_$(date +%Y%m%d)
```

## NCBI Source Map

Read and cite exact line ranges from these files before changing Rust:

| Concern | Required NCBI source to inspect |
| --- | --- |
| Preliminary gapped HSP creation and interval-tree precheck | `c++/src/algo/blast/core/blast_gapalign.c`, especially the ungapped-HSP precheck, `Blast_HSPInit`, score-only extension, and accepted-HSP insertion paths |
| Final traceback update and post-traceback cleanup | `c++/src/algo/blast/core/blast_traceback.c`, especially traceback update, identity/length test, common-endpoint purge, ambiguity re-evaluation, score resort, and final containment purge |
| Common endpoint purge and edit-script trimming | `c++/src/algo/blast/core/blast_hits.c`, especially `s_QueryOffsetCompareHSPs`, `s_QueryEndCompareHSPs`, `s_CutOffGapEditScript`, and `Blast_HSPListPurgeHSPsWithCommonEndpoints` |
| Final HSP coordinate conversion and sorting | `c++/src/algo/blast/core/blast_hits.c`, especially `Blast_HSPGetAdjustedOffsets`, `ScoreCompareHSPs`, and HSP-list comparators |
| Interval containment semantics | `c++/src/algo/blast/core/blast_itree.c`, especially `s_HSPIsContained`, `BlastIntervalTreeContainsHSP`, and megablast close-HSP logic |
| Megablast greedy traceback and edit scripts | `c++/src/algo/blast/core/greedy_align.c` and the megablast call sites in `blast_gapalign.c` |
| Final score, bit score, and E-value assignment | `c++/src/algo/blast/core/blast_stat.c`, `blast_hits.c`, and traceback call sites that choose per-context effective search space |
| Nucleotide lookup and ungapped seed behavior, if a traced HSP diverges before gapped extension | `c++/src/algo/blast/core/na_ungapped.c` and nucleotide lookup files used by the task |

Use `/mnt/c/Users/genom/GitHub/ncbi-blast/` first. If a file is missing there,
try `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/`.

## Phase 0: Freeze Scope And Comparison Discipline

1. Confirm the manifest rows are the complete BLASTN supported-candidate scope.
2. Confirm every row names the exact NCBI task and option set:
   - task-blastn: reward `2`, penalty `-3`, gap open `5`, gap extend `2`,
     word size `11`, dust on, soft masking on, both strands.
   - megablast: reward `1`, penalty `-2`, gap open `0`, gap extend `0`,
     word size `28`, dust on, soft masking on, both strands.
3. Stop treating any stale default `blastn` output as a task-blastn oracle.
4. Extend the saved summary, if needed, so each case records:
   - hit count.
   - common coordinate keys.
   - NCBI-only and LOSAT-only keys.
   - shared-coordinate bit score, E-value, pident, length, mismatch, and
     gapopen differences.
   - order and formatting differences for the exact release output mode.

Exit criteria:

- The current failing cases and exact defect classes are reproducible from the
  manifest.
- The artifact directory contains all commands and versions needed to reproduce
  the baseline.

## Phase 1: Diagnose The Task-BLASTN Survivor-Set Delta

Start with one high-impact NCBI-only HSP and the one LOSAT-only HSP:

- `MjeNMV.MelaMJNV.task_blastn` NCBI-only:
  `q=147883-148031 s=257891-257755`, length `149`, bit score `46.4`.
- `MjeNMV.MelaMJNV.task_blastn` LOSAT-only:
  `q=144023-144052 s=270709-270737`, length `30`, bit score `33.7`.

Trace each target through the BLASTN pipeline:

1. lookup seed.
2. ungapped extension.
3. preliminary gapped score-only extension.
4. preliminary HSP creation and immediate interval-tree precheck.
5. preliminary HSP insertion.
6. preliminary common-endpoint purge.
7. final traceback.
8. identity/length test.
9. ambiguity re-evaluation.
10. post-traceback common-endpoint purge.
11. score resort.
12. final interval-tree containment purge.
13. hitlist pruning.
14. output coordinate conversion.

Use existing diagnostics first:

```bash
LOSAT_TRACE_BLASTN_HSP="qstart,qend,sstart,send"
LOSAT_TRACE_BLASTN_CONTEXT=<context_idx>
LOSAT_TRACE_BLASTN_SUBJECT=<subject_id_or_index>
LOSAT_TRACE_BLASTN_STAGE=all
```

Add new diagnostics only if an existing trace cannot identify the first
divergent stage. Diagnostic code must still cite the NCBI source stage it is
mirroring.

Exit criteria:

- For the selected NCBI-only HSP, the first LOSAT stage where it disappears is
  known.
- For the selected LOSAT-only HSP, the first stage where NCBI would reject or
  purge it is known.
- The explanation identifies one NCBI function or comparator to port or fix.

## Phase 2: Patch The First Task-BLASTN Offender

Patch only the first proven offender from Phase 1. Likely areas include, but
must not be assumed:

- preliminary HSP interval-tree precheck or insertion order.
- score-only gapped extension acceptance timing.
- final traceback start/end reconstruction.
- common-endpoint trimming or deletion.
- ambiguity re-evaluation deletion decisions.
- final containment purge.
- final score/e-value ordering or hitlist pruning.

For every Rust edit:

1. Read the exact NCBI function.
2. Copy the relevant C/C++ snippet into a nearby comment.
3. Include source path and line numbers.
4. Preserve NCBI timing and order even if a Rust-local simplification looks
   equivalent.
5. Add or update a focused unit test only when the behavior is a ported helper
   with a small deterministic boundary.

Verification after the patch:

```bash
python3 tests/compare_blastn_parity.py \
  --manifest tests/blastn_parity_manifest.tsv \
  --limit 3
```

Exit criteria:

- The targeted HSP delta is gone or the trace proves a later NCBI source path is
  now the first offender.
- No previously exact BLASTN manifest case regresses.
- Shared-coordinate score and E-value differences remain zero.

## Phase 3: Sweep Remaining Task-BLASTN Cases

After the first offender is fixed, rerun the full manifest and group any
remaining task-blastn deltas by first divergent stage. Work in this order:

1. high-impact NCBI-only HSPs by `length * bitscore`.
2. LOSAT-only HSPs.
3. lower-score NCBI-only tail HSPs.

Do not tune thresholds globally to remove the tail. Each retained or removed HSP
must follow the same NCBI stage decision.

Exit criteria for task-blastn:

- All 9 task-blastn manifest rows have equal NCBI and LOSAT hit counts.
- `NCBI-only = 0`.
- `LOSAT-only = 0`.
- Shared-coordinate bit score, E-value, pident, length, mismatch, and gapopen
  diffs are all zero.
- Output order and formatting are verified for the release output mode.

## Phase 4: Diagnose The Megablast Statistic Diffs

Focus on `Sakai.MG1655.megablast`, because hit count, coordinates, bit score,
and E-value already match.

Trace the shared-coordinate HSPs reported by the comparison script:

- pident diff:
  `q=495212-495350 s=490068-489930`
- pident diff:
  `q=4944798-4944906 s=3330459-3330569`
- length/mismatch/gapopen diff:
  `q=48267-48485 s=1781063-1781281`
- length/mismatch/gapopen diff:
  `q=4663939-4664042 s=2684219-2684116`

Compare LOSAT and NCBI at the edit-script level:

1. greedy traceback raw edit operations.
2. edit-script trimming after common-endpoint purge.
3. ambiguity re-evaluation edits and deletion decisions.
4. identity, mismatch, gap-open, and alignment-length accumulation.
5. output formatting of `pident`.

Required NCBI source checks:

- `greedy_align.c` for greedy traceback script construction.
- `blast_gapalign.c` for megablast greedy call timing.
- `blast_traceback.c` for post-traceback update, identity/length testing, and
  ambiguity re-evaluation order.
- `blast_hits.c` for edit-script trimming and adjusted output offsets.

Exit criteria:

- The first edit-script/statistic divergence is known for each representative
  HSP.
- The fix target is tied to an NCBI source function and exact line range.

## Phase 5: Patch Megablast Edit-Script Statistics

Patch the proven megablast offender without changing task-blastn behavior unless
NCBI uses the same shared helper for both modes.

Verification after each patch:

```bash
python3 tests/compare_blastn_parity.py \
  --manifest tests/blastn_parity_manifest.tsv \
  --limit 3
```

Exit criteria for megablast:

- `EDL933.Sakai.megablast` remains exact.
- `Sakai.MG1655.megablast` has:
  - equal hit count.
  - equal coordinate keys.
  - zero shared-coordinate bit-score diffs.
  - zero shared-coordinate E-value diffs.
  - zero shared-coordinate pident diffs.
  - zero shared-coordinate length/mismatch/gapopen diffs.
- Output order and formatting are verified for the release output mode.

## Phase 6: Strict Output And Release Evidence

Once the manifest summary is clean, verify release-facing output:

1. Compare exact output lines for every BLASTN manifest row.
2. If the field-wise comparator ignores comments, add a separate formatting
   check for outfmt `7` comments and field ordering.
3. Confirm outfmt `6` and outfmt `7` behavior for BLASTN if both are claimed in
   release scope.
4. Confirm multithreaded LOSAT output reduces to the same order if BLASTN
   threaded output is claimed.

Required local gate from a clean release branch:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
cargo build --release
cargo package --list
cargo publish --dry-run --locked \
  --config 'build.target-dir="/tmp/losat-cargo-publish-target"'
```

Required BLASTN release evidence:

- Updated `docs/release/v0.1.0.md` scope table, changing BLASTN from
  experimental to supported only after all checks pass.
- NCBI BLAST+ version.
- LOSAT commit.
- Manifest command lines.
- Input checksums.
- Output checksums.
- Diff summaries showing zero differences.
- Any approved deviations. For BLASTN, no deviation is currently approved.

## Final Acceptance Checklist

BLASTN may be promoted to supported only when all are true:

- `tests/blastn_parity_manifest.tsv` has 11/11 exact cases.
- task-blastn has `NCBI-only = 0` and `LOSAT-only = 0` for all 9 cases.
- megablast has zero pident, length, mismatch, and gapopen diffs for both cases.
- bit score, E-value, coordinates, hit order, and formatting match NCBI for the
  release output modes.
- All Rust behavior changes include NCBI source comments with path, line
  numbers, and C/C++ snippets.
- No LOSAT runtime, build, fallback, or unsupported path invokes NCBI BLAST.
- The final evidence is produced from a clean release branch.
- The release note, README, CHANGELOG, and CLI help agree on the promoted BLASTN
  scope.

If any item remains false, keep BLASTN experimental for v0.1.0.
