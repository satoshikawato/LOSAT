# AGENTS.md

Instructions for coding agents working in this repository. This file is the
authoritative, current guidance for agent behavior in LOSAT.

---

## Mandatory Compliance Requirements (Absolute)

1. NCBI BLAST is the only source of truth.
   - Use the NCBI C/C++ implementation as the sole authoritative reference.
   - Do not assume or guess behavior; always refer to the NCBI source code.
   - If the corresponding NCBI code cannot be found, the feature must not exist.

2. NCBI BLAST is reference-only, never an implementation dependency.
   - NCBI BLAST+ binaries, libraries, FFI bindings, and subprocess calls may be
     used only for source inspection, comparison tests, diff analysis, and
     validation-oracle runs.
   - Do not invoke `blastp`, `blastn`, `tblastx`, or any other NCBI executable
     from LOSAT runtime code, build scripts, feature code paths, or fallback
     paths.
   - Do not delegate unsupported, unported, or partially ported behavior to
     NCBI BLAST. If LOSAT does not implement a feature in Rust yet, it must
     fail fast with an explicit unsupported/unimplemented error and remain on
     the Rust porting backlog.
   - Do not link against, embed, or otherwise depend on NCBI BLAST code as a
     runtime implementation component. Reference the source; do not ship the
     behavior by calling out to NCBI.

3. Bit-perfect output parity is required.
   - Output must match NCBI BLAST+ exactly.
   - Do not simplify algorithms if it changes output.
   - Use the same floating-point precision as NCBI.
   - Approved project exception: for TBLASTX local `-s/--subject` searches,
     LOSAT must respect `--db-gencode` for subject translation/search/reporting
     for every non-default genetic code. This local-subject non-default
     `db_gencode` behavior is intentional even where NCBI BLAST+ local
     `-subject` semantics differ. Do not count those subject-genetic-code-only
     differences as LOSAT parity defects. This exception is narrow and does not
     permit any other deviation from NCBI timing, ordering, scoring, filtering,
     statistics, pruning, or output formatting.

4. NCBI code comments are mandatory for modifications.
   - Every code change must include NCBI C/C++ reference comments with file path
     and line numbers.
   - Include the relevant NCBI snippet immediately above the Rust code.
   - If you cannot add the NCBI reference comments, do not write the code.

5. No unauthorized features.
   - Never introduce functionality that does not exist in NCBI.
   - If a feature lacks an NCBI equivalent, delete it immediately.

6. Exhaustive difference resolution.
   - Identify and fix every discrepancy between LOSAT and NCBI BLAST.
   - Patch every offender in one sweep, then fix compile issues.

7. No assumptions or guessing.
   - Read the NCBI source; never speculate.

8. Correct timing, order, and context.
   - Call algorithms at the exact same timing and order as NCBI.
   - Input/output data must match NCBI exactly.

9. Algorithmic fidelity over aesthetics.
   - Do not refactor for readability if it changes logic.
   - Rust-specific deviations are allowed only to satisfy the memory model or
     performance while preserving behavior and Big O complexity.

10. Testing discipline.
   - Do not run "maybe helpful" tests while known NCBI divergences remain.
   - Run tests only when requested or after completing a parity sweep.

---

## Current Status and Active Work (Keep Updated)

### Active parity work
- BLASTN task parity remains active. Re-run current comparison fixtures before
  diagnosing; do not rely on old hit-count percentages. Recent work aligned
  query-context offsets, output coordinate adjustment, x-drop clamping, common
  endpoint purge ordering, ambiguity re-evaluation, edit-script trimming, and
  hitlist pruning with NCBI. Remaining BLASTN work should focus only on a
  reproducible current fixture and the corresponding NCBI source path.
- BLASTP is implemented but remains secondary to TBLASTX/BLASTN. Treat BLASTP
  parity and performance work as ongoing unless a current comparison fixture
  proves the exact behavior being touched. Unsupported or incomplete BLASTP
  behavior must fail fast; never delegate to external BLASTP.
- Wasm/threading work is implementation-level only. Plain `wasm32-wasip1` builds
  are intentionally serial; real command-Wasm threading requires the
  `wasm32-wasip1-threads` target and `wasm-threads` feature. Any native or Wasm
  work partitioning must reduce results back into the NCBI order.

### Resolved regression targets
- The TBLASTX long-sequence AP027131/AP027133 gencode-4 local `-subject` 2x-hit
  issue is resolved for the approved local-subject `db_gencode` behavior. Keep
  this as a regression fixture; do not reopen the old extension-boundary or
  reevaluation hypothesis without a new current diff.
- The short LC738874/LC738875 TBLASTX threshold sweep reached exact NCBI parity
  for hit counts, coordinates, E-values, and bit scores at `-evalue 10`, `100`,
  and `10000`. Do not create a generic active TBLASTX chaining issue unless a
  new fixture reproduces one.
- Chain member filtering remains a critical parity rule, not an open issue:
  filter `linked_set && !start_of_chain` during output, not during linking.

### Performance work boundaries
- External implementations such as DIAMOND may be read for data layout,
  scheduling, cache locality, and SIMD implementation ideas only. They are not
  parity or behavior authorities.
- Do not import DIAMOND-style spaced seeds, minimizers, frequency masking,
  sensitivity rounds, ungapped prefilters, or candidate-pruning heuristics into
  LOSAT unless NCBI BLAST has the same behavior for the same program/task.
- Performance changes must preserve the same candidate set, HSP construction,
  pruning, ordering, statistics, and formatting as NCBI. If a faster path cannot
  be proven byte-identical on the relevant fixtures, keep it disabled or remove
  it.

---

## Critical Parity Notes (TBLASTX)

- Sequence encoding: nucleotides are 2-bit packed; amino acids use NCBISTDAA with
  sentinel byte 0 (NULLB); BLOSUM62 is the scoring matrix.
- Local `-s/--subject` TBLASTX searches intentionally honor `--db-gencode` for
  translated subject search/reporting for every non-default genetic code. Do not
  use NCBI BLAST+ local `-subject` non-default `db_gencode` output as a failure
  oracle for that subject-genetic-code behavior; use NCBI as the oracle for all
  other behavior in the same run.
- Frame concatenation shares boundary sentinels; frame offset advances by
  `aa_len + 1`, not `aa_len + 2`.
- Length adjustment asymmetry: query uses full adjustment, subject uses one third;
  effective search space uses full adjustment for both.
- Masking: SEG applies to query only; extension uses masked sequence while identity
  uses unmasked sequence.
- Subject frame sort order: negative frames come first (ascending frame value).
- Chain member filtering: filter `linked_set && !start_of_chain` during output,
  not during linking.
- Cutoff score capping: `min(BLAST_Cutoffs, gap_trigger, cutoff_score_max)`;
  tblastx uses `scale_factor = 1.0`.

---

## Critical Parity Notes (BLASTN)

- Query contexts: blastn uses plus and minus per query; context index is
  `q_idx * 2 + strand`, and context offset advances by `query_len + 1`.
- Total query length for interval tree and offsets is `2 * query_len + 1`.
- Subject is plus-only for blastn (no reverse-complement subject).
- Output coordinates must follow `Blast_HSPGetAdjustedOffsets` logic; when query
  is minus, flip subject coords while keeping internal subject offsets canonical.
- HSP pruning/comparisons use internal (contexted) offsets; adjust to output
  coords after pruning.
- Gapped DP x-drop uses `min(x_dropoff, ungapped_score)` for the trace cutoff.
- Hitlist pruning follows NCBI: trim by `max_hsps_per_subject`, then apply
  `min(hitlist_size, max_target_seqs)` with NCBI score/evalue ordering.
- `SCAN_RANGE_BLASTN` is 0 (no scan range for blastn tasks).
- Common-endpoint purge pass-1 uses `s_QueryOffsetCompareHSPs` tie-breaker
  (query/subject end DESC on score ties) to keep longer HSPs.
- Post-traceback filtering mirrors `blast_traceback.c`: re-sort by
  `ScoreCompareHSPs`, then interval-tree containment purge
  (`BlastIntervalTreeContainsHSP`).
- Re-evaluation uses canonical (plus) subject sequence even when output
  coordinates are on the minus strand.

---

## Debug/Diagnostics Environment Variables

- `LOSAT_TRACE_HSP="qstart,qend,sstart,send"` trace a specific TBLASTX HSP.
- `LOSAT_TRACE_HSP_MASKS=1` print mask coverage for the traced TBLASTX HSP.
- `LOSAT_TRACE_CHAIN_HSP="qstart,qend,sstart,send"` trace TBLASTX chain
  selection for a specific HSP.
- `LOSAT_TRACE_LINK_SELECTIONS=1` print TBLASTX link-selection details.
- `LOSAT_DUMP_TBLASTX_STAGE=<dir>` append TBLASTX stage snapshots as TSV files.
- `LOSAT_TRACE_BLASTN_HSP="qstart,qend,sstart,send"` trace a specific BLASTN
  HSP.
- `LOSAT_TRACE_BLASTN_SEED="q,s"` trace a specific BLASTN seed.
- `LOSAT_TRACE_BLASTN_CONTEXT=<context_idx>` restrict BLASTN tracing by context.
- `LOSAT_TRACE_BLASTN_SUBJECT=<subject_id_or_index>` restrict BLASTN tracing by
  subject.
- `LOSAT_TRACE_BLASTN_STAGE=<seed|ungapped|prelim|traceback|purge|hitlist|all>`
  restrict BLASTN tracing by stage.
- `LOSAT_DEBUG_CUTOFFS=1` cutoff calculations (tblastx + blastn).
- `LOSAT_DEBUG_CUTOFFS_ALL=1` verbose TBLASTX cutoff diagnostics.
- `LOSAT_DEBUG_CHAINING=1` chaining debug (legacy; tblastx).
- `LOSAT_DEBUG_EXTENSION=1` tblastx extension debug.
- `LOSAT_DEBUG_HSP_SAVING=1` TBLASTX HSP-save diagnostics.
- `LOSAT_DEBUG_OUTPUT_FILTER=1` TBLASTX output filter diagnostics.
- `LOSAT_DEBUG_BLASTN=1` blastn hit loss diagnostics.
- `LOSAT_DEBUG_COORDS=1` blastn coordinate transforms.
- `LOSAT_DEBUG_COORDS_START=<int>` narrow selected BLASTN coordinate diagnostics.
- `LOSAT_DEBUG_SCAN_SOFF=<int>` tblastx scan debug center subject offset.
- `LOSAT_DEBUG_SCAN_WINDOW=<int>` tblastx scan debug window size.
- `LOSAT_TIMING=1` timing breakdown.
- `LOSAT_DIAGNOSTICS=1` general diagnostics counters.
- `LOSAT_STARTUP_TRACE=1` startup trace.
- `LOSAT_WASI_THREADS_DEBUG=1` threaded-WASI scheduling diagnostics.
- `LOSAT_TBLASTX_PARALLEL_CHUNKS=1` force TBLASTX subject-chunk parallel path
  for diagnostics.
- `LOSAT_TBLASTX_PARALLEL_SCAN_CHUNKS=1` diagnostic-only TBLASTX scan-interior
  chunking; do not enable by default without parity proof.

---

## Build, Lint, Format (from repo root)

```bash
cd LOSAT && cargo build --release
cd LOSAT && cargo test
cd LOSAT && cargo clippy
cd LOSAT && cargo fmt
```

## Testing Expectations

- Add unit tests for NCBI-ported functions, including edge cases and boundaries.
- Reference NCBI unit tests when available:
  `ncbi-blast/c++/src/algo/blast/unit_tests/`.
- Integration tests must compare output with NCBI BLAST+ and verify hit counts,
  bit scores, E-values, and coordinates.
- NCBI BLAST+ execution is allowed only as a comparison oracle during testing
  and diagnostics; it must not be part of LOSAT's feature implementation,
  runtime execution, build pipeline, fallback handling, or unsupported-feature
  path.
- Hit-count deltas are a diagnostic only; if tracked, use <0.2% as a trend
  threshold, but acceptance still requires bit-perfect parity.

## Integration and Comparison Scripts

```bash
cd LOSAT/tests && bash run_comparison.sh
cd LOSAT/tests && bash run_all_tests.sh
bash compare_tblastx_results.sh
bash tests/compare_self_tblastx.sh
bash tests/compare_long_sequences_debug.sh
bash docs/compare_seg_mask.sh
```

---

## Repository Layout (Concise)

```
LOSAT/                     # Rust crate root
├── Cargo.toml
├── src/
│   ├── main.rs            # CLI entry; dispatch to blastn/tblastx
│   ├── algorithm/         # Core algorithms (tblastx, blastn, blastp, common)
│   ├── core/              # NCBI-ported primitives (stats, filters, encoding)
│   ├── stats/             # Karlin-Altschul and sum statistics
│   ├── align/             # Alignment utilities/traceback
│   ├── report/            # Output formatting (outfmt 0/6/7)
│   ├── format/            # Format helpers
│   ├── post/              # Post-processing (chaining/filtering)
│   ├── api/               # Options and API layer
│   ├── blastinput/        # CLI argument parsing
│   ├── config/            # Compatibility/config helpers
│   ├── seed/              # Word finding
│   ├── sequence/          # Sequence storage/encoding
│   └── utils/             # Shared utilities and tables
└── tests/
    ├── run_comparison.sh  # Compare vs NCBI BLAST+
    ├── run_all_tests.sh
    ├── unit/              # Rust unit tests
    ├── blast_out/ ncbi_out/ losat_out/ fasta/
    └── plots/
```

Additional scripts and datasets live at repo root: `tests/`, `losat_out/`, and
`compare_tblastx_results.sh`.

---

## Entry Points

- CLI: `LOSAT/src/main.rs`
- TBLASTX engine: `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`
- BLASTN engine: `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- BLASTP engine: `LOSAT/src/algorithm/blastp/blast_engine.rs`

---

## Key References

### NCBI Source Locations
- Primary NCBI repo: `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/`
- Alternate NCBI repo: `/mnt/c/Users/genom/GitHub/ncbi-blast/`

### Local NCBI Snapshot (subset, for quick lookup)
- `.ncbi_ref/blast_gapalign.c`
- `.ncbi_ref/blast_engine.c`
- `.ncbi_ref/blast_parameters.c`
- `.ncbi_ref/blast_setup.c`
- `.ncbi_ref/na_ungapped.c`
- `.ncbi_ref/greedy_align.c`
- `.ncbi_ref/blast_encoding.c`

### NCBI Source Files (examples)
- `c++/src/algo/blast/core/aa_ungapped.c`
- `c++/src/algo/blast/core/link_hsps.c`
- `c++/src/algo/blast/core/blast_parameters.c`
- `c++/src/algo/blast/core/blast_query_info.c`
- `c++/src/algo/blast/core/blast_stat.c`
- `c++/src/algo/blast/core/greedy_align.c`
- `c++/src/algo/blast/core/blast_gapalign.c`

---

## Conventions and Porting Notes

- Match NCBI function names when porting (example: `s_BlastAaExtendTwoHit` ->
  `extend_hit_two_hit`).
- Keep NCBI terminology in comments: query_offset, subject_offset, context, frame.
- Use `#[inline]` for hot-path functions; use `unsafe` only when necessary and
  always add safety comments.
- Prefer `cargo fmt` and `cargo clippy` conventions for Rust style.

---

## Other Guidance Files (Reference Only)

- `CLAUDE.md`
- `cursorrules.mdc`
- `.cursor/rules/cursorrules.mdc`
- `.cursor/rules/global_rules.mdc`

This file consolidates their requirements; if there is a conflict, follow this file.

---

## Project Summary

LOSAT is a Rust reimplementation of NCBI BLAST targeting bit-perfect parity.
Primary focus: TBLASTX and BLASTN. BLASTP support exists but remains secondary
and must be treated as ongoing parity/performance work unless the touched path is
covered by current comparison fixtures. The Rust crate root is `LOSAT/`.
Incomplete areas must not be filled by delegating to external NCBI BLAST
executables or libraries.
