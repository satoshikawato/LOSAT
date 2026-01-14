# AGENTS.md

Instructions for coding agents working in this repository. This file is the
authoritative, current guidance for agent behavior in LOSAT.

---

## Mandatory Compliance Requirements (Absolute)

1. NCBI BLAST is the only source of truth.
   - Use the NCBI C/C++ implementation as the sole authoritative reference.
   - Do not assume or guess behavior; always refer to the NCBI source code.
   - If the corresponding NCBI code cannot be found, the feature must not exist.

2. Bit-perfect output parity is required.
   - Output must match NCBI BLAST+ exactly.
   - Do not simplify algorithms if it changes output.
   - Use the same floating-point precision as NCBI.

3. NCBI code comments are mandatory for modifications.
   - Every code change must include NCBI C/C++ reference comments with file path
     and line numbers.
   - Include the relevant NCBI snippet immediately above the Rust code.
   - If you cannot add the NCBI reference comments, do not write the code.

4. No unauthorized features.
   - Never introduce functionality that does not exist in NCBI.
   - If a feature lacks an NCBI equivalent, delete it immediately.

5. Exhaustive difference resolution.
   - Identify and fix every discrepancy between LOSAT and NCBI BLAST.
   - Patch every offender in one sweep, then fix compile issues.

6. No assumptions or guessing.
   - Read the NCBI source; never speculate.

7. Correct timing, order, and context.
   - Call algorithms at the exact same timing and order as NCBI.
   - Input/output data must match NCBI exactly.

8. Algorithmic fidelity over aesthetics.
   - Do not refactor for readability if it changes logic.
   - Rust-specific deviations are allowed only to satisfy the memory model or
     performance while preserving behavior and Big O complexity.

9. Testing discipline.
   - Do not run "maybe helpful" tests while known NCBI divergences remain.
   - Run tests only when requested or after completing a parity sweep.

---

## Current High-Priority Issues (Keep Updated)

### TBLASTX long-sequence 2x hits
- Long sequences (600kb+, gencode=4) produce ~2x hits; short sequences are near parity.
- Excess hits cluster in ~22-30 bit scores, implying extension boundary or
  reevaluation differences.
- Investigate extension boundaries and `Blast_HSPListReevaluateUngapped` parity.

### TBLASTX chain formation differences
- Short HSPs sometimes link into higher-score chains, causing E-value mismatches.
- Focus on `link_hsps.c` parity and chain member filtering timing.

### BLASTN coverage gap (blastn task)
- Megablast parity is higher; blastn task shows ~90% coverage on some datasets.
- Suspect gapped DP boundary/coordinate handling (`Blast_SemiGappedAlign`,
  `extend_gapped_one_direction_ex`) in `blast_gapalign.c`.
- Recent parity work aligned query-context offsets, output coordinate adjustment,
  x-drop clamping, and hitlist pruning with NCBI; re-run comparisons before
  digging deeper.

---

## Critical Parity Notes (TBLASTX)

- Sequence encoding: nucleotides are 2-bit packed; amino acids use NCBISTDAA with
  sentinel byte 0 (NULLB); BLOSUM62 is the scoring matrix.
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

---

## Debug/Diagnostics Environment Variables

- `LOSAT_TRACE_HSP="qstart,qend,sstart,send"` trace a specific HSP.
- `LOSAT_DEBUG_CUTOFFS=1` cutoff calculations (tblastx + blastn).
- `LOSAT_DEBUG_CHAINING=1` chaining debug (legacy; tblastx).
- `LOSAT_DEBUG_EXTENSION=1` tblastx extension debug.
- `LOSAT_DEBUG_BLASTN=1` blastn hit loss diagnostics.
- `LOSAT_DEBUG_COORDS=1` blastn coordinate transforms.
- `LOSAT_DEBUG_SCAN_SOFF=<int>` tblastx scan debug center subject offset.
- `LOSAT_DEBUG_SCAN_WINDOW=<int>` tblastx scan debug window size.
- `LOSAT_TIMING=1` timing breakdown.
- `LOSAT_DIAGNOSTICS=1` general diagnostics counters.
- `LOSAT_STARTUP_TRACE=1` startup trace.

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
│   ├── algorithm/         # Core algorithms (tblastx, blastn, common)
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
Primary focus: TBLASTX and BLASTN. The Rust crate root is `LOSAT/`.
