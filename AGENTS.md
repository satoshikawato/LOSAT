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
   - Approved TBLASTN-only product decision (`PD-TLOSAN-LOCAL-GENCODE-32`):
     local `-subject` searches must apply the selected `-db_gencode` to subject
     translation, candidate search, HSP re-evaluation, scoring, statistics,
     coordinates, and `outfmt 0/6/7` reporting for all 27 NCBI `gc.prt`
     genetic-code IDs. The TBLASTN CLI must accept ID 32 even though the pinned
     NCBI BLAST+ CLI rejects it; verify ID 32 with a comparison-only NCBI C++
     API oracle using `FindGeneticCode(32)` through search and formatting.
     Reject invalid IDs explicitly. Differences from NCBI local `-subject`
     caused solely by honoring a non-default subject code are permitted; no
     difference in call timing, ordering, candidate rules, linking, filtering,
     statistical formulas, or formatting is permitted. This does not expand
     the TBLASTX exception or authorize NCBI as a runtime/build dependency.
   - Approved CLI exceptions outside the search results, for every program
     (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`): argument-parser syntax errors and
     `-help`/`-h` use LOSAT's parser text and exit code 2 (NCBI: USAGE, exit
     1); with `-subject`, LOSAT honors `-num_threads` and does not print NCBI's
     two thread warnings; a write failure in outfmt 6/7, where NCBI aborts,
     is reported with a non-zero exit; an allocation failure aborts (NCBI:
     "BLAST ran out of memory", exit 4); an outfmt 0 write to a closed pipe,
     where NCBI is ended by SIGPIPE, is reported as "BLAST failed to write
     output", exit 6, and a closed pipe that LOSAT's completed writes do not
     reach ends LOSAT with exit 0 (timing; NCBI's later flush gets SIGPIPE);
     a standard output closed at the start (`>&-`), where NCBI's first write
     fails (outfmt 0 exit 6, outfmt 6/7 abort), is the `/dev/null` that the
     Rust runtime opens before `main`, so LOSAT discards the report and exits
     as the search ends. Errors raised after argument parsing, the outfmt 0 write
     failure ("BLAST failed to write output", exit 6), all other warnings, and
     every search result must still match NCBI. Non-UTF-8
     file names, `.ncbirc` keys that change output, and NCBI C++ Toolkit
     words (`-version`, `-dryrun`, `-logfile`, ...; also in an option's value)
     are explicit rejections, not exceptions.
   - Approved exceptions for NCBI defects (`PD-LOSAT-NCBI-DEFECTS`): where
     NCBI BLAST+ fails (a crash, a debug-only assumption, a read past a buffer)
     and LOSAT's result was shown to equal NCBI's output for a nearby input
     that does not reach the defect. BLASTN: a query chunk that NCBI would
     split again (CHUNK_SIZE/OVERLAP_CHUNK_SIZE; NCBI stops with a
     CCoreException) is searched once. BLASTN, TBLASTX and TBLASTN: outfmt 0
     titles made only of punctuation stop at the end of the string (NCBI
     reads past it and crashes). BLASTP: a one-hit gapped start whose window
     of 11 letters passes the end of a sequence reads the letters past it as
     the sentinel (NCBI reads past its buffer and can crash; LOSAT's output
     equals NCBI's under valgrind). Deterministic NCBI results, even wrong-looking ones, are
     reproduced, not excepted; NCBI failures without a checkable valid result
     are explicit rejections.

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

## Parity boundaries and regression targets

### Program and target boundaries
- PR 6 native certification has two non-interchangeable contracts under
  `PD-NCBI-PLATFORM-VARIANCE`: Gate A requires every LOSAT result to match the
  frozen PR 5 raw bytes; Gate B requires each of six official NCBI searches to
  match its exact registered platform fingerprint. Never restore the old
  requirement that LOSAT equal platform-local NCBI bytes when those registered
  fingerprints differ. Native fingerprints never become LOSAT expected output,
  never authorize platform-specific LOSAT behavior, and never create a new
  parity exception. An unknown native fingerprint hard-fails and requires a new
  characterized authority version and review.
- Re-run current comparison fixtures before diagnosing BLASTN. Do not rely on
  old hit-count percentages or session notes. Start from one reproducible
  current fixture and the corresponding NCBI source path.
- BLASTP is secondary to TBLASTX and BLASTN. A current comparison fixture must
  prove the exact behavior being changed. Unsupported or incomplete BLASTP
  behavior must fail fast; never delegate it to external BLASTP.
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

### Benchmark target protocol
- Treat repository benchmark figures as compute-program benchmarks, not as
  biological interpretation or downstream bioinformatics analysis.
- While a long benchmark command is running, poll its status at ten-minute
  intervals to conserve agent/tool tokens. Short expected completions and
  immediate failure diagnosis are the only exceptions; do not busy-poll.
- The standard formal benchmark protocol is one untimed warmup followed by
  exactly three retained timed repetitions for every case and execution mode.
  Report the median and the full three-sample min-max range; never select the
  fastest sample. Add repetitions only when the user explicitly requests them
  or the three retained samples are demonstrably inconclusive.
- For NCBI BLAST+ execution-time measurements, prepare the subject with
  `makeblastdb` before timing and run every timed BLASTN, BLASTP, and TBLASTX
  search with `-db`. Exclude database-construction time from search wall time and
  retain the database-build command, version, input checksum, and elapsed time as
  separate provenance. Do not use `-subject` timings as multithreaded NCBI
  results: NCBI reduces local-subject searches to one thread.
- For hit-distribution figures, use NCBI `-subject` output for BLASTN and BLASTP,
  but use NCBI `-db` output for TBLASTX. TBLASTX database searches must pass the
  matching `-db_gencode`, including genetic code 4, so the plotted NCBI subject
  translation is comparable to LOSAT's intentional local-subject behavior.
- Keep timing targets and distribution targets in separately named outputs and
  metadata. Never relabel a local-subject NCBI output as a threaded timing or a
  BLASTN/BLASTP database output as the local-subject distribution oracle.

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

The `LOSAT_TRACE_*`, `LOSAT_DEBUG_*`, `LOSAT_TIMING`, `LOSAT_DIAGNOSTICS` and related
variables are listed in [`docs/agents/diagnostics.md`](docs/agents/diagnostics.md).

---

## Build, Lint, Format (from repo root)

```bash
cd LOSAT && cargo build --release
cd LOSAT && cargo test
cd LOSAT && cargo clippy
cd LOSAT && cargo fmt
```

## Testing Expectations

- Use `$verify-ncbi-parity-and-speed` for parity, benchmark, native/Wasm, and
  release-evidence work (Claude Code: skill `verify-ncbi-parity-and-speed`).
- Ask the `ncbi_parity_auditor` custom agent for an independent read-only check
  before accepting a release-facing parity or performance claim (Claude Code:
  agent `losat-reviewer`, independent audit).
- Choose verification by what the change touches; Claude Code's tiers (quick,
  standard, full) are in the skill `losat-gates`. Rule 10 above governs parity test
  runs; `web/AGENTS.md` governs checks for application commits.
- Add unit tests for NCBI-ported functions, including edge cases and boundaries.
- Reference NCBI unit tests when available:
  `ncbi-blast/c++/src/algo/blast/unit_tests/`.
- Port the NCBI unit-test cases of a ported function one case at a time, inside the parity
  sweep of the module that owns it (rule 10). The test lives next to the implementation, in
  its `#[cfg(test)]` module or a sibling `<name>_tests.rs` declared with
  `#[cfg(test)] mod <name>_tests;`, with the line
  `// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/<dir>/<file>:<lines> <Case>`
  above it and a row in `docs/evidence/ncbi_unit_cases/LEDGER.tsv` (statuses: ported,
  partial, e2e, to-port, n-a, superseded; `python3 LOSAT/tests/ncbi_unit_case_ledger.py`).
  Do not add integration binaries under `LOSAT/tests/` for unit tests (each links the whole
  crate and sees only `pub` items) and do not make private items `pub` for a test. Cases that
  need the object manager, GenBank, BLAST databases, or PSI/RPS/PHI are `n-a`; whole-search
  cases (bl2seq) become fixture rows, not unit tests.
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

```
web/                       # LOSAT Web (browser application)
├── AGENTS.md              # Rules for code under web/
├── adapter/               # Rust crate: Wasm ABI v2 over the LOSAT library (added in session S05)
└── app/                   # Vite + Vue 3 + TypeScript application
```

### Scope of these rules for `web/`

Code under `web/` follows [`web/AGENTS.md`](web/AGENTS.md) under
[`PD-LOSAT-WEB-APP-BOUNDARY`](docs/product_decisions/PD-LOSAT-WEB-APP-BOUNDARY.md).
Application code there must not compute or change BLAST behavior or compatibility
outputs; it therefore does not carry NCBI reference comments, and application features
permitted by that decision are not "unauthorized features" under requirement 5. Every
change under `LOSAT/`, including one made for LOSAT Web, remains subject to all the
requirements above.

---

## Entry Points

- CLI: `LOSAT/src/main.rs`
- TBLASTX engine: `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`
- BLASTN engine: `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- BLASTP engine: `LOSAT/src/algorithm/blastp/blast_engine.rs`

---

## Key References

### NCBI Source Locations
- Pinned NCBI C/C++ source: `satoshikawato/ncbi-blast` at commit
  `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4` (NCBI BLAST+ 2.17.0). The local
  checkout path is machine-specific (Claude Code: `$NCBI_SRC` in `CLAUDE.local.md`).

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
