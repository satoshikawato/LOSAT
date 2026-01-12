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
   - If you cannot add the NCBI reference comments, do not write the code.

4. No unauthorized features.
   - Never introduce functionality that does not exist in NCBI.
   - If a feature lacks an NCBI equivalent, delete it immediately.

5. Exhaustive difference resolution.
   - Identify and fix every discrepancy between LOSAT and NCBI BLAST.
   - Verify bit-perfect matching after changes.

6. No assumptions or guessing.
   - Read the NCBI source; never speculate.

7. Correct timing, order, and context.
   - Call algorithms at the exact same timing and order as NCBI.
   - Input/output data must match NCBI exactly.

---

## Current High-Priority Issues (Keep Updated)

### BLASTN parity gap (Session 3, 2026-01-09)
- Query coordinates off by 1 in self-comparison after removing non-NCBI
  workaround.
- Root suspect: `extend_gapped_one_direction_ex` asymmetry between query and
  subject handling.
- Must fix by matching NCBI `blast_gapalign.c:735-962`.

### Non-NCBI code to remove (must delete)
- `LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs`:
  `filter_contained_hsps` is non-NCBI and caused large hit loss.
- Remove any imports/exports and disabled calls in `blast_engine/mod.rs`.

### TBLASTX 2x hit issue (long sequences)
- LOSAT produces ~2x hits for long sequences (gencode=4), near parity for short
  sequences.
- LOSAT bit scores are ~7 bits higher at identical coordinates, implying longer
  alignments or score accumulation differences.
- Most likely root: extension boundary/termination or reevaluation step.
- Investigate `Blast_HSPListReevaluateUngapped` and extension boundaries.

---

## Debug/Diagnostics Environment Variables

- `LOSAT_TRACE_HSP="qstart,qend,sstart,send"` trace a specific HSP.
- `LOSAT_DEBUG_CUTOFFS=1` cutoff calculations.
- `LOSAT_DEBUG_CHAINING=1` chaining debug.
- `LOSAT_DEBUG_EXTENSION=1` TBLASTX extension debug.
- `LOSAT_DEBUG_BLASTN=1` BLASTN hit loss diagnostics.
- `LOSAT_DEBUG_COORDS=1` BLASTN coordinate transforms.
- `LOSAT_TIMING=1` timing breakdown.
- `LOSAT_NO_SIMD=1` force scalar processing.
- `LOSAT_DIAGNOSTICS=1` general diagnostics.
- `LOSAT_STARTUP_TRACE=1` startup trace.

---

## Build and Test (from repo root)

```bash
cd LOSAT && cargo build --release
cd LOSAT && cargo test
cd LOSAT && cargo clippy
cd LOSAT && cargo fmt
```

## Integration Tests

```bash
cd LOSAT/tests && bash run_comparison.sh
cd LOSAT/tests && bash run_all_tests.sh
```

---

## Codebase Overview (Concise)

```
LOSAT/src/
├── main.rs                # CLI entry; dispatch to blastn/tblastx
├── algorithm/             # Core algorithms
│   ├── tblastx/            # TBLASTX pipeline
│   ├── blastn/             # BLASTN/megablast pipeline
│   └── common/             # Shared algorithm utilities
├── core/                  # NCBI-ported primitives (stats, filters, encoding)
├── report/                # Output formatting (outfmt0/6/7)
├── blastinput/            # CLI argument parsing
├── sequence/              # Sequence storage/encoding
└── utils/                 # Shared utilities and tables
```

---

## Key References

### NCBI Source Files (examples)
- `c++/src/algo/blast/core/aa_ungapped.c`
- `c++/src/algo/blast/core/na_ungapped.c`
- `c++/src/algo/blast/core/link_hsps.c`
- `c++/src/algo/blast/core/blast_parameters.c`
- `c++/src/algo/blast/core/blast_query_info.c`
- `c++/src/algo/blast/core/blast_stat.c`
- `c++/src/algo/blast/core/greedy_align.c`
- `c++/src/algo/blast/core/blast_gapalign.c`

### LOSAT Files of Interest
- `LOSAT/src/algorithm/blastn/alignment/gapped.rs`
- `LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs`
- `LOSAT/src/algorithm/tblastx/utils.rs`
- `LOSAT/src/algorithm/tblastx/lookup.rs`
- `LOSAT/src/algorithm/tblastx/extension.rs`
- `LOSAT/src/algorithm/tblastx/ncbi_cutoffs.rs`

---

## Project Summary

LOSAT is a Rust reimplementation of NCBI BLAST targeting bit-perfect parity.
Primary focus: TBLASTX and BLASTN.
