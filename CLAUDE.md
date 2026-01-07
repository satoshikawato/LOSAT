# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

```bash
# Build (from LOSAT/ subdirectory)
cd LOSAT && cargo build --release

# Run tests
cd LOSAT && cargo test

# Run a single test
cd LOSAT && cargo test test_name

# Run tests in a specific module
cd LOSAT && cargo test tblastx::

# Lint
cd LOSAT && cargo clippy

# Format
cd LOSAT && cargo fmt
```

## CLI Usage

```bash
# TBLASTX (translated nucleotide vs translated nucleotide)
./target/release/LOSAT tblastx -q query.fasta -s subject.fasta -o output.txt

# BLASTN (nucleotide vs nucleotide)
./target/release/LOSAT blastn -q query.fasta -s subject.fasta -o output.txt

# Common options:
#   -e, --evalue <EVALUE>       E-value threshold (default: 10.0)
#   -n, --num-threads <N>       Number of threads (0 = auto)
#   --query-gencode <N>         Genetic code for query (default: 1)
#   --db-gencode <N>            Genetic code for subject (default: 1)
#   --outfmt <FORMAT>           Output format (6=tabular, 7=tabular with headers)
```

## Integration Tests

```bash
# Run comparison tests against NCBI BLAST (from LOSAT/tests/ directory)
cd LOSAT/tests && bash run_comparison.sh

# Run all tests including plots
cd LOSAT/tests && bash run_all_tests.sh
```

## Project Overview

LOSAT (LOcal Sequence Alignment Tool) is a Rust reimplementation of NCBI BLAST targeting **bit-perfect parity** with NCBI BLAST+ output. Primary focus is TBLASTX.

## Architecture

```
LOSAT/src/
├── algorithm/
│   ├── tblastx/           # TBLASTX implementation (primary focus)
│   │   ├── utils.rs       # Main execution: run(), run_with_neighbor_map()
│   │   ├── lookup.rs      # Lookup table, context management
│   │   ├── extension.rs   # Two-hit extension, X-drop termination
│   │   ├── sum_stats_linking.rs  # HSP chaining, E-value calculation
│   │   ├── ncbi_cutoffs.rs       # Cutoff score calculations
│   │   ├── hsp_culling.rs        # Interval tree HSP culling
│   │   └── translation.rs        # 6-frame translation
│   ├── blastn/            # BLASTN implementation
│   └── common/            # Shared algorithm utilities
├── stats/                 # Karlin-Altschul statistics, E-value, sum statistics
├── sequence/              # 2-bit packed nucleotide sequences
├── utils/                 # BLOSUM62 matrix, genetic codes, SEG masking
└── report/                # Output formatting (outfmt6, outfmt7)
```

## Core Principle: NCBI Parity

- **NCBI C/C++ implementation is the ONLY source of truth**
- Output must match NCBI BLAST+ **bit-perfectly**
- Never simplify algorithms for "readability" if it affects output
- When porting NCBI code, include verbatim C/C++ code as comments:

```rust
// NCBI reference: blast_parameters.c:219-221
// Int4 x_dropoff = (Int4)(sbp->scale_factor * ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
let x_dropoff = (scale_factor * (x_drop_bits * NCBIMATH_LN2 / ungapped_params.lambda).ceil()) as i32;
```

## NCBI Reference Location

- **NCBI BLAST source**: `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/` or `/mnt/c/Users/genom/GitHub/ncbi-blast/`
- Key files:
  - `c++/src/algo/blast/core/aa_ungapped.c` - Extension logic
  - `c++/src/algo/blast/core/link_hsps.c` - Sum-statistics linking
  - `c++/src/algo/blast/core/blast_parameters.c` - Cutoff calculations
  - `c++/src/algo/blast/core/blast_stat.c` - Karlin parameters

## Sequence Encoding

1. **Nucleotides**: 2-bit packed (A=00, C=01, G=10, T=11), 4 bases per byte
2. **Amino Acids**: NCBISTDAA (0-27), sentinel byte = 0 (NULLB)
3. **Scoring**: BLOSUM62 matrix (25×25)
4. **Frames**: 6 reading frames (+1,+2,+3,-1,-2,-3), each with sentinel bytes at boundaries

## Key Constants

| Constant | Value | Description |
|----------|-------|-------------|
| `TWO_HIT_WINDOW` | 40 | Window for two-hit requirement |
| `X_DROP_UNGAPPED_BITS` | 7.0 | X-drop threshold in bits |
| `GAP_TRIGGER_BIT_SCORE` | 22.0 | Gap trigger threshold |
| `CUTOFF_E_TBLASTX` | 1e-300 | Fixed E-value for cutoffs |
| `GAP_SIZE` | 40 | Sum-stats linking gap |
| `OVERLAP_SIZE` | 9 | Sum-stats linking overlap |
| BLOSUM62 ungapped | λ=0.3176, K=0.134 | Karlin parameters |
| BLOSUM62 gapped | λ=0.267, K=0.041 | Karlin parameters |

## Critical Implementation Details

### Coordinate Systems
- **Extension**: Frame-specific buffers, coordinates are frame-relative
- **Linking**: 0-indexed frame-relative coordinates
- **Output**: 1-indexed nucleotide positions

### Length Adjustment Asymmetry (tblastx)
- Query: Apply full `length_adjustment`
- Subject: Apply 1/3 of `length_adjustment`
- Reference: `link_hsps.c:560-571`

### Cutoff Score Three-Stage Capping
```
cutoff = MIN(BLAST_Cutoffs, gap_trigger, cutoff_score_max)
```

### SEG Masking
- Query only (subject never masked in tblastx)
- Extension uses masked sequence (`X=21`)
- Identity calculation uses unmasked sequence

### Chain Member Filtering
- Filter `linked_set && !start_of_chain` during OUTPUT phase, not linking phase
- Reference: `link_hsps.c:1014-1020`

### Subject Frame Sort Order
- Negative frames come FIRST (ascending by frame value)
- Reference: `link_hsps.c:351-357`

### Scale Factor
- tblastx always uses `scale_factor = 1.0`
- Only RPS-BLAST uses `scale_factor > 1.0`

## Debug Environment Variables

```bash
LOSAT_TRACE_HSP="qstart,qend,sstart,send"  # Trace specific HSP
LOSAT_DEBUG_CUTOFFS=1                       # Show cutoff calculations
LOSAT_DEBUG_CHAINING=1                      # Show chaining debug
LOSAT_TIMING=1                              # Print timing breakdown
```

## Performance Patterns

- Use `#[inline]` / `#[inline(always)]` for hot-path functions
- SIMD (AVX2/SSE2) for k-mer scanning and offset pair copying
- `unsafe` only with safety comments for bounds-checked hot paths
- Rayon for parallel processing where algorithm permits

## Naming Conventions

- Match NCBI function names when porting: `s_BlastAaExtendTwoHit` → `extend_hit_two_hit`
- Keep NCBI terminology in comments: "query_offset", "subject_offset", "context", "frame"

## Common Pitfalls

1. Don't simplify NCBI logic for readability if it changes output
2. Don't use different floating-point precision than NCBI
3. Don't skip edge cases (empty sequences, boundary conditions)
4. Always verify bit-perfect output matching after changes

## Known Issues

1. **Long sequences (600kb+)**: Excessive hits (2x NCBI count)
2. **Chain formation differences**: E-value mismatches with short HSPs
