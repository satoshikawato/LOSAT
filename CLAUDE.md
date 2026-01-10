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

# BLASTN (nucleotide vs nucleotide, megablast default)
./target/release/LOSAT blastn -q query.fasta -s subject.fasta -o output.txt

# BLASTN (traditional blastn task)
./target/release/LOSAT blastn --task blastn -q query.fasta -s subject.fasta -o output.txt
```

### TBLASTX Options
```bash
-e, --evalue <EVALUE>              E-value threshold (default: 10.0)
-t, --threshold <THRESHOLD>        Neighbor threshold (default: 13)
-w, --word-size <WORD_SIZE>        Word size (default: 3)
-n, --num-threads <N>              Number of threads (0 = auto)
-v, --verbose                      Verbose logging to stderr
--query-gencode <N>                Genetic code for query (default: 1)
--db-gencode <N>                   Genetic code for subject (default: 1)
--outfmt <FORMAT>                  Output format (0=pairwise, 6=tabular, 7=tabular+headers)
--seg [true|false]                 SEG low-complexity masking (default: true)
--window-size <N>                  Two-hit window size (default: 40, 0=one-hit mode)
--neighbor-map                     Use pre-computed neighbor map (faster scanning)
--ncbi-compat                      Use NCBI-compatible parameters
--culling-limit <N>                HSP culling limit (default: 0=disabled)
--percent-identity <N>             Minimum percent identity filter
--min-hit-length <N>               Minimum alignment length filter
--only-qframe <1|2|3|-1|-2|-3>     Restrict to single query frame
--only-sframe <1|2|3|-1|-2|-3>     Restrict to single subject frame
```

### BLASTN Options
```bash
--task <megablast|blastn>          Algorithm (default: megablast)
-w, --word-size <N>                Word size (default: 28 for megablast)
--reward <N>                       Match reward (default: 1)
--penalty <N>                      Mismatch penalty (default: -2)
--gap-open <N>                     Gap open penalty (default: 0=auto)
--gap-extend <N>                   Gap extend penalty (default: 0=auto)
--dust                             Enable DUST low-complexity filter
--scan-step <N>                    Scan stride (default: auto)
--hitlist-size <N>                 Maximum hits to save (default: 500)
--max-hsps-per-subject <N>         Max HSPs per subject (0=unlimited)
```

## Integration Tests

```bash
# Run comparison tests against NCBI BLAST (from LOSAT/tests/ directory)
cd LOSAT/tests && bash run_comparison.sh

# Run all tests including plots
cd LOSAT/tests && bash run_all_tests.sh
```

## Project Overview

LOSAT (LOcal Sequence Alignment Tool) is a Rust reimplementation of NCBI BLAST targeting **bit-perfect parity** with NCBI BLAST+ output.

Supported algorithms:
- **TBLASTX**: Translated nucleotide vs translated nucleotide (primary focus)
- **BLASTN**: Nucleotide vs nucleotide (megablast/blastn tasks)

## Architecture

```
LOSAT/src/
├── algorithm/
│   ├── tblastx/                   # TBLASTX implementation (primary focus)
│   │   ├── blast_engine/          # Main execution orchestration
│   │   │   ├── run_impl.rs        # Standard run() implementation
│   │   │   └── run_neighbor_map.rs # Optimized neighbor-map mode
│   │   ├── extension/             # Hit extension
│   │   │   ├── two_hit.rs         # Two-hit extension logic
│   │   │   ├── ungapped.rs        # Ungapped extension
│   │   │   └── gapped.rs          # Gapped extension
│   │   ├── lookup/                # Lookup table management
│   │   │   ├── backbone.rs        # Core lookup table
│   │   │   └── neighbor_map.rs    # Pre-computed neighbor relationships
│   │   ├── scan/                  # Sequence scanning
│   │   │   ├── offset_pairs.rs    # Offset pair generation
│   │   │   └── simd_helpers.rs    # AVX2/SSE2 SIMD optimizations
│   │   ├── sum_stats_linking/     # HSP chaining and E-value
│   │   │   ├── linking.rs         # Sum-statistics linking algorithm
│   │   │   ├── cutoffs.rs         # Cutoff calculations
│   │   │   └── params.rs          # Linking parameters
│   │   ├── filtering/             # HSP filtering
│   │   │   └── purge_endpoints.rs # Endpoint-based filtering
│   │   ├── constants.rs           # TBLASTX-specific constants
│   │   ├── ncbi_cutoffs.rs        # NCBI cutoff score calculations
│   │   ├── hsp_culling.rs         # Interval tree HSP culling
│   │   ├── chaining.rs            # HSP chaining logic
│   │   ├── translation.rs         # 6-frame translation
│   │   └── tracing.rs             # Debug tracing infrastructure
│   ├── blastn/                    # BLASTN implementation
│   │   ├── blast_engine/          # Main execution
│   │   │   └── run.rs             # blastn/megablast run logic
│   │   ├── alignment/             # Gapped alignment
│   │   │   ├── gapped.rs          # Semi-global gapped alignment
│   │   │   ├── greedy.rs          # Greedy alignment (megablast)
│   │   │   └── statistics.rs      # Alignment statistics
│   │   ├── filtering/             # HSP filtering
│   │   │   ├── purge_endpoints.rs # Endpoint filtering
│   │   │   └── subject_best_hit.rs # Subject-best-hit filtering
│   │   ├── constants.rs           # BLASTN-specific constants
│   │   ├── ncbi_cutoffs.rs        # NCBI cutoff calculations
│   │   ├── lookup.rs              # Nucleotide lookup table
│   │   ├── extension.rs           # Ungapped extension
│   │   └── coordination.rs        # Multi-threading coordination
│   └── common/                    # Shared algorithm utilities
│       ├── chaining.rs            # Shared chaining logic
│       ├── diagnostics.rs         # Debug diagnostics
│       └── evalue.rs              # E-value calculations
├── core/                          # Core BLAST primitives (NCBI ports)
│   ├── blast_stat/                # Karlin-Altschul statistics
│   │   ├── karlin_params.rs       # Karlin parameter lookup
│   │   ├── length_adjustment.rs   # Effective length calculation
│   │   ├── search_space.rs        # Search space computation
│   │   └── sum_statistics.rs      # Sum statistics
│   ├── blast_filter.rs            # Low-complexity filtering
│   ├── blast_seg.rs               # SEG algorithm
│   ├── blast_encoding.rs          # Sequence encoding
│   ├── blast_util.rs              # Utility functions
│   └── gencode_singleton.rs       # Genetic code tables
├── stats/                         # Statistical calculations
├── sequence/                      # 2-bit packed nucleotide sequences
├── utils/                         # BLOSUM62 matrix, genetic codes, SEG masking
├── report/                        # Output formatting (outfmt0, outfmt6, outfmt7)
├── post/                          # Post-processing (chain, filter)
├── blastinput/                    # CLI argument parsing
├── api/                           # API interfaces
└── format/                        # Output formatting utilities
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
  - `c++/src/algo/blast/core/na_ungapped.c` - Nucleotide extension
  - `c++/src/algo/blast/core/link_hsps.c` - Sum-statistics linking
  - `c++/src/algo/blast/core/blast_parameters.c` - Cutoff calculations
  - `c++/src/algo/blast/core/blast_stat.c` - Karlin parameters
  - `c++/src/algo/blast/core/greedy_align.c` - Greedy alignment

## Sequence Encoding

1. **Nucleotides**: 2-bit packed (A=00, C=01, G=10, T=11), 4 bases per byte
2. **Amino Acids**: NCBISTDAA (0-27), sentinel byte = 0 (NULLB)
3. **Scoring**: BLOSUM62 matrix (25x25)
4. **Frames**: 6 reading frames (+1,+2,+3,-1,-2,-3), each with sentinel bytes at boundaries

## Key Constants

### TBLASTX Constants
| Constant | Value | Description |
|----------|-------|-------------|
| `TWO_HIT_WINDOW` | 40 | Window for two-hit requirement |
| `X_DROP_UNGAPPED_BITS` | 7.0 | X-drop threshold in bits |
| `X_DROP_UNGAPPED` | 16 | Pre-calculated X-drop raw score |
| `GAP_TRIGGER_BIT_SCORE` | 22.0 | Gap trigger threshold |
| `CUTOFF_E_TBLASTX` | 1e-300 | Fixed E-value for cutoffs |
| `MIN_UNGAPPED_SCORE` | 14 | Minimum ungapped score |
| `SENTINEL_BYTE` | 0 | Sequence boundary marker |
| `SENTINEL_PENALTY` | -4 | Penalty for sentinel hits |
| `STOP_CODON` | 24 | Stop codon encoding |
| BLOSUM62 ungapped | lambda=0.3176, K=0.134 | Karlin parameters |
| BLOSUM62 gapped | lambda=0.267, K=0.041 | Karlin parameters |

### BLASTN Constants
| Constant | Value | Description |
|----------|-------|-------------|
| `X_DROP_UNGAPPED` | 20 | Ungapped X-drop |
| `X_DROP_GAPPED_NUCL` | 30 | Gapped X-drop (blastn) |
| `X_DROP_GAPPED_GREEDY` | 25 | Gapped X-drop (megablast) |
| `X_DROP_GAPPED_FINAL` | 100 | Final traceback X-drop |
| `TWO_HIT_WINDOW` | 0 | One-hit mode (NCBI default) |
| `SCAN_RANGE_BLASTN` | 4 | Off-diagonal scan range |
| `GREEDY_MAX_COST` | 1000 | Greedy alignment max cost |

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
# HSP Tracing (trace specific HSP through pipeline)
LOSAT_TRACE_HSP="qstart,qend,sstart,send"  # Example: "111880,111860,163496,163476"

# Debug Flags
LOSAT_DEBUG_CUTOFFS=1                       # Show cutoff calculations
LOSAT_DEBUG_CHAINING=1                      # Show chaining debug
LOSAT_DEBUG_EXTENSION=1                     # Show extension debug (tblastx)
LOSAT_DEBUG_BLASTN=1                        # Show BLASTN hit loss diagnostics
LOSAT_DEBUG_COORDS=1                        # Show coordinate transformations (blastn)

# Performance
LOSAT_TIMING=1                              # Print timing breakdown
LOSAT_NO_SIMD=1                             # Force scalar processing (disable AVX2)

# Diagnostics
LOSAT_DIAGNOSTICS=1                         # Enable general diagnostics
LOSAT_STARTUP_TRACE=1                       # Trace startup
```

## Performance Patterns

- Use `#[inline]` / `#[inline(always)]` for hot-path functions
- SIMD (AVX2/SSE2) for k-mer scanning and offset pair copying
- `unsafe` only with safety comments for bounds-checked hot paths
- Rayon for parallel processing where algorithm permits

## Naming Conventions

- Match NCBI function names when porting: `s_BlastAaExtendTwoHit` -> `extend_hit_two_hit`
- Keep NCBI terminology in comments: "query_offset", "subject_offset", "context", "frame"

## Common Pitfalls

1. Don't simplify NCBI logic for readability if it changes output
2. Don't use different floating-point precision than NCBI
3. Don't skip edge cases (empty sequences, boundary conditions)
4. Always verify bit-perfect output matching after changes

## Known Issues

1. **Long sequences (600kb+)**: Excessive hits (2x NCBI count) - under investigation
2. **Chain formation differences**: E-value mismatches with short HSPs in some edge cases
