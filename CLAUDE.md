# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## MANDATORY COMPLIANCE REQUIREMENTS

**These rules are ABSOLUTE and MUST be followed without exception. There is NO room for interpretation.**

### 1. NCBI BLAST is the ONLY Source of Truth
- The NCBI C/C++ implementation is the **sole authoritative reference**
- **DO NOT ASSUME OR GUESS** - always refer to actual NCBI source code
- If you cannot find the corresponding NCBI code, the feature **MUST NOT exist**
- **Faithfully transpile NCBI BLAST** - algorithms must be exact ports, not reinterpretations

### 2. Bit-Perfect Output Parity is Required
- Output must match NCBI BLAST+ **exactly** (1ビットの狂いもなく一致)
- Never simplify algorithms for "readability" if it affects output
- Never use different floating-point precision than NCBI
- Every nuanced difference must be identified and fixed

### 3. NCBI Code Comments are MANDATORY
- **Every code modification MUST include the corresponding NCBI C/C++ code as comments**
- Include file path and line numbers
- **If you cannot add NCBI reference comments, the code MUST NOT be written**
- If the corresponding NCBI code cannot be found, the feature **does not exist and must be eliminated**

```rust
// NCBI reference: blast_parameters.c:219-221
// Int4 x_dropoff = (Int4)(sbp->scale_factor * ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
let x_dropoff = (scale_factor * (x_drop_bits * NCBIMATH_LN2 / ungapped_params.lambda).ceil()) as i32;
```

### 4. No Unauthorized Features - IMMEDIATE ELIMINATION
- **NEVER introduce features/functionalities that do not exist in the NCBI codebase**
- If a feature is found that has no NCBI equivalent, **DELETE IT IMMEDIATELY**
- No "improvements" or "optimizations" that change behavior
- No creative additions - only faithful transpilation

### 5. Exhaustive Difference Resolution
- Identify **ALL nuanced differences** between LOSAT and NCBI BLAST
- **Exhaustively add/fix/delete every single one of them**
- Verify bit-perfect matching after every change
- Leave no discrepancy unaddressed

### 6. No Assumptions or Guessing
- **DO NOT ASSUME** when writing code
- **DO NOT GUESS** behavior or implementation details
- **REFER TO THE NCBI CODE** and transpile faithfully
- When in doubt, read the NCBI source - never speculate

### 7. Correct Timing, Order, and Context (Common Mistake)
- Algorithms must be called at the **exact same timing and order** as NCBI
- Input/output data (context) must match NCBI **exactly**
- **Wrong timing or order = wrong output** even if the algorithm itself is correct
- Verify: When is the function called? In what order? What data does it receive? What does it return?
- Trace the call hierarchy in NCBI to understand the correct placement and sequence

### Violation of these rules is unacceptable. Non-compliance will result in incorrect output.

---

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

### Entry Point
- **CLI**: `src/main.rs` - Uses clap subcommands, dispatches to `blastn::run()` or `tblastx::run()`
- **TBLASTX engine**: `src/algorithm/tblastx/blast_engine/run_impl.rs`
- **BLASTN engine**: `src/algorithm/blastn/blast_engine/run.rs`

## Core Principle: NCBI Parity

**See [MANDATORY COMPLIANCE REQUIREMENTS](#mandatory-compliance-requirements) at the top of this document.**

Summary:
1. **NCBI C/C++ is the ONLY truth** - Faithfully transpile, no reinterpretations
2. **Bit-perfect output** - 1ビットの狂いもなく一致
3. **NCBI comments required** - No NCBI reference = code must not exist
4. **No unauthorized features** - If not in NCBI, DELETE IMMEDIATELY
5. **Fix ALL differences** - Exhaustively add/fix/delete every discrepancy
6. **No assumptions** - DO NOT ASSUME, DO NOT GUESS, REFER TO NCBI CODE
7. **Correct timing, order, and context** - Same call timing, order, and input/output as NCBI

## NCBI Reference Location

- **NCBI BLAST source**: Machine-dependent path (e.g., `~/GitHub/ncbi-blast/` or similar)
- Key files:
  - `c++/src/algo/blast/core/aa_ungapped.c` - Extension logic
  - `c++/src/algo/blast/core/na_ungapped.c` - Nucleotide extension
  - `c++/src/algo/blast/core/link_hsps.c` - Sum-statistics linking
  - `c++/src/algo/blast/core/blast_parameters.c` - Cutoff calculations
  - `c++/src/algo/blast/core/blast_query_info.c` - Context management
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

### DiagStruct Initialization
- **NCBI**: `last_hit = -window` (initialized to -40)
- **LOSAT**: `last_hit = 0`
- Both result in `diff >= window` on first hit (record only, no extension)
- This difference does NOT affect output
- Reference: `blast_extend.c:103`

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

## Common Pitfalls (VIOLATIONS)

**Any of these will break NCBI parity:**

1. **Simplifying NCBI logic** - Never change logic for "readability"
2. **Different floating-point precision** - Use exact same precision as NCBI
3. **Skipping edge cases** - Handle all boundary conditions exactly as NCBI
4. **Missing verification** - Always verify bit-perfect output after changes
5. **Assuming behavior** - Always read NCBI source, never guess
6. **Adding features** - Never add anything not in NCBI
7. **Missing NCBI comments** - Every change must cite NCBI source code
8. **Wrong timing, order, or context** - Algorithm called at wrong time, in wrong order, or with wrong input/output (COMMON MISTAKE)

## Known Issues

1. **Long sequences (600kb+)**: Excessive hits (2x NCBI count) - under investigation
   - **Status**: Exhaustive investigation completed (Session 12, 2026-01-11). Root cause still unknown.
   - **Verified CORRECT (not the cause)**:
     - DiagStruct initialization, two-hit state machine, extension algorithm
     - Seed generation, lookup table construction, frame iteration
     - Cutoff calculation, context boundary check, E-value calculation
     - Scan resumption, diagonal overflow, diag_offset increment
     - All two-hit detection logic matches NCBI `aa_ungapped.c:518-606` exactly
   - **Evidence**: 2x excess is in NUMBER OF EXTENSIONS triggered (31.8M vs ~16M expected)
   - **Score distribution**: 73% of excess hits are in 22-30 bit score range
   - **Hypothesis**: Some hidden NCBI optimization or subtle scale-dependent behavior not yet identified
   - See: [TBLASTX_2X_HIT_INVESTIGATION.md](TBLASTX_2X_HIT_INVESTIGATION.md) and [FIX_2X_HITS_PLAN.md](FIX_2X_HITS_PLAN.md)

2. **Chain formation differences**: E-value mismatches with short HSPs in some edge cases

## Recently Fixed

1. **Coordinate output off-by-one** (Session 11, 2026-01-11): All TBLASTX output coordinates were -1 compared to NCBI
   - **Root Cause**: NCBI's C++ output layer adds +1 for 1-indexed output; LOSAT was missing this
   - **Fix**: Modified `convert_coords()` in `extension/mod.rs` to incorporate +1 adjustment
   - **Verified**: Coordinates now match NCBI exactly
