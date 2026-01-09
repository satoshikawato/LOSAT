# LOSAT BLASTN Implementation Status

**Last Updated**: 2026-01-09
**Branch**: fix2

## Overview

This document summarizes the current state of the BLASTN implementation in LOSAT, including known issues, completed fixes, and next steps.

## What's Working

### Megablast (default task)
- Core algorithm functional
- Hit coverage: **98-99%** compared to NCBI BLAST
- Performance: ~3x slower than NCBI (acceptable for now)

### Blastn (task blastn)
- Core DP-based gapped alignment implemented
- Follows NCBI's `Blast_SemiGappedAlign` algorithm
- Hit coverage: 70-115% depending on dataset (needs investigation)

## Completed Fixes

### 1. Diagonal Tracking (Session 1)
**Problem**: blastn was getting stuck on self-comparison (infinite loop)
**Root Cause**: Diagonal tracking used seed start position instead of ungapped extension end
**Fix**: Updated `last_hit` tracking to use ungapped extension endpoint
**File**: `src/algorithm/blastn/blast_extend.rs`

### 2. Score-Only DP Implementation (Session 1-2)
**Problem**: DP gapped alignment was 6x slower than NCBI
**Root Cause**: LOSAT tracked alignment statistics inline (40 bytes/cell vs NCBI's 8 bytes)
**Fix**: Reimplemented `extend_gapped_one_direction` to match NCBI's `Blast_SemiGappedAlign`
- Uses minimal `BlastGapDP` struct (8 bytes: `best` + `best_gap`)
- Statistics estimated from final score/positions
- Dynamic memory reallocation without fixed limits

**File**: `src/algorithm/blastn/alignment/gapped.rs:209-485`

**NCBI Reference**: `blast_gapalign.c:735-962`

### 3. Left Extension Optimization (Session 2)
**Problem**: Left extension still used old 40-byte cell implementation
**Fix**: Updated `extend_gapped_one_direction_ex` to use same minimal approach
**Result**: blastn reduced from 42s → 26s (38% improvement)

**File**: `src/algorithm/blastn/alignment/gapped.rs:498-731`

### 4. Coordinate Off-by-One Fix (Session 3) - FIXED
**Problem**: Top hit started at position 2 instead of 1
```
BEFORE: 2-657100 (length: 657099)
AFTER:  1-657101 (length: 657101) ✓
```

**Root Cause**: Three issues in gapped extension (gapped.rs):
1. Band contraction too aggressive: `b_size = last_b_index + 1` → `b_size = last_b_index + 2`
2. Missing final score comparison after inner loop
3. Incorrect s_consumed indexing: `b_offset` → `b_offset + 1`

**NCBI Reference**: `blast_gapalign.c:933-952` (band contraction)

**Files Modified**:
- `src/algorithm/blastn/alignment/gapped.rs` (lines 435-453, 500-505, 715-742, 778-783)

### 5. Cutoff Score Calculation Fix (Session 3) - FIXED
**Problem**: Using same Karlin params for ungapped and gapped calculations
**Fix**: Now correctly uses:
- **UNGAPPED params (kbp_std)** for gap_trigger calculation
- **GAPPED params (kbp_gap)** for cutoff_score_max calculation

**NCBI Reference**: `blast_parameters.c:340-344`

**Files Modified**:
- `src/algorithm/blastn/blast_engine/run.rs` (lines 62-89, 330-349)

### 6. Removed Non-NCBI Containment Filter (Session 3)
**Problem**: Post-processing containment filter removed too many valid HSPs
**Fix**: Removed entirely - NCBI's containment filtering (`BlastIntervalTreeContainsHSP`) is applied DURING alignment, not as post-processing

**NCBI Reference**: `blast_gapalign.c:3918`, `blast_itree.c:810-847`

## Known Issues

### Issue 1: Missing Hits in Self-Comparison
**Symptom**: NZ_CP006932 self-comparison shows only 70% of NCBI hits (8,609 vs 12,340)

**Possible Causes**:
1. Differences in seed finding logic
2. Extension algorithm differences
3. Scoring calculation differences

**Priority**: High

### Issue 2: Extra Hits in Some Cases
**Symptom**: Some datasets show 105-115% of NCBI hits

**Possible Causes**:
1. Missing HSP filtering step in NCBI that LOSAT doesn't implement
2. Different scoring thresholds

**Priority**: Medium

## Test Results Summary

| Dataset | Algorithm | LOSAT | NCBI | Coverage | Status |
|---------|-----------|-------|------|----------|--------|
| EDL933.Sakai | megablast | 5,634 | 5,718 | 98.5% | ✓ Excellent |
| Sakai.MG1655 | megablast | 6,434 | 6,476 | 99.3% | ✓ Excellent |
| NZ_CP006932 self | blastn | 8,609 | 12,340 | 69.7% | ✗ Investigate |
| MjeNMV.MelaMJNV | blastn | 2,588 | 2,668 | 97.0% | ✓ Good |
| SiNMV.ChdeNMV | blastn | 4,192 | 4,367 | 95.9% | ✓ Good |
| MelaMJNV.PemoMJNVA | blastn | 3,088 | 2,729 | 113.1% | ⚠ Extra hits |
| PemoMJNVA.PeseMJNV | blastn | 3,094 | 2,940 | 105.2% | ⚠ Extra hits |
| PeseMJNV.PemoMJNVB | blastn | 13,454 | 11,668 | 115.3% | ⚠ Extra hits |

## Key Files

### Core BLASTN Implementation
```
src/algorithm/blastn/
├── alignment/
│   ├── gapped.rs         # DP gapped extension (Blast_SemiGappedAlign)
│   ├── greedy.rs         # Greedy alignment for megablast
│   └── statistics.rs     # Alignment statistics
├── blast_engine/
│   ├── mod.rs            # Engine coordination, HSP filtering
│   └── run.rs            # Main execution loop
├── blast_extend.rs       # Ungapped extension, diagonal tracking
├── extension.rs          # Two-hit extension logic
├── lookup.rs             # Word lookup table
├── ncbi_cutoffs.rs       # Cutoff score calculations
└── filtering/
    └── purge_endpoints.rs # HSP endpoint purging
```

### NCBI Reference Files
```
/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/
├── blast_gapalign.c      # Blast_SemiGappedAlign (lines 735-962)
├── blast_parameters.c    # Cutoff calculations (lines 340-344)
├── blast_hits.c          # HSP filtering (lines 2455-2535)
├── na_ungapped.c         # Nucleotide ungapped extension
└── blast_itree.c         # Interval tree containment check
```

## Next Steps

### Priority 1: Investigate Missing Hits (NZ_CP006932)
1. Compare seed detection between NCBI and LOSAT
2. Add diagnostic logging at each stage
3. Identify where hits are being filtered out

### Priority 2: Investigate Extra Hits
1. Compare HSP scoring between NCBI and LOSAT
2. Check for missing filtering steps in NCBI

### Priority 3: Performance Optimization (Later)
1. Profile to identify bottlenecks
2. Consider SIMD for DP inner loop
3. Optimize memory allocation patterns

## Debug Commands

```bash
# Run blastn with timing
cd LOSAT && time ./target/release/LOSAT blastn -q tests/fasta/NZ_CP006932.fasta -s tests/fasta/NZ_CP006932.fasta -o /tmp/test.out --task blastn -n 1

# Compare hit counts
wc -l tests/losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.out
grep -v "^#" tests/blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out | wc -l

# Compare top hits
head -5 tests/losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.out
grep -v "^#" tests/blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out | head -5

# Run with coordinate debug
LOSAT_DEBUG_COORDS=1 ./target/release/LOSAT blastn -q tests/fasta/AP027280.fasta -s tests/fasta/AP027280.fasta -o /tmp/test.out --task blastn -n 1 2>&1 | head -20
```

## NCBI Code References

### Key Functions Implemented
1. `Blast_SemiGappedAlign` - Score-only DP ✓
2. `gap_trigger_raw_score` - Gap trigger calculation ✓
3. `cutoff_score_max_from_evalue` - Cutoff from E-value ✓
4. `Blast_HSPListPurgeHSPsWithCommonEndpoints` - Endpoint purging ✓

### Important Constants
```c
// From blast_options.h
#define BLAST_GAP_OPEN_NUCL 5
#define BLAST_GAP_EXTN_NUCL 2
#define BLAST_WORD_SIZE_NUCL 28  // megablast
#define BLAST_WORD_SIZE_NUCL_BLASTN 11  // blastn
#define BLAST_GAP_TRIGGER_NUCL 27.0  // Gap trigger bit score
#define CUTOFF_E_BLASTN 0.05  // Ungapped cutoff E-value
```
