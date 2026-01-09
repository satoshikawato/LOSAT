# LOSAT BLASTN Implementation Status

**Last Updated**: 2026-01-09
**Branch**: fix2

## Overview

This document summarizes the current state of the BLASTN implementation in LOSAT, including known issues, completed fixes, and next steps.

## What's Working

### Megablast (default task)
- Core algorithm functional
- Hit coverage: 89-99% compared to NCBI BLAST
- Performance: ~3x slower than NCBI (acceptable for now)

### Blastn (task blastn)
- Core DP-based gapped alignment implemented
- Follows NCBI's `Blast_SemiGappedAlign` algorithm
- Hit coverage varies: 75-111% depending on dataset

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

## Known Issues

### Issue 1: Coordinate Off-by-One
**Symptom**: Top hit starts at position 2 instead of 1
```
NCBI:  1	657101	1	657101  (length: 657101)
LOSAT: 2	657100	2	657100  (length: 657099)
```

**Likely Cause**: Off-by-one in coordinate conversion during gapped alignment or output formatting

**Files to Investigate**:
- `src/algorithm/blastn/alignment/gapped.rs` - coordinate handling in DP
- `src/algorithm/blastn/blast_engine/run.rs` - HSP coordinate conversion
- `src/report/outfmt6.rs` - output coordinate formatting

**Priority**: Medium (affects output accuracy but not significantly)

### Issue 2: Missing Hits in Self-Comparison
**Symptom**: NZ_CP006932 self-comparison shows only 75% of NCBI hits (9,283 vs 12,340)

**Possible Causes**:
1. HSP filtering differences
2. Two-hit extension threshold differences
3. Ungapped extension X-drop differences
4. Score cutoff differences

**Files to Investigate**:
- `src/algorithm/blastn/blast_extend.rs` - ungapped extension
- `src/algorithm/blastn/extension.rs` - two-hit logic
- `src/algorithm/blastn/filtering/` - HSP filtering

**Priority**: High (significant hit loss)

### Issue 3: Extra Hits in Some Cases
**Symptom**: MelaMJNV.PemoMJNVA shows 111% of NCBI hits (3,049 vs 2,729)

**Possible Causes**:
1. Missing HSP filtering/culling step
2. Different E-value cutoff application
3. Missing containment filtering

**Priority**: Medium (extra hits are less problematic than missing hits)

### Issue 4: Performance Gap
**Current**: blastn ~4x slower than NCBI (26s vs 6.7s for NZ_CP006932 self)
**Current**: megablast ~3x slower than NCBI (0.86s vs 0.28s)

**Possible Optimizations**:
1. SIMD in DP inner loop
2. Better memory allocation strategy
3. Reduce unnecessary gapped extensions
4. Profile to find actual bottlenecks

**Priority**: Low (functionality first, optimize later)

## Test Results Summary

| Dataset | Algorithm | LOSAT | NCBI | Coverage | Status |
|---------|-----------|-------|------|----------|--------|
| EDL933.Sakai | megablast | 5,611 | 5,718 | 98.1% | ✓ Good |
| Sakai.MG1655 | megablast | 6,422 | 6,476 | 99.2% | ✓ Good |
| NZ_CP006932 self | megablast | 407 | 454 | 89.6% | ⚠ Acceptable |
| NZ_CP006932 self | blastn | 9,283 | 12,340 | 75.2% | ✗ Investigate |
| MjeNMV.MelaMJNV | blastn | 2,628 | 2,668 | 98.5% | ✓ Good |
| MelaMJNV.PemoMJNVA | blastn | 3,049 | 2,729 | 111.7% | ⚠ Extra hits |

## Key Files

### Core BLASTN Implementation
```
src/algorithm/blastn/
├── alignment/
│   ├── gapped.rs         # DP gapped extension (Blast_SemiGappedAlign)
│   ├── greedy.rs         # Greedy alignment for megablast
│   └── statistics.rs     # Alignment statistics
├── blast_engine/
│   ├── mod.rs            # Engine coordination
│   └── run.rs            # Main execution loop
├── blast_extend.rs       # Ungapped extension, diagonal tracking
├── extension.rs          # Two-hit extension logic
├── lookup.rs             # Word lookup table
└── filtering/            # HSP filtering
```

### NCBI Reference Files
```
/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/
├── blast_gapalign.c      # Blast_SemiGappedAlign (lines 735-962)
├── na_ungapped.c         # Nucleotide ungapped extension
├── blast_extend.c        # Extension algorithms
└── blast_parameters.c    # Score parameters
```

## Next Steps

### Priority 1: Investigate Missing Hits
1. Add diagnostic logging to count hits at each stage
2. Compare NCBI vs LOSAT at:
   - Seed detection
   - Ungapped extension
   - Two-hit trigger
   - Gapped extension
   - Final filtering
3. Identify where hits are being lost

### Priority 2: Fix Coordinate Off-by-One
1. Trace coordinate flow from DP to output
2. Compare with NCBI coordinate handling
3. Verify 0-indexed vs 1-indexed conversions

### Priority 3: Investigate Extra Hits
1. Check if LOSAT is missing a filtering step
2. Compare E-value calculations
3. Check HSP containment/overlap filtering

### Priority 4: Performance Optimization (Later)
1. Profile to identify bottlenecks
2. Consider SIMD for DP inner loop
3. Optimize memory allocation patterns

## Debug Commands

```bash
# Run blastn with timing
cd LOSAT && time ./target/release/LOSAT blastn -q tests/fasta/NZ_CP006932.fasta -s tests/fasta/NZ_CP006932.fasta -o /tmp/test.out --task blastn -n 1

# Compare hit counts
grep -v "^#" tests/blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out | wc -l
grep -v "^#" tests/losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.out | wc -l

# Compare top hits
head -5 tests/blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out
head -5 tests/losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.out

# Run with NCBI for comparison
blastn -task blastn -query tests/fasta/NZ_CP006932.fasta -subject tests/fasta/NZ_CP006932.fasta -out /tmp/ncbi.out -outfmt 7
```

## NCBI Code References

### Key Functions to Match
1. `Blast_SemiGappedAlign` - Score-only DP (implemented ✓)
2. `s_BlastNaExtend` - Nucleotide ungapped extension
3. `BLAST_DiagUpdate` - Diagonal tracking
4. `s_NuclUngappedExtend` - Core extension logic
5. `BlastHitSavingParametersNew` - Hit filtering parameters

### Important Constants
```c
// From blast_options.h
#define BLAST_GAP_OPEN_NUCL 5
#define BLAST_GAP_EXTN_NUCL 2
#define BLAST_WORD_SIZE_NUCL 28  // megablast
#define BLAST_WORD_SIZE_NUCL_BLASTN 11  // blastn
```
