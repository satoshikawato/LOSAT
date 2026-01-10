# LOSAT BLASTN Implementation Status

**Last Updated**: 2026-01-10
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

**Root Cause**: Issues in gapped extension (gapped.rs) - original analysis was incorrect.

**NCBI Reference**: `blast_gapalign.c:933-952` (band contraction)

**Files Modified**:
- `src/algorithm/blastn/alignment/gapped.rs`

### 7. Alignment Length Off-by-Two Fix (Session 4) - FIXED
**Problem**: LOSAT produced 21-base alignments while NCBI produced 19-base alignments
```
NCBI:  180 hits at length 19 (100% identity)
LOSAT: 186 hits at length 21 (90.476% identity = 19/21 matches)
```

**Root Cause**: In LOSAT's DP loop, `score_val` at loop index `b_index=k` contains the score from the PREVIOUS iteration (`next_score` at `b_index=k-1`), which represents position `k-1`, not `k`. When we set `b_offset=k`, it's 1 higher than the actual consumed position.

**Fix**: Three changes to match NCBI behavior:
1. Removed post-loop score comparison (didn't exist in NCBI, was over-extending alignments)
2. Changed band contraction from `last_b_index + 2` back to `last_b_index + 1` (matches NCBI)
3. Changed `s_consumed = b_offset` (not `b_offset + 1`) to account for the offset shift

**Result**: Alignment lengths now match NCBI exactly (19 bases with 100% identity)
```
AFTER:  8,542 hits with correct 19-base alignments (69.2% coverage)
```

**Files Modified**:
- `src/algorithm/blastn/alignment/gapped.rs` (s_consumed calculation, band contraction, removed post-loop comparison)

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
**Symptom**: NZ_CP006932 self-comparison shows 72% of NCBI hits (8,843 vs 12,340)

**Status**: Partially fixed - endpoint purging improved from 8,489 to 8,843 hits.

**Investigation (2026-01-10)**:
- Tested lut_word_length=11 (scan_step=1) but hit count remained at ~8,500
- Alignment length distribution analysis shows missing hits are in 30-99 bp range:
  - NCBI 30-49: 4,245 | LOSAT: 2,570 (-1,675)
  - NCBI 50-99: 3,146 | LOSAT: 1,454 (-1,692)
  - Short alignments (<30) are similar between NCBI and LOSAT
- Debug output shows 13,602 hits before endpoint purge → 8,843 after
- Endpoint purging removes 35% of raw hits due to duplicate endpoints

**Fixes Applied (2026-01-10)**:
- Fixed endpoint purge to match NCBI's `Blast_HSPListPurgeHSPsWithCommonEndpoints`:
  1. Added per-subject grouping (NCBI's BlastHSPList is per-subject)
  2. Added query_id to criteria (emulates NCBI's context)
  3. Changed to canonical coordinates (subject.offset = min, subject.end = max)
  4. Removed strand from sort key (matches NCBI's s_QueryOffsetCompareHSPs)
  5. Added tie-breaker fields to match NCBI exactly (blast_hits.c:2310-2318, 2376-2384):
     - Pass 1: Added q_end DESC, s_end DESC tie-breakers
     - Pass 2: Added q_start DESC, s_offset DESC tie-breakers
- Improvement: 8,489 → 8,843 hits (4% improvement)
- Note: Tie-breaker addition didn't change hit count (same 8,843) but ensures NCBI parity

**NCBI Traceback Trimming Discovery**:
NCBI's blastn has a sophisticated 2-pass endpoint purge in blast_traceback.c:
1. First pass (line 638) with `purge=FALSE`: Instead of deleting duplicates, NCBI TRIMS them
   - Uses `s_CutOffGapEditScript()` to extract non-overlapping fragments
   - Requires traceback (GapEditScript) information
2. Re-evaluate trimmed fragments (lines 647-665): `Blast_HSPReevaluateWithAmbiguitiesGapped()`
3. Second pass (line 668) with `purge=TRUE`: Final cleanup of remaining duplicates

LOSAT only implements the `purge=TRUE` behavior (delete duplicates entirely).
Implementing the trimming would require traceback infrastructure not currently available.

**Remaining Gap Analysis**:
- LOSAT generates 13,602 raw hits, NCBI has 12,340 final hits
- LOSAT purges 35% (4,759 hits) due to duplicate endpoints
- Root cause: LOSAT generates more HSPs with same endpoints (before purge)
- Potential sources: seed finding, diagonal tracking, extension overlap
- The trimming feature above may also contribute to NCBI retaining more HSPs

**Files Modified**:
- `src/algorithm/blastn/filtering/purge_endpoints.rs`

**Priority**: Medium (remaining gap likely in seed/extension phase)

### Issue 2: Extra Hits in Some Cases
**Symptom**: Some datasets show 105-116% of NCBI hits

**Status**: Subject Best Hit filtering is NOT enabled by default in NCBI.

**Investigation (2026-01-10)**:
- Implemented subject_best_hit filter (blast_hits.c:2537-2606)
- Discovered it's OPTIONAL in NCBI - only applied when `-subject_besthit` is specified
- Filter is available but disabled by default to match NCBI behavior
- Extra hits distribution (PeseMJNV.PemoMJNVB):
  - NCBI <30: 3,046 | LOSAT: 4,016 (+970)
  - NCBI 30-49: 4,638 | LOSAT: 5,455 (+817)
- Extra hits are predominantly short alignments

**Root Cause**: Different seed finding or extension behavior, possibly:
1. More seeds found due to different diagonal tracking
2. Different X-drop termination
3. Different cutoff score calculations

**Priority**: Medium

## Test Results Summary

| Dataset | Algorithm | LOSAT | NCBI | Coverage | Status |
|---------|-----------|-------|------|----------|--------|
| EDL933.Sakai | megablast | 5,634 | 5,718 | 98.5% | ✓ Excellent |
| Sakai.MG1655 | megablast | 6,434 | 6,476 | 99.3% | ✓ Excellent |
| NZ_CP006932 self | blastn | 8,843 | 12,340 | 71.7% | ⚠ Partially fixed, see Issue 1 |
| MjeNMV.MelaMJNV | blastn | 2,588 | 2,668 | 97.0% | ✓ Good |
| SiNMV.ChdeNMV | blastn | 4,192 | 4,367 | 95.9% | ✓ Good |
| MelaMJNV.PemoMJNVA | blastn | 3,088 | 2,729 | 113.1% | ⚠ Extra short hits |
| PemoMJNVA.PeseMJNV | blastn | 3,094 | 2,940 | 105.2% | ⚠ Extra short hits |
| PeseMJNV.PemoMJNVB | blastn | 13,612 | 11,668 | 116.6% | ⚠ Extra short hits |

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
**Investigation complete**: Scan_step is not the issue. Missing hits are 30-99 bp alignments.

Next actions:
1. Compare gapped extension output between NCBI and LOSAT for same seeds
2. Check if endpoint purging is removing valid HSPs (13,520 → 8,489 seems aggressive)
3. Investigate why LOSAT produces fewer medium-length alignments

### Priority 2: Investigate Extra Hits
**Investigation complete**: Subject Best Hit filtering is optional in NCBI.

Next actions:
1. Add `-subject_besthit` CLI option to enable optional filtering
2. Compare X-drop termination thresholds
3. Compare cutoff score calculations for short alignments

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
