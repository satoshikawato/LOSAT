# Two-Stage Lookup Table Implementation for BLASTN/Megablast

## Summary

This PR implements a two-stage lookup table optimization for BLASTN/megablast searches, following NCBI BLAST's approach. This provides significant performance improvements for large word sizes (e.g., `word_size=28` for megablast).

## Key Changes

### 1. Two-Stage Lookup Table Implementation

- **`lut_word_length = 8`**: Used for indexing (direct array access, 65,536 entries)
- **`word_length = 28`**: Used for extension triggering (for megablast)
- Implements `TwoStageLookup` struct that combines:
  - Direct array lookup using `PvDirectLookup` with 8-mer indexing
  - 28-mer match verification before triggering extensions

### 2. 28-mer Match Check

- Verifies full 28-mer match before triggering ungapped extension
- Dramatically reduces unnecessary extensions (8-mer matches that don't have full 28-mer matches)
- Critical optimization: ~65% performance improvement for EDL933 vs Sakai

### 3. Optimized Scan Step

- Implements NCBI BLAST formula: `scan_step = word_length - lut_word_length + 1 = 21` for megablast
- Reduces lookup operations by ~5x compared to previous `scan_step = 4`

### 4. Default Task Fix

- Removed redundant `--task megablast` from test scripts (default is already megablast)
- Updated output file names to match BLAST+ convention (`*.blastn.megablast.out`)

## Performance Results

| Test Case | Before | After | Improvement |
|-----------|--------|-------|-------------|
| NZ_CP006932 Self | 1.047s | **0.729s** | ~30% |
| EDL933 vs Sakai | 7.136s | **2.496s** | **~65%** |
| Sakai vs MG1655 | 9.848s | **5.312s** | **~46%** |

## Technical Details

### Architecture

1. **Lookup Table Construction**:
   - `build_two_stage_lookup()`: Creates lookup table with `lut_word_length=8`
   - Uses `PvDirectLookup` for O(1) direct array access
   - Presence Vector (PV) for fast filtering

2. **Scanning Loop**:
   - Uses 8-mer rolling k-mer for subject scanning
   - Performs 8-mer lookup (O(1) direct array access)
   - Verifies 28-mer match before extension
   - Applies two-hit filter and ungapped extension only for valid 28-mer matches

3. **Algorithm Selection**:
   - Automatically selects two-stage lookup for `word_size >= 16`
   - Falls back to direct lookup (word_size <= 13) or hash lookup as needed

## Reference

Based on NCBI BLAST implementation:
- `ncbi-blast/c++/src/algo/blast/core/blast_nascan.c`
- `ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h`

## Testing

All tests pass with correct output generation. Output files verified for:
- NZ_CP006932 Self
- EDL933 vs Sakai  
- Sakai vs MG1655

## Future Work

- Further optimization: 28-mer match check can be optimized with SIMD instructions
- Target: Achieve <1 second for all megablast test cases (currently 2.5-5.3 seconds)

