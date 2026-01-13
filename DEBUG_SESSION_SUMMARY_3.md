# BLASTN Debug Session Summary - Session 3

**Date**: 2026-01-09
**Branch**: fix2

## What Was Investigated

This session focused on fixing BLASTN to achieve NCBI parity:
- Target: 95%+ hit coverage (was 75%)
- Target: Fix coordinate off-by-one issue
- Target: Remove non-NCBI features

## Changes Made This Session

### 1. Removed Coverage Comparison Workaround (gapped.rs)

**File**: `LOSAT/src/algorithm/blastn/alignment/gapped.rs`

**What was removed**: A workaround in `extend_final_traceback` that compared preliminary vs final traceback coverage and kept whichever was better.

**Why removed**: This was NOT in NCBI code. User explicitly demanded: "DO NOT ASSUME. REFER TO THE NCBI CODE AND TRANSPILE FAITHFULLY"

**Impact**: After removal, query coordinates became off by 1:
- Before removal: `q=1-657101, s=1-657101` (both correct)
- After removal: `q=2-657100, s=1-657101` (query off by 1)

### 2. Added filter_contained_hsps (SHOULD BE REMOVED)

**File**: `LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs`

**What was added**: A containment filter to remove HSPs that are nearly contained within higher-scoring HSPs.

**Problem**: This is NOT in NCBI BLAST. It caused massive hit loss:
- 90% threshold: 9,385 → 3,789 hits (60% loss)
- 99% threshold: 9,385 → 3,904 hits (58% loss)

**Current state**: Function exists but is disabled in `blast_engine/mod.rs`

**Action needed**: DELETE this function entirely - it's non-NCBI

### 3. Previous Fixes (From Earlier in Session)

These fixes were correct and should be kept:

1. **Right-extension X-drop** (`extension.rs:76`): Removed `current_score_r <= 0` check. NCBI only uses X-drop termination.

2. **s_consumed calculation** (`gapped.rs:454, 707`): Fixed to `b_offset + 1` because b_index is 0-indexed.

## Current State

### Test Results
```
NCBI top hit:  1     657101  1       657101  (q_len=657101)
LOSAT top hit: 2     657100  1       657101  (q_len=657099)
```

### Hit Counts
- NCBI: 12,340 hits
- LOSAT: 9,360 hits (76% coverage)

### Coordinate Asymmetry Bug

For self-comparison (same sequence), query and subject should produce identical coordinates. But:
- Query: `2-657100` (off by 1 at both ends)
- Subject: `1-657101` (correct)

This asymmetry indicates a bug in `extend_gapped_one_direction_ex` where query and subject are handled differently.

## Root Cause Analysis (Incomplete)

The asymmetry was traced to `extend_gapped_one_direction_ex` in `gapped.rs`. Possible causes:

1. **get_q vs get_s helper functions**: May have different indexing logic
2. **score_val initialization**: May initialize query/subject indices differently
3. **a_index vs b_index handling**: a_index is 1-indexed, b_index is 0-indexed - this asymmetry may not be handled correctly

### Key NCBI Reference
`blast_gapalign.c:735-962` - `Blast_SemiGappedAlign`

The NCBI code uses symmetric handling for both sequences. LOSAT's implementation may have introduced asymmetry.

## Files That Need Attention

### Must Fix
1. `gapped.rs` - `extend_gapped_one_direction_ex`: Fix query/subject asymmetry

### Must Remove (Non-NCBI)
1. `purge_endpoints.rs` - `filter_contained_hsps`: Delete entire function
2. `blast_engine/mod.rs` - Remove disabled filter call and unused import

## What Was NOT Broken

These are working correctly:
1. Right-extension X-drop termination (matches NCBI)
2. s_consumed calculation (matches NCBI)
3. Subject coordinate output (correct)
4. Basic DP algorithm structure

## Next Session Tasks

### Priority 1: Fix Query Coordinate Bug
1. Read NCBI `blast_gapalign.c:735-962` line by line
2. Compare with `extend_gapped_one_direction_ex`
3. Find where query indexing diverges from subject
4. Look specifically at:
   - Initial `a_index` vs `b_index` setup
   - How `q_consumed` is calculated vs `s_consumed`
   - The `get_q` vs `get_s` helper closures

### Priority 2: Remove Non-NCBI Code
1. Delete `filter_contained_hsps` function
2. Remove its import/export from `filtering/mod.rs`
3. Remove disabled call from `blast_engine/mod.rs`

### Priority 3: Investigate Missing 24% Hits
1. Add diagnostic logging at each pipeline stage
2. Compare seed counts, ungapped extension counts, gapped extension counts
3. Identify where hits are being lost

## Commands for Next Session

```bash
# Build
cd LOSAT && cargo build --release

# Test
./target/release/LOSAT blastn -q tests/fasta/NZ_CP006932.fasta \
    -s tests/fasta/NZ_CP006932.fasta -o /tmp/test.out --task blastn -n 1

# Compare
echo "NCBI: $(grep -v '^#' tests/blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out | wc -l)"
echo "LOSAT: $(grep -v '^#' /tmp/test.out | wc -l)"
head -5 tests/blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out
head -5 /tmp/test.out
```

## Lesson Learned

**Never add features that don't exist in NCBI.** Even if they seem to "fix" an issue, they mask the underlying bug and create non-parity behavior. Always trace the root cause in the actual algorithm implementation.
