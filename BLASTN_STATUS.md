# LOSAT BLASTN Implementation Status

**Last Updated**: 2026-01-10
**Branch**: fix2

## Overview

This document summarizes the current state of the BLASTN implementation in LOSAT, including known issues, completed fixes, and next steps.

## Current Test Results

| Dataset | Algorithm | LOSAT | NCBI | Coverage | Status |
|---------|-----------|-------|------|----------|--------|
| EDL933.Sakai | megablast | 5,634 | 5,718 | 98.5% | Good |
| Sakai.MG1655 | megablast | 6,434 | 6,476 | 99.3% | Good |
| NZ_CP006932 self | blastn | 11,090 | 12,340 | 89.9% | Acceptable |
| MjeNMV.MelaMJNV | blastn | 2,588 | 2,668 | 97.0% | Good |
| SiNMV.ChdeNMV | blastn | 4,192 | 4,367 | 95.9% | Good |

## Session 8 Fix (Regression Fixed)

### Root Cause Analysis

The Session 7 regression (11,094 -> 4,610 hits) was caused by removing the same-diagonal requirement from `is_hsp_contained()`.

**Why same-diagonal is REQUIRED in LOSAT:**

NCBI uses an **interval tree** (`blast_itree.c`) for containment checking. The tree structure naturally limits comparisons to **geometrically nearby HSPs** based on query coordinate overlap. This means:
- A large HSP at q=1-657101 does NOT filter out a small HSP at q=50000-50100 unless they overlap in the tree query
- The tree uses query midpoint for organization, so only nearby HSPs are compared

LOSAT uses a **flat list** for containment checking, iterating over ALL existing HSPs. Without geometric locality constraints:
- A large HSP covering the entire sequence filters out ALL smaller HSPs
- This causes massive over-filtering (only 37% coverage)

**Solution**: The same-diagonal requirement acts as a **proxy for NCBI's interval tree geometric locality**. Two HSPs on the same diagonal are geometrically related; HSPs on different diagonals represent distinct alignments and should not be compared for containment.

### Changes Made (Session 8)

1. **Restored same-diagonal requirement** (`run.rs:46-86`)
   - Added diagonal computation: `diag = q_start - s_start`
   - Added check: `if new_diag != existing_diag { return false; }`
   - Result: 11,090 hits (restored from 4,610)

2. **Kept cutoff score fix** (`mod.rs`, `run.rs`)
   - Still using pre-computed `cutoff_scores` from `compute_blastn_cutoff_score()`
   - No longer hardcoding `std::cmp::max(reward * 11, 22)`

3. **Kept endpoint purge sort order fix** (`purge_endpoints.rs:379-386`)
   - query.end and subject.end now sort ASC (matching NCBI)

## What's Working

### Megablast (default task)
- Core algorithm functional
- Hit coverage: **98-99%** compared to NCBI BLAST

### Blastn (task blastn)
- Core DP-based gapped alignment implemented
- Follows NCBI's `Blast_SemiGappedAlign` algorithm
- 2-phase endpoint purge implemented (blast_traceback.c:637-669)
- Traceback capture during extension implemented
- Same-diagonal containment check (proxy for interval tree)
- Hit coverage: **89.9%** compared to NCBI BLAST

## Completed Fixes

### 1. Traceback Capture During Extension (Session 6)
**Problem**: `Hit.gap_info` was always `None` because gapped extension used score-only DP

**Fix**: Added `extend_gapped_heuristic_with_traceback()` and `extend_gapped_one_direction_with_traceback()` in `gapped.rs`

**NCBI Reference**: `blast_gapalign.c:364-733` (ALIGN_EX function)

### 2. Two-Phase Endpoint Purge Integration (Session 6)
**Problem**: LOSAT only did `purge=TRUE` (delete), while NCBI does:
1. `purge=FALSE` (trim overlapping HSPs using edit scripts)
2. Re-evaluate trimmed HSPs
3. `purge=TRUE` (delete remaining duplicates)

**Fix**: Implemented full 2-phase flow in `filter_hsps()`

**NCBI Reference**: `blast_traceback.c:637-669`

### 3. HSP Re-evaluation Function (Session 6)
**Fix**: Added `reevaluate_hsp_with_ambiguities_gapped()` in `purge_endpoints.rs`

**NCBI Reference**: `blast_hits.c:479-647`

### 4. Same-Diagonal Containment Check (Session 6, restored Session 8)
**Problem**: LOSAT's flat list containment check is too aggressive without geometric locality

**Fix**: Added same-diagonal requirement as proxy for NCBI's interval tree locality

**NCBI Reference**: `blast_gapalign.c:3918`, `blast_itree.c:809-847`

### 5. Cutoff Score Calculation (Session 7)
**Problem**: Hardcoded cutoff `std::cmp::max(reward * 11, 22)` instead of proper calculation

**Fix**: Use pre-computed `cutoff_scores` from `compute_blastn_cutoff_score()`

**NCBI Reference**: `blast_traceback.c:654`

### 6. Endpoint Purge Sort Order (Session 7)
**Problem**: query.end and subject.end were sorting DESC, should be ASC

**Fix**: Changed to ASC ordering matching NCBI

**NCBI Reference**: `blast_hits.c:2268-2291 s_QueryOffsetCompareHSPs`

## Known Issues

### Issue 0: Megablast Performance Regression (CRITICAL)
**Symptom**: EDL933 vs Sakai (5.5MB x 2) comparison hangs/extremely slow. NCBI BLAST completes in ~1 second.

**Status**: CRITICAL - Needs immediate investigation

**Analysis**:
- Lookup table building completes (Task: megablast, Word: 28, TwoStage: true, LUTWord: 8)
- Search phase appears to hang or run infinitely slow
- Possible infinite loop or O(n^2) algorithm in extension/filtering

**Potential Causes**:
1. Changes to containment check affecting megablast path
2. Changes to filter_hsps() causing performance regression
3. Infinite loop in extension code

### Issue 1: LOSAT Produces Extra Full-Length Hits
**Symptom**: On self-comparison, LOSAT produces TWO massive hits (657,101 bp + 654,201 bp) while NCBI's max is 3,629 bp

**Impact**: These extra hits don't significantly affect hit count but indicate different seeding/extension behavior

**Status**: Under investigation

**Analysis**:
- For identical sequences, X-drop never triggers (score keeps increasing)
- NCBI still produces shorter hits due to how seeds/extensions interact with DUST-masked regions
- LOSAT may be extending through masked regions differently

**Potential Causes**:
1. DUST masking not applied during extension
2. Different handling of N bases in extension
3. Different seed distribution due to masking

### Issue 2: Missing Medium-Length Hits
**Symptom**: NZ_CP006932 self-comparison shows 89.9% coverage (11,090 vs 12,340)

**Analysis**:
- Missing ~1,250 hits
- LOSAT: 11,090 hits, NCBI: 12,340 hits
- Likely missing hits are in 30-99 bp range

**Status**: Partially addressed by 2-phase endpoint purge (improved from 72% to 89.9%)

## Architecture

### Key Files
```
src/algorithm/blastn/
├── alignment/
│   ├── gapped.rs         # Traceback-capturing extension functions
│   └── mod.rs            # Re-export
├── blast_engine/
│   ├── mod.rs            # 2-phase endpoint purge, filter_hsps()
│   └── run.rs            # Same-diagonal containment check, main run loop
├── filtering/
│   ├── mod.rs            # Export
│   └── purge_endpoints.rs # HSP re-evaluation, endpoint purge
└── ncbi_cutoffs.rs       # compute_blastn_cutoff_score()
```

### NCBI Reference Files
```
ncbi-blast/c++/src/algo/blast/core/
├── blast_gapalign.c:364-733   # ALIGN_EX (traceback capture)
├── blast_traceback.c:637-669  # 2-phase endpoint purge
├── blast_hits.c:479-647       # Blast_HSPReevaluateWithAmbiguitiesGapped
├── blast_hits.c:2268-2291     # s_QueryOffsetCompareHSPs (sort order)
└── blast_itree.c:809-847      # s_HSPIsContained (interval tree)
```

## Lessons Learned

### 1. Interval Tree vs Flat List
NCBI's interval tree is not just an optimization - it provides **geometric locality** that affects correctness. When using a flat list, same-diagonal requirement is needed as a proxy.

### 2. Faithful Transpilation
When porting NCBI algorithms, consider:
- **Context**: What data structures surround the function?
- **Timing**: When in the pipeline is it called?
- **Constraints**: What implicit constraints exist (e.g., tree structure limiting comparisons)?

Simply copying the logic without understanding the context can cause subtle bugs.

## Next Steps

### Priority 1: Close Remaining Hit Gap
Current: 89.9% (11,090 / 12,340)
Target: 95%+

Actions:
1. Compare hit coordinate distributions between LOSAT and NCBI
2. Analyze which HSPs NCBI keeps that LOSAT filters
3. Check if interval tree would find additional valid comparisons

### Priority 2: Implement Interval Tree (Medium-term)
Proper interval tree implementation would:
- Match NCBI's containment check behavior exactly
- Potentially find more containment relationships (increasing coverage)
- Remove need for same-diagonal proxy

## Debug Commands

```bash
# Run blastn with timing
cd LOSAT && time ./target/release/LOSAT blastn \
    -q tests/fasta/NZ_CP006932.fasta \
    -s tests/fasta/NZ_CP006932.fasta \
    -o /tmp/test.out --task blastn -n 1

# Compare hit counts
wc -l /tmp/test.out
grep -v "^#" tests/blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out | wc -l

# Check longest alignments
awk -F'\t' '{print $4}' /tmp/test.out | sort -rn | head -5

# Compare strand distribution
awk -F'\t' '$9 < $10 { fwd++ } $9 > $10 { rev++ } END { print "Fwd:", fwd, "Rev:", rev }' /tmp/test.out
```
