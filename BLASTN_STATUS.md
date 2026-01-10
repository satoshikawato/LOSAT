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
| NZ_CP006932 self | blastn | **4,610** | 12,340 | **37.4%** | **REGRESSION** |
| MjeNMV.MelaMJNV | blastn | 2,588 | 2,668 | 97.0% | Good |
| SiNMV.ChdeNMV | blastn | 4,192 | 4,367 | 95.9% | Good |

## Session 7 Changes (CAUSED REGRESSION)

### Changes Made

1. **Cutoff score hardcode removed** (`mod.rs:127-134`)
   - Removed `std::cmp::max(reward * 11, 22)` hardcode
   - Now uses pre-computed cutoff_scores from `compute_blastn_cutoff_score()`
   - NCBI reference: `blast_traceback.c:654`

2. **Endpoint purge sort order fixed** (`purge_endpoints.rs:379-386`)
   - Fixed: query.end and subject.end were DESC, should be ASC
   - NCBI reference: `blast_hits.c:2268-2291 s_QueryOffsetCompareHSPs`

3. **Containment check same-diagonal requirement removed** (`run.rs:36-88`)
   - Removed: `if new_diag != existing_diag { return false; }`
   - NCBI reference: `blast_itree.c:809-847 s_HSPIsContained`
   - **This change likely caused the regression**

### Regression Analysis

Before session 7: 11,094 hits (89.9%)
After session 7: **4,610 hits (37.4%)**

The containment check now filters too aggressively. While NCBI does not require same-diagonal in `s_HSPIsContained`, there may be other constraints (e.g., query context, interval tree structure) that prevent over-filtering.

### Recommended Rollback

The same-diagonal requirement in `is_hsp_contained()` should be restored until further investigation. The NCBI interval tree implementation has additional constraints that prevent the containment check from being too aggressive.

## What's Working

### Megablast (default task)
- Core algorithm functional
- Hit coverage: **98-99%** compared to NCBI BLAST

### Blastn (task blastn)
- Core DP-based gapped alignment implemented
- Follows NCBI's `Blast_SemiGappedAlign` algorithm
- 2-phase endpoint purge implemented (blast_traceback.c:637-669)
- Traceback capture during extension implemented

## Completed Fixes (Session 6)

### 1. Traceback Capture During Extension
**Problem**: `Hit.gap_info` was always `None` because gapped extension used score-only DP

**Fix**: Added `extend_gapped_heuristic_with_traceback()` and `extend_gapped_one_direction_with_traceback()` in `gapped.rs` that capture edit scripts during alignment

**NCBI Reference**: `blast_gapalign.c:364-733` (ALIGN_EX function)

**Files Modified**:
- `src/algorithm/blastn/alignment/gapped.rs`
- `src/algorithm/blastn/alignment/mod.rs`

### 2. Two-Phase Endpoint Purge Integration
**Problem**: LOSAT only did `purge=TRUE` (delete), while NCBI does:
1. `purge=FALSE` (trim overlapping HSPs using edit scripts)
2. Re-evaluate trimmed HSPs
3. `purge=TRUE` (delete remaining duplicates)

**Fix**: Implemented full 2-phase flow in `filter_hsps()`:
```rust
// Phase 1: purge=FALSE - trim overlapping HSPs
let (mut result_hits, extra_start) = purge_hsps_with_common_endpoints_ex(result_hits, false);

// Phase 2: Re-evaluate trimmed HSPs
for i in extra_start..result_hits.len() {
    let delete = reevaluate_hsp_with_ambiguities_gapped(...);
    if delete { hit.raw_score = i32::MIN; }
}
result_hits.retain(|h| h.raw_score != i32::MIN);

// Phase 3: purge=TRUE for BLASTN
let (result_hits, _) = purge_hsps_with_common_endpoints_ex(result_hits, true);
```

**NCBI Reference**: `blast_traceback.c:637-669`

**Files Modified**:
- `src/algorithm/blastn/blast_engine/mod.rs`
- `src/algorithm/blastn/filtering/purge_endpoints.rs`
- `src/algorithm/blastn/filtering/mod.rs`

### 3. HSP Re-evaluation Function
**Problem**: Missing `Blast_HSPReevaluateWithAmbiguitiesGapped()` equivalent

**Fix**: Added `reevaluate_hsp_with_ambiguities_gapped()` in `purge_endpoints.rs`

**NCBI Reference**: `blast_hits.c:479-647`

### 4. Same-Diagonal Containment Check
**Problem**: Duplicate hits on same diagonal (e.g., q=1-657101 and q=2-657100 on diagonal 0)

**Fix**: Added `is_hsp_contained_same_diagonal()` check DURING extension to filter hits that are contained within existing HSPs on the SAME diagonal only. This matches NCBI's interval tree behavior where geometrically close HSPs are compared.

**NCBI Reference**: `blast_gapalign.c:3918`, `blast_itree.c:809-847`

**Files Modified**:
- `src/algorithm/blastn/blast_engine/run.rs`

## Known Issues

### Issue 1: LOSAT Produces Full-Length Diagonal Hit
**Symptom**: On self-comparison, LOSAT produces a 657,101 bp full-length hit while NCBI's max is 3,629 bp

**Impact**: The massive full-length hit affects hit count comparison

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
**Symptom**: NZ_CP006932 self-comparison shows 89.9% coverage (11,094 vs 12,340)

**Analysis**:
- Missing hits are in 30-99 bp range
- LOSAT: 11,094 hits, NCBI: 12,340 hits
- Gap of ~1,246 hits

**Status**: Partially addressed by 2-phase endpoint purge (improved from 72% to 89.9%)

## Architecture

### Files Modified in Session 6
```
src/algorithm/blastn/
├── alignment/
│   ├── gapped.rs         # Added traceback-capturing extension functions
│   └── mod.rs            # Re-export new functions
├── blast_engine/
│   ├── mod.rs            # 2-phase endpoint purge flow
│   └── run.rs            # Same-diagonal containment check, traceback usage
├── filtering/
│   ├── mod.rs            # Export new functions
│   └── purge_endpoints.rs # reevaluate_hsp_with_ambiguities_gapped()
└── common.rs             # GapEditOp enum (already existed)
```

### NCBI Reference Files
```
ncbi-blast/c++/src/algo/blast/core/
├── blast_gapalign.c:364-733   # ALIGN_EX (traceback capture)
├── blast_traceback.c:637-669  # 2-phase endpoint purge
├── blast_hits.c:479-647       # Blast_HSPReevaluateWithAmbiguitiesGapped
└── blast_itree.c:809-847      # s_HSPIsContained
```

## Next Steps

### Priority 1: Investigate Full-Length Hit Issue
Why does LOSAT extend to full sequence length while NCBI stops at ~3,600 bp?

Actions:
1. Check if DUST-masked sequence is used for extension
2. Compare seed positions and extension boundaries
3. Check X-drop handling with masked regions

### Priority 2: Close Remaining Hit Gap
Current: 89.9% (11,094 / 12,340)
Target: 98%+

Actions:
1. Compare hit coordinate distributions
2. Analyze which HSPs NCBI keeps that LOSAT filters
3. Fine-tune endpoint purge logic

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
