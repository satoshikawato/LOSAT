# TBLASTX 2x Hit Count Investigation

## Problem Statement

LOSAT produces approximately twice as many alignment hits as NCBI BLAST for long sequences (~600kb).

| Metric | LOSAT | NCBI | Ratio |
|--------|-------|------|-------|
| Hit count (long seq) | 29,756 | 14,877 | **2.00x** |
| Hit count (short seq) | 23,766 | 23,715 | **1.002x** |
| Long test case | AP027131 vs AP027133 (gencode 4, 671kb vs 615kb) | | |
| Short test case | MjeNMV vs MelaMJNV (gencode 1, ~125kb) | | |

**Key observation**: Short sequences have near-parity (~1.002x). The issue scales with sequence length.

---

## BUGS FOUND AND FIXED (2026-01-08 Session 5)

### 1. Neighbor Word Expansion - Now Matches NCBI Exactly

**File:** `lookup.rs:431-445, 489-507`

**Issue:** Previous "fix" was WRONG. NCBI blast_aalookup.c s_AddWordHits (lines 504-519) DOES add explicit entries when `self_score < threshold`. This is correct because neighbor search won't find exact matches for low-scoring words.

**NCBI logic:**
1. If threshold==0 OR self_score < threshold: Add explicit entry for query word
2. If threshold==0: Return (skip neighbor search)
3. Otherwise: Do neighbor search for all s where score(w,s) >= threshold

**Current LOSAT (now matches NCBI):**
```rust
// Compute self-score of query word
let self_score = blosum62_score(w0, w0) + blosum62_score(w1, w1) + blosum62_score(w2, w2);

// NCBI: If threshold==0 OR self_score < threshold, add exact match explicitly
if threshold == 0 || self_score < threshold {
    entry_counts[idx] += num_offsets;
}
// NCBI: If threshold==0, skip neighbor search entirely
if threshold == 0 {
    continue;
}
// Neighbor search
```

### 2. Flag Update Condition - Now Matches NCBI Exactly

**File:** `utils.rs:2001-2019`

**Bug:** LOSAT was setting `flag=1` based on `right_extend` (whether left extension reached first hit), but NCBI sets it based on `score >= cutoff_score`.

**Impact:** Affects diagonal suppression behavior for subsequent hits on the same diagonal.

**NCBI reference (aa_ungapped.c:575-591):**
```c
if (score >= cutoffs->cutoff_score) {
    diag_array[diag_coord].flag = 1;
    diag_array[diag_coord].last_hit = s_last_off - (wordsize - 1) + diag_offset;
} else {
    diag_array[diag_coord].last_hit = subject_offset + diag_offset;
}
```

**Current LOSAT (now matches NCBI):**
```rust
if score >= cutoff {
    diag_entry.flag = 1;
    diag_entry.last_hit = s_last_off - (wordsize - 1) + diag_offset;
} else {
    diag_entry.last_hit = subject_offset + diag_offset;
}
```

---

## REMAINING MYSTERY (As of 2026-01-08 Session 5)

**Short sequence test (MjeNMV vs MelaMJNV):** LOSAT 23,766 vs NCBI 23,715 = **1.002x** (near parity)

**Long sequence test (AP027131 vs AP027133):** LOSAT 29,756 vs NCBI 14,877 = **2.00x** (double)

The issue is clearly **LENGTH-DEPENDENT**. With fixes applied, short sequences have near-parity.

### Key Observations
1. Short sequences: Near-parity (~1.002x)
2. Long sequences: 2x excess
3. Excess concentrated in low-score HSPs (0-30 bits)
4. Both neighbor expansion and flag update now match NCBI exactly
5. SIMD vs scalar produces identical results
6. No duplicate offset pairs in scan output

### What's Still Different?
The remaining difference must scale with sequence length. Possibilities:
1. Diagonal array state management over many frames
2. Integer precision issues in long sequences
3. Some NCBI-specific optimization not yet identified
4. Different handling of frame boundaries or sentinel bytes

---

## VERIFIED CORRECT (Eliminated as Root Cause)

| Component | Verification | NCBI Reference |
|-----------|--------------|----------------|
| Threshold | 13 in both | `blast_aalookup.c` |
| BLOSUM62 matrix | Verbatim copy | `src/utils/matrix.rs` |
| Frame generation | 6 frames each | `translation.rs` |
| Lookup table entries | No duplicates | `lookup.rs:333-367` |
| Output writing | NOT duplicated | User confirmed |
| Two-hit logic (flag branch) | Matches NCBI | `utils.rs:1844-1858` vs `aa_ungapped.c:519-530` |
| Two-hit logic (else branch) | Matches NCBI | `utils.rs:1861-1901` vs `aa_ungapped.c:533-573` |
| Diagonal formula | `(q - s) & mask` | `aa_ungapped.c:516` |
| Extension algorithm | Score matches | `extension.rs` |
| Cutoff calculation | Three-stage cap | `ncbi_cutoffs.rs` |
| Karlin params | lambda=0.3176, K=0.134 | |
| X-drop | 16 (7.0 bits) | |
| HSP culling | disabled (0) | |
| Reevaluation trimming | Matches NCBI | `reevaluate.rs` |
| Scan function | Matches NCBI | `utils.rs:840-1160` |
| get_context_idx | Binary search | `lookup.rs:178-191` |
| Offset pair generation | Matches NCBI | `utils.rs:908-943` |
| Per-frame diag_offset increment | Matches NCBI verbatim | `utils.rs:2088-2097` |
| diag_array sizing/collision | Tested 2x size, no change | `utils.rs:1603-1611` |
| Neighbor expansion (exact match bug) | **FIXED** | `lookup.rs:431-489` |

---

## NOT YET VERIFIED (Remaining Possibilities)

| Component | Status | File |
|-----------|--------|------|
| Neighbor word score calculation | **VERIFIED** (matches NCBI) | `lookup.rs:505-524` |
| Two-hit window update logic | **VERIFIED** (matches NCBI) | `utils.rs:1982-1998` |
| Strand distribution | **VERIFIED** (ratios consistent 1.75x-2.17x across all strands) | N/A |
| Scan range initialization | **VERIFIED** (1 to len-4) | `utils.rs:1798` |

### Remaining Mystery

After fixing the self_score < threshold bug, the ratio dropped from 2.00x to 1.94x. However, the remaining ~94% excess is still unaccounted for. The issue:

1. Affects all strand combinations proportionally (++, +-, -+, --)
2. Scales with sequence length (short sequences now have 0.97x, large have 1.94x)
3. Is NOT caused by:
   - Hash table collisions (tested with 2x table size)
   - Diagonal overflow handling
   - Two-hit window logic
   - Cutoff calculations
   - E-value filtering

---

## Current Diagnostic Output (After Fix)

```
Exact entries added (self_score < threshold): 0  (was 126,185)
Neighbor entries added: 10,969,813
K-mer matches found: 1,124,805,734  (was 1,200,397,107)
Total extensions: 28,465,892
Total linked HSPs: 329,196
Final hits after E-value filter: 28,819  (NCBI: 14,877)
```

---

## Debug Output Analysis

```
Total linked HSPs: 329,196  (was 338,859)
Final hits after E-value filter: 28,819  (was 29,766)
Filtering ratio: ~11x
```

- NCBI has 14,877 final hits
- If same filtering ratio: NCBI starts with ~163k linked HSPs
- LOSAT has ~2x the initial HSPs BEFORE filtering
- **Conclusion: Problem is in HSP GENERATION, not filtering**

---

## Strand Distribution (After Fix - 1.94x)

| Strand | LOSAT | NCBI | Ratio |
|--------|-------|------|-------|
| ++ | 7,881 | 4,039 | 1.95x |
| -+ | 7,272 | 3,766 | 1.93x |
| -- | 6,875 | 3,932 | 1.75x |
| +- | 6,791 | 3,134 | 2.17x |

---

## Key Files

| Component | LOSAT File | NCBI Reference |
|-----------|------------|----------------|
| Scan function | `utils.rs:840-1160` | `blast_aalookup.c:BlastAaScanSubject` |
| Two-hit detection | `utils.rs:1817-1974` | `aa_ungapped.c:502-610` |
| Lookup building | `lookup.rs:235-600` | `blast_aalookup.c` |
| Extension | `extension.rs` | `aa_ungapped.c` |

---

## Composition-Based Statistics

LOSAT implements amino acid composition calculation for Karlin parameters (`lookup.rs:252-275`), but does NOT implement:
- Composition-based scoring adjustment during alignment
- NCBI's `blast_kappa.c` post-alignment adjustments

**Status:** Composition-based adjustments affect E-values, not initial hit generation. Unlikely to be the 2x cause.

## Diagnostic Log Saved

Full diagnostics: `tests/losat_out/AP027131.AP027133.tlosatx.diag.log`

Key metrics:
- Exact positions indexed: 862,165
- Unique words: 9,255
- Neighbor entries: 10,969,813
- Expansion factor: 12.72x

## CRITICAL FINDING: Score Distribution Analysis (2026-01-08 Session 4)

**THE 2x EXCESS IS CONCENTRATED IN LOW-SCORE HITS**

| Bit Score Range | NCBI | LOSAT | Ratio | Difference |
|-----------------|------|-------|-------|------------|
| **0-30** | 8,476 | 20,862 | **2.46x** | +12,386 |
| 30-35 | 1,743 | 2,853 | 1.63x | +1,110 |
| 35-40 | 1,107 | 1,210 | 1.09x | +103 |
| 40-45 | 734 | 791 | 1.07x | +57 |
| 45-50 | 474 | 487 | 1.02x | +13 |
| 50-60 | 614 | 651 | 1.06x | +37 |
| 60-80 | 644 | 763 | 1.18x | +119 |
| 80-100 | 291 | 333 | 1.14x | +42 |
| 100+ | 788 | 869 | 1.10x | +81 |
| **TOTAL** | 14,871 | 28,819 | **1.94x** | +13,948 |

**Key insight:** 88% of the excess hits (12,386/13,948) are in the 0-30 bit score range.

### Coordinate Comparison

| Category | Count |
|----------|-------|
| Exact matches (same coords) | 10,885 (73% of NCBI) |
| NCBI only (missing in LOSAT) | 3,986 (27% of NCBI) |
| LOSAT only (extra) | 17,934 (62% of LOSAT) |

**Observation:** LOSAT has both missing hits AND extra hits. The extra hits dominate.

---

## VERIFIED: NOT Duplicate Generation

- SIMD disabled (LOSAT_NO_SIMD=1): Same 28,819 hits
- Offset pair uniqueness check: 0 duplicates in scan output
- Scan function verified: No duplicate (q_off, s_off) pairs returned

**Conclusion:** The issue is NOT in k-mer scanning or SIMD boundary handling.

---

## CHANGES MADE (2026-01-08 Session 6)

### 3. Per-Subject Diagonal Array Reset

**File:** `utils.rs:1637-1655`

**Change:** Reset diagonal array for EACH subject sequence to match NCBI's Blast_ExtendWordNew behavior.

**NCBI behavior:** Creates a new ewp (diagonal table) for each subject via calloc (zeros all entries).

**LOSAT now matches:**
```rust
// Reset diagonal state for each subject to match NCBI behavior:
for d in st.diag_array.iter_mut() {
    *d = DiagStruct::default();
}
st.diag_offset = window;
```

**Result:** No change in hit count. Still 29,756 vs 14,877 = 2.0x

### 4. Per-Frame init_hsps Reset

**File:** `utils.rs:1787-1789, 2135-2150`

**Change:** Reset init_hsps between subject frames to match NCBI's BlastInitHitListReset behavior.

**NCBI behavior:**
- blast_engine.c:491 calls `BlastInitHitListReset(init_hitlist)` before each frame
- blast_engine.c:561-584 calls `BLAST_GetUngappedHSPList` per-frame, then `Blast_HSPListsMerge`

**LOSAT now matches:**
```rust
for (s_f_idx, s_frame) in s_frames.iter().enumerate() {
    // NCBI: BlastInitHitListReset(init_hitlist) - reset per-frame
    let mut init_hsps: Vec<InitHSP> = Vec::new();

    // ... two-hit detection and extension ...

    // NCBI: BLAST_GetUngappedHSPList - per-frame conversion
    if !init_hsps.is_empty() {
        let frame_ungapped_hits = get_ungapped_hsp_list(init_hsps, ...);
        // NCBI: Blast_HSPListsMerge - merge per-frame results
        combined_ungapped_hits.extend(frame_ungapped_hits);
    }
}
```

**Result:** No change in hit count. Still 29,756 vs 14,877 = 2.0x

---

## ROOT CAUSE HYPOTHESIS (Updated)

The excess is in LOW-SCORE HSPs (0-30 bits), not high-score ones. This suggests:

1. **Different HSP generation at the boundary** - LOSAT may be generating HSPs that NCBI discards early
2. **Different extension termination** - LOSAT's x-drop may differ slightly
3. **Different E-value/filtering for low-score HSPs** - Though E-values look similar
4. **Different sum-statistics linking** - Chain formation may differ for low-score HSPs

### Verified NOT the cause:
- Per-subject diagonal array reset (now matches NCBI)
- Per-frame init_hsps reset (now matches NCBI)
- Neighbor word expansion (now matches NCBI)
- Flag update condition (now matches NCBI)

The 2x issue persists despite matching all verified NCBI behaviors. The root cause is still unknown.

## Next Investigation Steps

1. **Compare extension scores** - Same (q_start, s_start) should produce same raw score
2. **Check for chain-related duplicates** - Low-score HSPs may be fragments of chains
3. **Compare X-drop termination** - Different x-drop could produce different HSP lengths
4. **Focus on 20-30 bit score range** - This is where the biggest ratio difference is
5. **Compare two-hit trigger rates** - Verify same seeds trigger two-hit detection

---

*Last updated: 2026-01-08 (Session 6)*
