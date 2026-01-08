# Plan: Fix LOSAT vs NCBI BLAST 2x Hit Count Discrepancy

## Problem Summary

LOSAT produces ~2x more hits than NCBI BLAST for long sequences (600kb+), while short sequences (125kb) have near-parity (1.002x). 88% of the excess is in low-score HSPs (0-30 bits).

## CRITICAL BUGS FOUND (Thorough Code Comparison)

### BUG 1: Two-Hit Extension Returns Wrong Score (HIGH PRIORITY)
**File:** `extension.rs:596`

**NCBI (aa_ungapped.c:1157):**
```c
return MAX(left_score, right_score);  // Returns BEST of both extensions
```

**LOSAT (line 596):**
```rust
return (..., max_score_total, ...);  // Only returns right extension score
```

**Impact:** When right extension degrades score, LOSAT returns lower value than NCBI, potentially accepting/rejecting different HSPs.

### BUG 2: Missing Final Output Sorts in Sum-Stats Linking (HIGH PRIORITY)
**File:** `sum_stats_linking.rs`

**NCBI (link_hsps.c:990-1000):** Two qsort passes after all linking:
```c
qsort(..., s_RevCompareHSPsTransl);  // Reverse sort by context
qsort(..., s_FwdCompareHSPsTransl);  // Forward sort
```

**LOSAT:** Only does ONE initial sort. Missing final two sorts entirely.

**Impact:** Output order differs, may affect downstream hit selection and filtering.

### BUG 3: Extra Boundary Checks in Word Scanning (MEDIUM PRIORITY)
**File:** `extension.rs:428-430, 520-521`

**LOSAT (not in NCBI):**
```rust
if q_pos + i >= q_seq.len() || s_pos + i >= s_seq.len() {
    break;  // Early termination
}
```

**NCBI:** No such checks - assumes caller provides valid coordinates.

**Impact:** LOSAT may terminate word scanning early, missing valid scoring positions.

### BUG 4: Hardcoded Word Size (LOW PRIORITY)
**File:** `extension.rs:512`

**LOSAT:** `let k_size = 3;` (hardcoded)
**NCBI:** Uses `word_size` parameter

**Impact:** Less flexible, but 3 is correct for TBLASTX.

### BUG 5: Missing Positives Calculation (LOW PRIORITY for hit count)
**File:** `reevaluate.rs:480-519`

**NCBI (blast_hits.c:786-787):**
```c
if (matrix[*q][*s] > 0)
    num_pos++;  // Counts BLOSUM62 positives
```

**LOSAT:** Returns `(num_ident, align_length)` only - no positives count.

**Impact:** Output format may have incorrect positive counts, but doesn't affect hit generation.

---

## Root Cause Analysis (Original)

### Discrepancy 1: Missing Hit List Capacity Management (MEDIUM PRIORITY)

**NCBI implements heap-based capacity limits at two levels:**

1. **Per-Subject HSP Limit (`hsp_num_max`)**
   - `blast_hits.c:Blast_HSPListSaveHSP` uses a MIN-HEAP when HSP count exceeds capacity
   - Low-scoring HSPs are evicted when capacity is reached
   - For tblastx with `gapped_calculation=TRUE`: `hsp_num_max = INT4_MAX` (unlimited)

2. **Hit List Capacity (`prelim_hitlist_size`)**
   - `blast_hits.c:Blast_HitListUpdate` limits subjects with hits
   - For `hitlist_size=500` with gapped: `prelim_hitlist_size = MIN(MAX(2*500, 10), 500+50) = 550`
   - When full, uses E-value heap to evict worst subjects

**LOSAT Status:** No capacity limits - accumulates ALL hits that pass E-value threshold.

**Important Note:** For single-subject comparisons (one query vs one subject), the hit list capacity at subject level (550) doesn't apply since there's only 1 subject. The 2x issue in the test case (AP027131 vs AP027133) is a single-subject comparison, so capacity limits may NOT be the primary cause.

**NCBI Defaults for tblastx:**
- `hitlist_size`: 500 (user-configurable via `-max_target_seqs`)
- `prelim_hitlist_size`: 550 (for gapped, 500 for ungapped) - limits SUBJECTS, not HSPs
- `hsp_num_max`: INT4_MAX (unlimited per-subject HSPs for gapped searches)

For the 2x issue with single-subject comparisons, we need to look elsewhere.

### Discrepancy 1b: Sum-Statistics Linking Differences (HIGH PRIORITY - Single Subject)

The investigation document notes:
> **Known Issues:** Chain formation differences - E-value mismatches with short HSPs

For single-subject comparisons, the difference must be in:
1. **Chain formation** - How HSPs are grouped into chains
2. **Chain member handling** - Which chain members are kept vs filtered
3. **E-value calculation for chains** - Different sum-stats could filter differently

LOSAT may be:
- Creating more/different chains than NCBI
- Keeping chain members that NCBI filters
- Calculating different E-values for chains

The investigation noted:
> LOSAT has both **missing hits** AND **extra hits**. The extra hits dominate.

This suggests fundamental differences in HSP generation, not just filtering.

### Discrepancy 2: Missing Preliminary E-value Filtering

**NCBI:** `blast_engine.c:s_Blast_HSPListReapByPrelimEvalue()` filters HSPs using `hit_params->prelim_evalue` per-subject BEFORE adding to results.

**LOSAT:** Only filters by E-value at final output stage.

### Discrepancy 3: Potential Scan Resumption Differences

**NCBI:** When offset array fills (4096 entries), scan suspends and resumes from that position. The resumption logic in `blast_aalookup.c` may differ subtly.

**LOSAT:** Implements similar pattern but may have subtle differences in boundary handling.

---

## Implementation Plan

### FIX 1: Two-Hit Extension Score Bug (HIGH PRIORITY)

**File:** [extension.rs](LOSAT/src/algorithm/tblastx/extension.rs)

**Current Code (line ~596):**
```rust
if reached_first_hit {
    let (new_max_score, new_right_disp, j) = extend_right_ungapped(...);
    max_score_total = new_max_score;
}
return (..., max_score_total, ...);
```

**Fix:**
```rust
let mut left_score = max_score;  // Score from left extension
let mut right_score = 0i32;

if reached_first_hit {
    let (new_max_score, new_right_disp, j) = extend_right_ungapped(...);
    right_score = new_max_score;
}

// NCBI: return MAX(left_score, right_score)
let final_score = left_score.max(right_score);
return (..., final_score, ...);
```

### FIX 2: Add Final Output Sorts in Sum-Stats Linking (HIGH PRIORITY)

**File:** [sum_stats_linking.rs](LOSAT/src/algorithm/tblastx/sum_stats_linking.rs)

**NCBI Reference (link_hsps.c:990-1000):**
```c
if (kTranslatedQuery) {
    qsort(link_hsp_array, num, sizeof(BlastHSP*), s_RevCompareHSPsTransl);
    qsort(link_hsp_array, num, sizeof(BlastHSP*), s_FwdCompareHSPsTransl);
}
```

**Add at end of `link_hsp_chains` function:**
```rust
// NCBI: Two final sorts for translated queries (link_hsps.c:990-1000)
// Sort 1: Reverse by context/frame
result.sort_by(|a, b| {
    let a_ctx = a.ctx_idx / 3;  // Query context group
    let b_ctx = b.ctx_idx / 3;
    match b_ctx.cmp(&a_ctx) {  // Reverse order
        std::cmp::Ordering::Equal => {
            // Secondary: by subject frame (descending)
            b.s_frame.cmp(&a.s_frame)
        }
        other => other,
    }
});

// Sort 2: Forward stable sort within groups
result.sort_by(|a, b| {
    let a_ctx = a.ctx_idx / 3;
    let b_ctx = b.ctx_idx / 3;
    a_ctx.cmp(&b_ctx)  // Forward order
});
```

### FIX 3: Remove Extra Boundary Checks (MEDIUM PRIORITY)

**File:** [extension.rs](LOSAT/src/algorithm/tblastx/extension.rs)

**Remove these checks (lines ~428-430, ~520-521):**
```rust
// DELETE THIS:
if q_pos + i >= q_seq.len() || s_pos + i >= s_seq.len() {
    break;
}
```

**Rationale:** NCBI assumes caller provides valid coordinates. The lookup table and offset generation should guarantee valid positions. These checks cause early termination.

### FIX 4: Add Positives Calculation (LOW PRIORITY)

**File:** [reevaluate.rs](LOSAT/src/algorithm/tblastx/reevaluate.rs)

**Update `get_num_identities_and_positives_ungapped`:**
```rust
pub fn get_num_identities_and_positives_ungapped(
    query_nomask: &[u8],
    subject: &[u8],
    q_off: usize,
    s_off: usize,
    hsp_len: usize,
) -> Option<(usize, usize, usize)> {  // (num_ident, num_pos, align_length)
    let mut num_ident = 0usize;
    let mut num_pos = 0usize;

    for i in 0..hsp_len {
        let q = query_nomask[q_off + i];
        let s = subject[s_off + i];
        if q == s {
            num_ident += 1;
        } else if blosum62_score(q, s) > 0 {
            num_pos += 1;
        }
    }

    // NCBI: num_pos_ptr = num_pos + num_ident
    Some((num_ident, num_pos + num_ident, hsp_len))
}
```

### FUTURE: Hit List Capacity Management (if still needed after FIX 1-3)

**Only needed for multi-subject database searches.** Single-subject comparisons won't be affected.

**Files:** `args.rs`, `utils.rs`

**Add parameters:**
- `hitlist_size`: 500 (NCBI: max subjects with hits)
- `max_hsps_per_subject`: 0 (unlimited for gapped)

**Implement:** Heap-based hit list culling per NCBI `blast_hits.c:3243-3299`

---

## Implementation Order

1. **FIX 1** - Two-hit extension score bug (most likely cause of 2x)
2. **FIX 2** - Final output sorts in sum-stats linking
3. **FIX 3** - Remove extra boundary checks in extension
4. **Test** - Verify ratio improvement after each fix
5. **FIX 4** - Add positives calculation (if output format needs it)
6. **Phase 2** - Hit list capacity management (if still needed for database searches)

---

## Files to Modify

| File | Changes | Priority |
|------|---------|----------|
| `LOSAT/src/algorithm/tblastx/extension.rs` | FIX 1: MAX(left,right) score; FIX 3: Remove boundary checks | HIGH |
| `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs` | FIX 2: Add final output sorts | HIGH |
| `LOSAT/src/algorithm/tblastx/reevaluate.rs` | FIX 4: Add positives calculation | LOW |
| `LOSAT/src/algorithm/tblastx/args.rs` | Add `hitlist_size`, `max_hsps_per_subject` | MEDIUM |
| `LOSAT/src/algorithm/tblastx/utils.rs` | Integrate capacity-limited accumulation | MEDIUM |

---

## Success Criteria

### After FIX 1 (Two-Hit Score Bug)
- [ ] `extend_hit_two_hit` returns `MAX(left_score, right_score)`
- [ ] Run test: `cargo test tblastx::`
- [ ] Compare hit count: expect improvement toward 1.x ratio

### After FIX 2 (Final Sorts)
- [ ] `link_hsp_chains` applies two final sorts
- [ ] Output order matches NCBI for test cases

### After FIX 3 (Boundary Checks)
- [ ] No early termination in word scanning
- [ ] Score calculation matches NCBI for edge cases

### Overall Verification
- [ ] Long sequence ratio drops from 2.0x to <1.2x
- [ ] Short sequence ratio remains <1.01x
- [ ] No regression in existing unit tests (`cargo test`)
- [ ] Score distribution analysis shows reduced low-score excess
- [ ] Coordinate comparison shows fewer "LOSAT only" hits

### Test Commands
```bash
# Build
cd LOSAT && cargo build --release

# Unit tests
cargo test tblastx::

# Long sequence test
./target/release/LOSAT tblastx \
    -q tests/fasta/AP027131.fasta \
    -s tests/fasta/AP027133.fasta \
    --query-gencode 4 --db-gencode 4 \
    -o tests/losat_out/long_test.out
wc -l tests/losat_out/long_test.out  # Compare with 14,877 (NCBI count)

# Short sequence test (regression)
./target/release/LOSAT tblastx \
    -q tests/fasta/MjeNMV.fasta \
    -s tests/fasta/MelaMJNV.fasta \
    -o tests/losat_out/short_test.out
wc -l tests/losat_out/short_test.out  # Should stay ~23,715
```
