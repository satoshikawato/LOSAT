# TBLASTX 2x Hit Count Investigation

## Problem Statement

LOSAT produces approximately twice as many alignment hits as NCBI BLAST for long sequences (~600kb).

| Metric | LOSAT (before) | LOSAT (after fix) | NCBI | Ratio |
|--------|----------------|-------------------|------|-------|
| Hit count | 29,766 | 28,819 | 14,877 | 1.94x |
| Test case | AP027131 vs AP027133 (gencode 4) | | | |
| Sequence lengths | 671kb vs 615kb | | | |

**Key observation**: Short sequences have near-parity (ratio ~1.006x). The issue scales with sequence length.

---

## BUGS FOUND AND FIXED

### 1. Neighbor Word Expansion - Extra Entries for Low-Scoring Words (2026-01-08)

**File:** `lookup.rs:431-438, 481-489`

**Bug:** LOSAT was adding exact match entries to lookup table when `self_score < threshold`. NCBI does NOT do this - it only adds entries through the neighbor loop when `score >= threshold`.

**Impact:** 126,185 extra lookup entries were added, generating ~1M extra k-mer matches.

**Fix:** Removed the incorrect `self_score < threshold` condition. Now only adds entries for `threshold == 0` (exact-only matching mode).

**Result:** Hits reduced from 29,766 to 28,819 (reduction of 947 hits).

```rust
// BEFORE (WRONG):
if threshold == 0 || self_score < threshold {
    entry_counts[idx] += num_offsets;  // Added 126k extra entries!
}

// AFTER (CORRECT - matches NCBI blast_aalookup.c s_AddWordHitsCore):
if threshold == 0 {
    entry_counts[idx] += num_offsets;
    continue;
}
```

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

## Next Investigation Steps

1. **Compare k-mer match counts** - 1.1B matches seems high; verify if NCBI has similar count
2. **Single position trace** - Follow one specific hit through the entire pipeline
3. **Check for duplicate offsets** - Verify no offset is indexed twice in lookup table

---

*Last updated: 2026-01-08 (Session 3)*
