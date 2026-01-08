# TBLASTX 2x Hit Issue - Debug Session Summary

## Problem Statement
LOSAT's TBLASTX produces approximately 2x the number of hits compared to NCBI BLAST for long sequences (~600kb), while shorter sequences produce near-parity.

## Current Test Results (January 2025)

### Long Sequences (AP027131 vs AP027133, ~600kb mollicute bacteria, gencode=4)
| Tool | Hit Count |
|------|-----------|
| NCBI | 14,871 |
| LOSAT | 29,766 |
| **Ratio** | **2.00x** |

### Short Sequences (AP027152 vs AP027155, ~280kb eukaryotic viruses, gencode=1)
| Tool | Hit Count |
|------|-----------|
| NCBI | 779 |
| LOSAT | 784 |
| **Ratio** | **1.01x** (near-parity) |

## Key Findings from Investigation

### Verified Correct (NOT the cause)
1. **Neighbor word generation** - Matches NCBI algorithm (triple nested loop with score threshold)
2. **Alphabet size** - Uses BLASTAA_SIZE=28 (matches NCBI)
3. **Threshold** - Default 13 (matches NCBI)
4. **Diagonal array sizing** - Power-of-2 based on qlen + window (matches NCBI formula)
5. **Two-hit detection** - Same logic including context boundary check
6. **Frame ordering** - +1,+2,+3,-1,-2,-3 (matches NCBI)
7. **Offset increment after frames** - `diag_offset += s_aa_len + window` (matches NCBI)
8. **Cutoff calculation**:
   - gap_trigger = 41 (from 22.0 bits)
   - cutoff_score_max = 64 (from E-value 10)
   - final_cutoff = 41 (MIN of above)
9. **Query context check** - `query_offset - diff < ctx.frame_base` (matches NCBI line 566-573)
10. **Frame base calculation** - `base += frame.aa_seq.len() - 1` (matches NCBI offset += length + 1)

### Diagnostic Data from Long Sequences (gencode=4)
```
=== Lookup Table Statistics ===
Total exact positions indexed: 829,398
Exact entries added (self_score < threshold): 159,741
Neighbor entries added: 6,900,411
Expansion factor: 8.51x
Non-empty cells: 13,668
Total entries: 7,060,152
Longest chain: (varies)

=== Extension Statistics ===
Total extensions: 17,391,241
Total linked HSPs: 129,945
Final output: 29,766 hits
```

### Characteristics of Extra LOSAT Hits (from earlier analysis)
- **97.8% have bit score < 50** (mostly 22-30 range)
- **Mean length: 26.9 AA** (short alignments)
- **Mean identity: 43.4%** (weak matches)
- **86.4% have E-value ≤ 0.001** (statistically significant)
- **LOSAT finds hits on 16,938 more diagonals than NCBI**
- **90.7% of extra hits on new diagonals have bitscore < 30**
- **All frame combinations show uniform ~2x ratio**

## Areas Still Under Investigation

### 1. Presence Vector Population
The number of non-empty lookup cells (13,668) and the expansion factor need to be compared directly with NCBI to verify k-mer indexing parity.

### 2. K-mer Match Counts
Need to verify that LOSAT finds the same number of initial seed hits (before two-hit filtering) as NCBI.

### 3. Two-Hit Window Passes
Need to count how many seed pairs pass the two-hit window check in LOSAT vs NCBI.

### 4. Extension Trigger Difference
LOSAT performs 17M extensions - need to verify NCBI performs a similar number.

## Hypothesis Status

### Hypothesis A: Missing HSP Overlap/Containment Filtering
**Status: NEEDS INVESTIGATION**
- LOSAT outputs 29,766 HSPs but has 129,945 "linked HSPs"
- Chain member filtering removes ~100k HSPs
- Need to verify if NCBI has additional overlap filtering

### Hypothesis B: Cutoff Score Difference
**Status: RULED OUT**
- Cutoff values match NCBI (gap_trigger=41, cutoff_score_max=64, final=41)
- Both accept the same minimum score threshold

### Hypothesis C: Diagonal Suppression
**Status: PARTIALLY INVESTIGATED**
- LOSAT finds 16,938 more diagonals with hits
- s_last_off calculation matches NCBI
- Issue may be upstream (more seeds passing two-hit check)

## Files of Interest

### LOSAT Source Files
- `LOSAT/src/algorithm/tblastx/utils.rs` - Main scan/extension logic (lines 1700-2100)
- `LOSAT/src/algorithm/tblastx/lookup.rs` - Lookup table construction
- `LOSAT/src/algorithm/tblastx/extension.rs` - Extension and s_last_off
- `LOSAT/src/algorithm/tblastx/ncbi_cutoffs.rs` - Cutoff calculations

### NCBI Source Files
- `c++/src/algo/blast/core/aa_ungapped.c` - s_BlastAaWordFinder_TwoHit (lines 440-614)
- `c++/src/algo/blast/core/blast_aascan.c` - s_BlastAaScanSubject
- `c++/src/algo/blast/core/blast_extend.c` - Diagonal table management
- `c++/src/algo/blast/core/blast_engine.c` - Subject frame loop

## Next Steps for Investigation

### Priority 1: Compare Initial Seed Counts
Add diagnostic to count total seeds found before any filtering:
```rust
// In s_blast_aa_scan_subject or the hit processing loop
// Count: total PV hits, total offset pairs generated
```

### Priority 2: Compare Two-Hit Pass Counts
Add counters for:
- Seeds that trigger two-hit extension (diff < window && diff >= wordsize)
- Seeds filtered by context boundary check
- Seeds filtered by "already extended past" check

### Priority 3: Compare Extension Success Rates
Count extensions that:
- Pass cutoff score check
- Are saved to init_hitlist

### Priority 4: Investigate HSP Filtering After Linking
The gap between 129,945 linked HSPs and 29,766 output suggests significant filtering.
- Verify chain member filtering logic
- Check if NCBI has additional HSP filtering

## How to Enable Diagnostics

```bash
# Enable detailed diagnostics
LOSAT_DIAGNOSTICS=1 ./target/release/LOSAT tblastx ...

# Enable cutoff debug
LOSAT_DEBUG_CUTOFFS=1 ./target/release/LOSAT tblastx ...

# Enable timing
LOSAT_TIMING=1 ./target/release/LOSAT tblastx ...
```

## Test Commands

```bash
# Long sequences with correct gencode (mollicute bacteria)
./target/release/LOSAT tblastx -q tests/fasta/AP027131.fasta -s tests/fasta/AP027133.fasta -o /tmp/losat.out --query-gencode 4 --db-gencode 4 -n 8

tblastx -query tests/fasta/AP027131.fasta -subject tests/fasta/AP027133.fasta -outfmt 6 -db_gencode 4 -query_gencode 4 -num_threads 8 > /tmp/ncbi.out

# Short sequences (eukaryotic viruses)
./target/release/LOSAT tblastx -q tests/fasta/AP027152.fasta -s tests/fasta/AP027155.fasta -o /tmp/losat_short.out --query-gencode 1 --db-gencode 1 -n 8

tblastx -query tests/fasta/AP027152.fasta -subject tests/fasta/AP027155.fasta -outfmt 6 -db_gencode 1 -query_gencode 1 -num_threads 8 > /tmp/ncbi_short.out
```

## Key Code Sections to Review

### LOSAT Two-Hit Check (utils.rs:1830-1870)
```rust
let last_hit = diag_entry.last_hit - diag_offset;
let diff = subject_offset - last_hit;

if diff >= window { ... continue; }
if diff < wordsize { ... continue; }

// Context boundary check (NCBI parity)
if query_offset - diff < ctx.frame_base {
    diag_entry.last_hit = subject_offset + diag_offset;
    continue;
}

// Extension triggered here
```

### NCBI Two-Hit Check (aa_ungapped.c:533-573)
```c
last_hit = diag_array[diag_coord].last_hit - diag_offset;
diff = subject_offset - last_hit;

if (diff >= window) { ... continue; }
if (diff < wordsize) { continue; }

// Context boundary check
if (query_offset - diff < query_info->contexts[curr_context].query_offset) {
    diag_array[diag_coord].last_hit = subject_offset + diag_offset;
    continue;
}

// Extension triggered here
```

## Summary
The 2x hit issue persists for long sequences with gencode=4. The root cause is likely in the seeding/scanning phase where LOSAT finds seeds on ~17k more diagonals than NCBI. All individual code components have been verified to match NCBI, so the issue may be in the interaction between components or in subtle differences in how data flows through the pipeline. The next session should focus on adding detailed counters at each stage to pinpoint where the divergence occurs.
