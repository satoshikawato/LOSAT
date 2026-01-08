# TBLASTX 2x Hit Issue - Debug Session Summary (Part 2)

## Date: January 2026

## Key New Findings

### 1. Score Difference Discovery
When comparing hits at identical coordinates, LOSAT bit scores are systematically **~7 bits higher** than NCBI:

```
Coord             LOSAT   NCBI    Diff
101249-101148     50.1    43.2    6.9
101435-101512     55.6    48.7    6.9
101745-101596     87.2    79.9    7.3
106400-105471     395     352     43
```

A 7-bit difference corresponds to ~15 raw score points (7 * ln(2) / 0.3176 ≈ 15.3).

### 2. Hit Distribution Analysis
For long sequences (AP027131 vs AP027133, ~600kb):
- Common hit coordinates: 11,233
- LOSAT-unique coordinates: 18,533 (extra hits)
- NCBI-unique coordinates: 3,638 (hits LOSAT misses)

Of the LOSAT-unique hits:
- ~93% have bit score < 50
- ~90% have E-value <= 0.01 (statistically significant)

### 3. Verified Correct Components
The following have been verified to match NCBI:
- BLOSUM62 matrix values
- Karlin-Altschul parameters (lambda=0.3176, K=0.134)
- X-drop value (16 raw score = 7 bits)
- Diagonal array sizing and offset management
- Two-hit window check (40 AA)
- Cutoff score calculation

### 4. Short Sequence Parity
For shorter sequences (AP027152 vs AP027155):
- LOSAT: 784 hits
- NCBI: 779 hits
- Ratio: 1.006x (near-parity)

This confirms the issue scales with sequence length.

## Hypotheses for Root Cause

### Hypothesis A: Extension Boundary Difference (Most Likely)
The systematic score difference suggests LOSAT alignments are **longer** than NCBI's.
Possible causes:
1. Different handling of sequence boundaries/sentinels
2. Different X-drop termination behavior
3. Different handling of the initial word score

### Hypothesis B: Score Accumulation Bug
A subtle difference in how scores are accumulated during extension.

### Hypothesis C: Reevaluation Step Missing
NCBI has a `Blast_HSPListReevaluateUngapped` step that recalculates scores.
This might truncate alignments differently.

## Next Investigation Steps

### Priority 1: Compare Individual HSP Details
Pick specific LOSAT-unique hits and trace the exact alignment:
```bash
# Example hit to trace:
# q: 101249-101148, s: 508751-508852, bit: 50.1 (LOSAT) vs 43.2 (NCBI)
```
- Compare alignment length
- Compare raw scores
- Verify residue-by-residue scoring

### Priority 2: Add Extension Tracing
Add detailed logging to `extend_hit_two_hit`:
- Input positions
- Initial score
- Max score achieved
- Final alignment boundaries
- Termination reason (score<=0, X-drop, or end of sequence)

### Priority 3: Investigate Reevaluation
Check if NCBI's reevaluation step (`Blast_HSPListReevaluateUngapped`) changes scores.

### Priority 4: Check Translation
Verify that translation is producing identical AA sequences.

## Test Commands

```bash
# Run with trace for specific HSP
LOSAT_TRACE_HSP="101249,101148,508751,508852" ./target/release/LOSAT tblastx \
    -q tests/fasta/AP027131.fasta -s tests/fasta/AP027133.fasta \
    -o /tmp/trace.out --query-gencode 4 --db-gencode 4 -n 1

# Compare outputs
diff <(cut -f7-12 /tmp/losat.out | sort) <(cut -f7-12 /tmp/ncbi.out | sort)
```

## Files Modified
None - this was an investigation session.

## Summary
The root cause appears to be in alignment scoring/boundaries rather than hit filtering.
LOSAT produces higher raw scores than NCBI for the same coordinates, leading to:
1. More HSPs passing the cutoff threshold
2. Higher bit scores in output
3. ~2x more final hits for long sequences
