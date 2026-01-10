# BLASTN (LOSATN) Bug Fix Manual

**Purpose**: This document provides Claude with task-specific guidance for debugging and fixing BLASTN issues in LOSAT.

---

## Quick Reference

### Build & Test Commands
```bash
cd LOSAT && cargo build --release
./target/release/LOSAT blastn -q query.fasta -s subject.fasta -o out.txt --task blastn
./target/release/LOSAT blastn -q query.fasta -s subject.fasta -o out.txt  # megablast (default)
```

### Debug Environment Variables
```bash
LOSAT_DEBUG_BLASTN=1      # Hit loss statistics per subject
LOSAT_DEBUG_COORDS=1      # Coordinate transformations
LOSAT_DEBUG_CUTOFFS=1     # Cutoff score calculations
LOSAT_DIAGNOSTICS=1       # General pipeline diagnostics
BLEMIR_DEBUG=1            # Build config and settings
```

### Key Test Files
```
LOSAT/tests/fasta/NZ_CP006932.fasta     # Large self-comparison test
LOSAT/tests/blast_out/*.task_blastn.out # NCBI reference output
LOSAT/tests/losat_out/*.losatn.*.out    # LOSAT output for comparison
```

---

## BLASTN Architecture Overview

### Execution Flow
```
1. Build lookup table (two-stage, 8-bit LUT → full word verification)
2. Scan subject sequence (rolling k-mer extraction)
3. For each seed hit:
   a. Verify full word match (extend to word_size)
   b. Ungapped extension (X-drop termination)
   c. Check cutoff score
   d. Gapped extension (DP or Greedy)
   e. Containment filtering
4. Post-processing:
   a. Endpoint purging (remove duplicates)
   b. E-value calculation
   c. Output formatting
```

### Key Files
| File | Responsibility |
|------|----------------|
| `blast_engine/run.rs` | Main execution loop, parallel search |
| `alignment/gapped.rs` | DP-based gapped extension (blastn task) |
| `alignment/greedy.rs` | Greedy alignment (megablast task) |
| `extension.rs` | Ungapped extension, X-drop |
| `lookup.rs` | Two-stage lookup table construction |
| `ncbi_cutoffs.rs` | Cutoff score calculations |
| `filtering/purge_endpoints.rs` | HSP endpoint purging |

---

## Task Differences: Megablast vs Blastn

| Parameter | Megablast | Blastn |
|-----------|-----------|--------|
| Word Size | 28 | 11 |
| Reward/Penalty | 1/-2 | 2/-3 |
| Gap Penalties | 0,0 (none) | 5,2 |
| Gapped Algorithm | Greedy | DP (Smith-Waterman) |
| X-drop Gapped | 25 | 30 |
| Use Case | High identity (>90%) | Divergent sequences |

---

## Critical Implementation Details

### DO: Follow These Practices

1. **Always include NCBI reference comments**
   ```rust
   // NCBI reference: blast_gapalign.c:735-962
   // Original C code:
   // best = row[b_index].best;
   // best_gap = row[b_index].best_gap;
   ```

2. **Use the correct Karlin parameters**
   - **UNGAPPED params (kbp_std)**: For `gap_trigger` calculation
   - **GAPPED params (kbp_gap)**: For `cutoff_score_max` calculation
   - Reference: `blast_parameters.c:340-344`

3. **Match NCBI's coordinate systems exactly**
   - Extension phase: 0-indexed, frame-relative
   - Output phase: 1-indexed, nucleotide positions
   - Subject coordinates: `offset = min(start, end)`, `end = max(start, end)`

4. **Preserve NCBI's floating-point precision**
   ```rust
   // Use ceil/floor as NCBI does
   let gap_trigger = (gap_trigger_bits * NCBIMATH_LN2 + kbp.logK) / kbp.lambda;
   let gap_trigger_raw = gap_trigger.ceil() as i32;
   ```

5. **Keep sentinel byte handling**
   - SENTINEL_BYTE = 0
   - Query/subject sequences have sentinel bytes at boundaries
   - Extension must terminate at sentinels

### DON'T: Avoid These Mistakes

1. **DON'T simplify NCBI algorithms for readability**
   - NCBI's exact logic is the source of truth
   - If NCBI has a weird condition, replicate it exactly

2. **DON'T change post-loop comparison behavior**
   ```rust
   // WRONG: Extra comparison after loop (causes over-extension)
   if best_score > max_score {
       max_score = best_score;
       // This extends alignment beyond NCBI's boundary
   }

   // RIGHT: Match NCBI exactly - no post-loop comparison
   // NCBI's loop already handles this internally
   ```

3. **DON'T mix up band contraction values**
   ```rust
   // WRONG: Using +2 for band contraction
   last_b_index = last_b_index + 2;

   // RIGHT: NCBI uses +1
   // Reference: blast_gapalign.c:933
   last_b_index = last_b_index + 1;
   ```

4. **DON'T add non-NCBI filtering**
   - No "max hits per k-mer" filtering
   - No diagonal masking beyond `hit_level_array`
   - No post-processing containment filter (NCBI does it during alignment)

5. **DON'T guess at edge cases**
   - Check NCBI source code for boundary conditions
   - Empty sequences, single-base alignments, etc.

---

## Known Issues and Status

### Issue 1: Missing Hits in Self-Comparison (Priority: High)
**Symptom**: NZ_CP006932 shows 72% coverage (8,843 vs 12,340 NCBI hits)

**Root Cause**: Missing traceback information for endpoint trimming
- LOSAT generates 13,602 raw hits
- Endpoint purging removes 35% (to 8,843)
- NCBI uses traceback to TRIM overlapping HSPs, not delete them

**Current Status**: Partially implemented
- `purge_endpoints.rs` has trimming infrastructure
- `Hit.gap_info` exists but is always `None` (no traceback capture)

**To Fix**:
1. Add traceback capture to `extend_gapped_one_direction()` in `gapped.rs`
2. Populate `Hit.gap_info` with `GapEditOp` sequence
3. Enable 2-pass endpoint purge (trim → re-evaluate → delete)

**NCBI References**:
- `blast_hits.c:2392-2452` (s_CutOffGapEditScript)
- `blast_traceback.c:637-669` (three-phase purge)

### Issue 2: Extra Hits in Some Datasets (Priority: Medium)
**Symptom**: Some datasets show 105-116% of NCBI hits

**Root Cause**: Unknown - possibly:
- Different seed finding behavior
- Different ungapped extension termination
- Different cutoff calculations for short alignments

**To Investigate**:
1. Compare seed positions between NCBI and LOSAT
2. Check X-drop termination for short alignments
3. Verify cutoff score calculations match NCBI

---

## Debugging Workflow

### Step 1: Reproduce the Issue
```bash
cd LOSAT

# Build with release optimizations (debug builds are too slow)
cargo build --release

# Run with debug output
LOSAT_DEBUG_BLASTN=1 ./target/release/LOSAT blastn \
    -q tests/fasta/NZ_CP006932.fasta \
    -s tests/fasta/NZ_CP006932.fasta \
    -o /tmp/test.out --task blastn -n 1 2>&1 | tee debug.log
```

### Step 2: Compare with NCBI
```bash
# Count hits
wc -l /tmp/test.out
grep -v "^#" tests/blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out | wc -l

# Compare top hits (should match positions and scores)
head -10 /tmp/test.out
grep -v "^#" tests/blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out | head -10
```

### Step 3: Trace Specific Coordinates
```bash
# Enable coordinate tracing
LOSAT_DEBUG_COORDS=1 ./target/release/LOSAT blastn \
    -q tests/fasta/query.fasta -s tests/fasta/subject.fasta \
    -o /tmp/test.out --task blastn 2>&1 | grep COORDS
```

### Step 4: Check Cutoff Calculations
```bash
LOSAT_DEBUG_CUTOFFS=1 ./target/release/LOSAT blastn \
    -q tests/fasta/query.fasta -s tests/fasta/subject.fasta \
    -o /tmp/test.out --task blastn 2>&1 | grep CUTOFF
```

---

## NCBI Source Code Locations

### Primary References
```
/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/

blast_gapalign.c      # Gapped alignment (most important)
  - Blast_SemiGappedAlign: lines 735-962
  - ALIGN_EX: lines 364-733

blast_parameters.c    # Cutoff calculations
  - gap_trigger: lines 340-344
  - cutoff_score_max: lines 348-367

blast_hits.c          # HSP filtering
  - PurgeHSPsWithCommonEndpoints: lines 2455-2535
  - s_CutOffGapEditScript: lines 2392-2452
  - ReevaluateWithAmbiguities: lines 479-647

na_ungapped.c         # Nucleotide extension
  - Word verification: lines 508-607
  - Ungapped extension: lines 220-233

blast_stat.c          # Karlin parameters
blast_itree.c         # Interval tree containment
greedy_align.c        # Greedy alignment (megablast)
```

### Finding NCBI Code
```bash
# Search for function
grep -rn "Blast_SemiGappedAlign" /mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/

# Search for constant
grep -rn "BLAST_GAP_TRIGGER_NUCL" /mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/include/
```

---

## Common Bug Patterns

### Pattern 1: Off-by-One in Coordinates
**Symptom**: Alignment starts/ends at wrong position
**Check**: Compare `q_consumed` and `s_consumed` calculations
**NCBI Reference**: `blast_gapalign.c:933-952`

### Pattern 2: Wrong Karlin Parameters
**Symptom**: Cutoff score too high/low, wrong E-values
**Check**: Verify `kbp_std` (ungapped) vs `kbp_gap` (gapped) usage
**NCBI Reference**: `blast_parameters.c:340-374`

### Pattern 3: Band Contraction Issues
**Symptom**: Alignments truncated or over-extended
**Check**: `last_b_index` updates in DP loop
**NCBI Reference**: `blast_gapalign.c:858-876`

### Pattern 4: Endpoint Purge Too Aggressive
**Symptom**: Missing valid HSPs after filtering
**Check**: Sort order and comparison criteria in `purge_endpoints.rs`
**NCBI Reference**: `blast_hits.c:2310-2384`

---

## Testing After Fixes

### Quick Validation
```bash
# Build and test
cd LOSAT && cargo build --release

# Run self-comparison
./target/release/LOSAT blastn -q tests/fasta/NZ_CP006932.fasta \
    -s tests/fasta/NZ_CP006932.fasta -o /tmp/test.out --task blastn -n 1

# Compare hit count
echo "LOSAT: $(wc -l < /tmp/test.out)"
echo "NCBI:  $(grep -v '^#' tests/blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out | wc -l)"
```

### Full Test Suite
```bash
cd LOSAT/tests && bash compare_results.sh
```

### Specific Dataset Tests
```bash
# Megablast (should be ~99% coverage)
./target/release/LOSAT blastn -q tests/fasta/EDL933.fasta \
    -s tests/fasta/Sakai.fasta -o /tmp/mega.out

# Blastn with extra hits (should be ~100-115% coverage)
./target/release/LOSAT blastn --task blastn -q tests/fasta/MelaMJNV.fasta \
    -s tests/fasta/PemoMJNVA.fasta -o /tmp/blastn.out
```

---

## Performance Notes

- **Release builds only**: Debug builds are 10-50x slower
- **Single thread for debugging**: Use `-n 1` to avoid race conditions in output
- **Memory**: Large sequences (600kb+) need ~2GB RAM
- **Greedy alignment**: Uses thread-local memory pool (GREEDY_MEM)

---

## Summary: What to Remember

1. **NCBI is truth**: Match their output bit-perfectly
2. **Include NCBI references**: Every algorithm change needs C code comments
3. **Test with multiple datasets**: Both self-comparison and divergent sequences
4. **Check both megablast and blastn tasks**: They use different algorithms
5. **Debug incrementally**: Use environment variables to trace issues
6. **Don't over-engineer**: Only fix what's broken, don't refactor "for readability"
