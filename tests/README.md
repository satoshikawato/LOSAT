# LOSAT Test Suite

This directory contains test scripts and utilities for comparing LOSAT output with NCBI BLAST+.

## Self-Comparison Test

The self-comparison test (`compare_self_tblastx.sh`) is designed to identify discrepancies between LOSAT and NCBI BLAST+ when comparing a sequence against itself.

### Known Issue (Fixed)

The original issue was that LOSAT TBLASTX would produce a single genome-wide hit for self-comparison (query == subject), while NCBI BLAST+ produces many smaller hits. This was caused by:

1. **Diagonal suppression (mask) not covering the full aligned region** - Fixed by using `se_ungapped` instead of `s_last_off - (k_size - 1)`
2. **Two-hit window too large** - Added `--ncbi_compat` flag to use NCBI standard value (16 vs LOSAT's 40)
3. **X-drop too permissive** - Added `--ncbi_compat` flag to use NCBI standard value (7 vs LOSAT's 11)

### Running the Test

```bash
# Build LOSAT first
cd LOSAT && cargo build --release && cd ..

# Run the self-comparison test
./tests/compare_self_tblastx.sh path/to/sequence.fna [genetic_code]

# Example with bacterial genome (genetic code 4)
./tests/compare_self_tblastx.sh NZ_CP006932.fna 4
```

### Analyzing Results

```bash
# Analyze all results in the output directory
python tests/analyze_blast_results.py --dir test_results_NZ_CP006932/

# Compare specific files
python tests/analyze_blast_results.py losat.tsv blast.tsv
```

### Expected Behavior

After the fixes, LOSAT should produce:
- Similar hit counts to NCBI BLAST+ (±20%)
- Similar alignment length distribution
- No single genome-wide hits for self-comparison

## Diagnostics

Enable diagnostics to get detailed pipeline statistics:

```bash
LOSAT_DIAGNOSTICS=1 ./LOSAT/target/release/losat tblastx -query seq.fna -subject seq.fna
```

This will output:
- Seed stage statistics (k-mer matches, two-hit filtering)
- Diagonal mask statistics (updates, suppressed seeds)
- Self-comparison detection
- Extension length statistics (average, max, warnings)
- HSP chaining and filtering statistics

## NCBI Compatibility Mode

Use the `--ncbi_compat` flag to match NCBI BLAST+ parameters more closely:

```bash
./LOSAT/target/release/losat tblastx \
    -query seq.fna \
    -subject seq.fna \
    --ncbi_compat
```

This changes:
- `TWO_HIT_WINDOW`: 40 → 16 (stricter)
- `X_DROP_UNGAPPED`: 11 → 7 (stricter)

These stricter parameters typically produce output that matches NCBI BLAST+ more closely, at the cost of potentially missing some marginal hits.

## Files

- `compare_self_tblastx.sh` - Shell script for running self-comparison tests
- `analyze_blast_results.py` - Python script for analyzing and comparing results
- `README.md` - This file

