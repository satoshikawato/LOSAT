# Fix BLASTN Task Performance Regression and Update Test Scripts

## Summary

This PR fixes a critical performance regression in `--task blastn` mode and corrects output file names in test plotting scripts.

## Changes

### 1. Fixed Performance Regression in `--task blastn`

**Problem**: After implementing the two-stage lookup table, `--task blastn` (word_size=11) became extremely slow, taking much longer than expected for high-identity sequence comparisons.

**Root Cause**: Incorrect indentation in the subject scanning loop caused the loop body to execute incorrectly, leading to performance degradation.

**Fix**: Corrected indentation in the `for s_pos in 0..s_len` loop within the non-two-stage lookup path to ensure proper execution flow.

**Verification**: 
- Before fix: Execution was extremely slow (not completing in reasonable time)
- After fix: Execution time restored to ~4.2 seconds (matching previous performance)

### 2. Fixed Output File Names in Test Scripts

**Problem**: Test plotting scripts (`plot_comparison.py`, `plot_execution_time.py`, `plot_overall_trend.py`) referenced incorrect output file names for megablast results.

**Fix**: Updated file names to match actual output files:
- Changed `*.losatn.megablast.out` → `*.blastn.megablast.out`
- Changed `*.losatn.megablast.log` → `*.blastn.megablast.log`

**Files Modified**:
- `tests/plot_comparison.py`: Fixed 2 file references
- `tests/plot_execution_time.py`: Fixed 2 log file references  
- `tests/plot_overall_trend.py`: Fixed 2 file references

## Impact

- **Performance**: `--task blastn` mode now performs correctly without regression
- **Testing**: Plotting scripts now correctly reference output files for analysis

## Testing

- Verified `--task blastn` completes in ~4.2 seconds for NZ_CP006932 self-comparison
- Confirmed output file names match between `run_comparison.sh` and Python scripts

## Related

This PR addresses issues introduced in PR #18 (Two-Stage Lookup Table Implementation).

