- The variables `LOSAT_X_BNLUTDROP`, `LOSAT_X_BNUNGAPCOUNT` and `LOSAT_X_CTXO1` are gone: the blastn paths they
  chose were already the default and now run unconditionally, so `=0` no longer reaches the previous code. What
  went is that previous code: the query batch's lookup tables held until the batch returns, the loop-for-loop
  exact ungapped extension (now a test reference only) and the query-context lookup without the constant-time form (the binary
  search over the contexts still answers where the layout does not allow the constant-time form). Outputs do not
  change.
