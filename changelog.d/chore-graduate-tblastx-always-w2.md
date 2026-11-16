- The variables `LOSAT_X_DISPTABLE` and `LOSAT_X_TBXMASKSEEK` are gone: the display paths they chose were already
  the default and now run unconditionally, so `=0` no longer reaches the previous code. `DISPTABLE` is read by
  tblastx and tblastn (the displayed residues come from a per-code table of 4-bit codons); the scan of each
  codon's bases and the scan of every range of a frame's SEG masks that `=0` ran stay only as the references of
  unit tests. Outputs do not change.
