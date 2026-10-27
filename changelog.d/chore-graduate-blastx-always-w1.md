- The variable `LOSAT_X_BXSCANBF` is gone: BLASTX's seed scans (the direct loop and `LOSAT_X_BXDIAG`'s
  fused loop) already found the subject words with lookup hits by tblastn's branch-free candidate loop by
  default, and now do so unconditionally. The scan on the presence test that `LOSAT_X_BXSCANBF=0` ran in
  the direct loop and the fused loop's presence-test cell mode were deleted from BLASTX (blastp keeps its
  own `LOSAT_X_BPSCANBF` switch and the shared cell mode); the compressed lookup table keeps its scan.
  `LOSAT_X_BXHITBF`'s test no longer has a `LOSAT_X_BXSCANBF=0` arm. Outputs do not change.
