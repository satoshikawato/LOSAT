- TBLASTX's sum-statistics HSP linking follows NCBI BLAST+'s 16-bit set count in both of its
  linking kernels: the count of a set of more than 32,767 HSPs wraps as NCBI's does (the reference
  kernel, used by `LOSAT_LINKING_LEGACY` and the linking traces, stopped with an overflow in builds
  with overflow checks), the E-value of such a set is compared as NCBI's signed count, and the
  large-gap-only linking of a search space too small for small gaps is NCBI's statement without the
  division by `1 - gap_prob` (LOSAT divided by 1 there). Outputs do not change.
