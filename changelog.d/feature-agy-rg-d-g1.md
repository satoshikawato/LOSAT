- `blastn -task megablast`: gap costs above 32767 (up to 2147483647) are no longer refused;
  affine greedy alignment widens its gap-cost, GCD, X-drop offset and distance arithmetic to 64
  bits and steps only over non-empty distances instead of allocating dense distance tables,
  matching NCBI BLAST+ 2.17.0 wherever NCBI completes without 32-bit signed overflow or dense-table
  memory exhaustion, and matching NCBI's nearby non-overflowing run where NCBI overflows `Int4`
  (proposed approved exception 31 of `PD-LOSAT-NCBI-DEFECTS`).
