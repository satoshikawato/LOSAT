- `blastn -task megablast`: `-gapopen` and `-gapextend` above 32767 (up to 2147483647), with or
  without `-no_greedy`, now run instead of being refused. Where NCBI BLAST+ 2.17.0 finishes, LOSAT
  writes the same output; with very large values where NCBI's 32-bit greedy arithmetic overflows or
  its dense tables run out of memory, LOSAT finishes in low memory using 64-bit arithmetic over
  only the reached distances and writes the alignments that NCBI BLAST+ finds without overflow
  (approved exception 32 (megablast greedy gap costs past `Int4`) of
  `docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md`).
