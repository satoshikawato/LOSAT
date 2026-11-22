- `losat tblastx` accepts `-sum_stats` and `-max_intron_length`, which LOSAT refused before. With `-sum_stats false`
  each alignment keeps its own E-value instead of being linked into a group, and the pairwise report's description
  table has no `N` column (the alignments show `Expect` instead of `Expect(n)`), as in NCBI BLAST+ 2.17.0;
  `-sum_stats true` gives the default results. The TBLASTX page of LOSAT Web offers the option too.
  `-max_intron_length 0` changes nothing, and any other value stops with NCBI BLAST+'s error, which says that uneven
  gap linking is allowed for blastx, tblastn and psitblastn only.
