- BLAST database alias files with `MEMB_BIT` are now read as NCBI BLAST+ 2.17.0 reads them in
  `blastdbcmd`, `blastdbcheck` and the search programs (`blastp`, `blastn`, `blastx`, `tblastn`,
  `tblastx`, `psiblast`, `rpsblast`, `rpstblastn`): definition lines whose first membership bitmask
  lacks the bit are hidden, a value of `0` disables filtering, and conflicting bits on the same
  volume stop with NCBI's error. `MEMB_BIT` above 32 is refused where a definition line carries
  membership bits, and an RPS profile left without definition lines in a multi-volume database is
  refused where it would be reported. A database with GI-indexed masks (`MASKLIST`) now reports an
  alias file's OID list error at its first mask read, as NCBI does.
