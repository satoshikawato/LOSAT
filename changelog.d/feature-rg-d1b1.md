- Alias files with `OIDLIST` are read as NCBI BLAST+ 2.17.0 reads them (they were refused): the OID list
  file (`makeblastdb -oid_masks` writes such files as `.pxm` and `.nxm`) keeps the OIDs whose bits are set,
  in each volume below the alias file and after its `FIRST_OID` and `LAST_OID`; the search's statistics
  and the printed counts are NCBI's count of the kept sequences unless `NSEQ` and `LENGTH` are given,
  and a missing file gives NCBI's error where NCBI meets it, in every search program, `blastdbcmd`,
  `segmasker` and `dustmasker`. An empty file, or one of its 4-byte header alone, is refused.
  `blastdbcmd -info -exact_length` prints the exact total and longest length of the sequences an
  alias file's `OIDLIST` keeps, as NCBI does.
