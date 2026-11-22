- Alias files with `OIDLIST` are read as NCBI BLAST+ 2.17.0 reads them (they were refused): the OID list
  file (`makeblastdb -oid_masks` writes such files as `.pxm` and `.nxm`) keeps the OIDs whose bits are set,
  in each volume below the alias file and after its `FIRST_OID` and `LAST_OID` (also where NCBI's
  combination with the ranges of alias files below it keeps more sequences than the list names); the
  search's statistics and the printed counts are NCBI's count of the kept sequences unless `NSEQ` and
  `LENGTH` are given, and a file that is missing or cannot be opened gives NCBI's error where NCBI meets
  it, in every search program, `blastdbcmd`, `segmasker`, `dustmasker` and `makembindex` (whose indexes, as NCBI's, hold every
  sequence of the volumes). A file that is empty or
  shorter than 5 bytes, or whose bits NCBI would read past the memory page after its end, is refused.
  An alias file whose `FIRST_OID`-`LAST_OID` range lies outside the OIDs that the ranges of the alias
  files below it keep, where NCBI writes outside its array of those OIDs, is refused where the OID list
  is first used.
  `blastdbcmd -info -exact_length` prints the exact total and longest length of the sequences an
  alias file's `OIDLIST` keeps, as NCBI does.
- `blastdbcmd` reports a database file that cannot be opened (for example one without read permission) in
  NCBI BLAST+ 2.17.0's words, `Error: [blastdbcmd] CMemoryFile: Cannot memory map file '<path>':<reason>`;
  it printed `Database memory map file error - ` before them.
