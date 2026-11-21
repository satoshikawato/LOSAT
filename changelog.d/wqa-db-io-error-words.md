- A database file that cannot be opened (no read permission, a directory in its place, too many open files) gives
  NCBI BLAST+'s error where LOSAT said the database "is not supported" (exit status 1): `LOSAT database error:
  Database memory map file error - CMemoryFile: Cannot memory map file '<path>':Permission denied` (exit status 2;
  with too many open files `Too many open files, please raise the open file limit`; the sequence file in the search:
  `LOSAT engine error: Database memory map file error`, exit status 3; a header file after the outfmt 0 prolog, as
  NCBI). LOSAT reads the first byte of each file where it opens it, its counterpart of NCBI's mapping of the whole
  file, and a first byte that cannot be read gives the text NCBI gives where that mapping fails (`Cannot map file
  '<path>' to memory (offset=0, length=<bytes>): <reason>`). The lookup files that `-taxids`, `-negative_taxids` and
  `-negative_seqidlist` open give NCBI's own errors: the taxonomy LMDB file (`.ptf`/`.ntf`) that cannot be opened or
  is empty, and the LMDB file (`.pdb`/`.ndb`) with too many open files, say `LMDB runtime error: mdb_env_open:
  <reason>` (exit status 2); an empty `.pos`/`.nos` or `.pot`/`.not` says `Failed to open oid-to-seqid lookup file`
  or `Failed to open oid-to-taxids lookup file` (exit status 2); a `.pos`, `.pot` or `.pto` that cannot be opened
  gives NCBI's memory map exception (exit status 255). blastdbcmd writes NCBI's text for a header, sequence or lookup
  file that it cannot open while it retrieves sequences (`Error: [blastdbcmd] CMemoryFile: Cannot memory map file
  '<path>':Permission denied`, without the search programs' `Database memory map file error - `). A volume's file
  (index, header, sequence, ISAM) whose bytes stop being readable after it was opened (a failing disk), where NCBI
  BLAST+ crashes, gives LOSAT's error `<path> cannot be read: <reason>` (exit status 1) instead of "not supported".
