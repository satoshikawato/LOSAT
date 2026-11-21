- An LMDB or lookup file of `-seqidlist`, `-negative_seqidlist`, `-taxids` or `-negative_taxids` (`.pdb`, `.ptf`,
  `.pto`, `.pot`, `.pos` and their nucleotide forms) whose bytes stop being readable after it was opened (a failing
  disk), where NCBI BLAST+ crashes, gives LOSAT's error `<path> cannot be read: <reason>` (exit status 1) instead of
  saying the file is not supported; an LMDB file whose first two pages cannot be read gives NCBI's `LMDB runtime
  error: mdb_env_open: <reason>` (exit status 2).
- A `taxdb.bti` or `taxdb.btd` on the search path that cannot be opened (no read permission): blastdbcmd gives NCBI
  BLAST+'s bytes (`-entry` FASTA: `Error: [blastdbcmd] CMemoryFile: Cannot memory map file '<path>':Permission denied`,
  exit status 1; `-outfmt` names `N/A`; `-tax_info` the exception, exit status 255), and a search whose output
  formats read nothing of the database sequences (such as `-outfmt "6 qseqid evalue bitscore"`) writes NCBI's rows
  after `Warning: [<program>] Error pre-fetching sequence data`. Every other search says that NCBI BLAST+ fails to
  load each database sequence with a taxid there and that LOSAT does not support it (exit status 1), where LOSAT
  said the taxonomy database is not supported.
