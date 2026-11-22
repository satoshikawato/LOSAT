- An LMDB or lookup file of `-seqidlist`, `-negative_seqidlist`, `-taxids` or `-negative_taxids` (`.pdb`, `.ptf`,
  `.pto`, `.pot`, `.pos` and their nucleotide forms) whose bytes stop being readable after it was opened (a failing
  disk), where NCBI BLAST+ crashes, gives LOSAT's error `<path> cannot be read: <reason>` (exit status 1) instead of
  saying the file is not supported; an LMDB file whose first two pages cannot be read gives NCBI's `LMDB runtime
  error: mdb_env_open: <reason>` (exit status 2).
- A `taxdb.bti` or `taxdb.btd` on the search path that cannot be opened (no read permission): blastdbcmd gives NCBI
  BLAST+'s bytes (`-entry` FASTA: `Error: [blastdbcmd] CMemoryFile: Cannot memory map file '<path>':Permission denied`,
  exit status 1; `-outfmt` names `N/A`; `-tax_info` the exception, exit status 255), and a search writes NCBI's output
  after `Warning: [<program>] Error pre-fetching sequence data` where its output formats need no database sequence's
  ids, titles, taxonomy or letters (tables of fields such as `qseqid`, `slen`, `sstart`, `length` or `evalue`, the
  organism report `-outfmt 18`, which then lists no database sequence with a taxid and is unchanged with `-subject`,
  and the Seq-align formats).
  psiblast's PSSM then gets no letters of an included sequence with a taxid, warns `Failed to retrieve sequence <id>`
  and uses the query's letters in that row, as NCBI BLAST+ does. Every other search stops with exit status 1 and says
  to make the file readable or to move it off the BLASTDB path, where LOSAT said the taxonomy database is not
  supported.
- A `taxonomy4blast.sqlite3` whose pages stop being readable while `-taxids` or `-negative_taxids` is expanded gives
  NCBI BLAST+'s messages in LOSAT's words: its SQLite error (`disk I/O error`), then `LOSAT query/options error: Invalid taxidlist
  file` (exit status 1), or in blastdbcmd the exception (exit status 255).
