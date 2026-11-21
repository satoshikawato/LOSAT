- An LMDB or lookup file of `-seqidlist`, `-negative_seqidlist`, `-taxids` or `-negative_taxids` (`.pdb`, `.ptf`,
  `.pto`, `.pot`, `.pos` and their nucleotide forms) whose bytes stop being readable after it was opened (a failing
  disk), where NCBI BLAST+ crashes, gives LOSAT's error `<path> cannot be read: <reason>` (exit status 1) instead of
  saying the file is not supported; an LMDB file whose first two pages cannot be read gives NCBI's `LMDB runtime
  error: mdb_env_open: <reason>` (exit status 2).
