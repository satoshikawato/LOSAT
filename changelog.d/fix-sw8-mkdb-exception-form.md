- Fixed: `losat makeblastdb` now writes makeblastdb's exception report for the errors that
  makeblastdb ends with through its general handler (exit status 255): an output directory that
  cannot be made or written to, a `-max_file_sz` that is no size, and a `BLASTDB_LMDB_MAP_SIZE` or
  `MAX_LMDB_TXN_ENTRY` that is no number when the database gets a second sequence. The text is
  `Error: NCBI C++ Exception:` and the line of the exception, with makeblastdb's source file,
  line, class and function, where LOSAT used to print `Error:` and the sentence only.
- Fixed: `losat makeblastdb -max_file_sz` is read where makeblastdb reads it, after the head of
  the log, the removal of an existing database of the same name and the checks of the output
  directory: a wrong value now prints that head and removes the old database, and an output
  directory makeblastdb refuses is the error reported, as in makeblastdb.
- Fixed: a `BLASTDB_LMDB_MAP_SIZE` of more digits than a 64-bit number holds no longer adds
  `, overflow` to makeblastdb's message.
