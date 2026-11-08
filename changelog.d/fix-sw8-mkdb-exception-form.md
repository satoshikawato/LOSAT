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
- Fixed: `losat makeblastdb` of version 5 with an output name of 247 bytes or more, whose LMDB lock
  file name would pass 255 bytes, now stops as makeblastdb does: `BLAST Database error: LMDB
  runtime error: mdb_env_open: File name too long` (exit status 2), or that creation error with
  one sequence (exit status 1), with no file left. It used to build the database at 247 and 248
  bytes and fail at the first volume file from 249 bytes on.
- Fixed: an LMDB environment variable of `losat makeblastdb` with a byte that is not UTF-8 is
  quoted in makeblastdb's error as makeblastdb quotes it (`'\377'`).
- Changed: `losat makeblastdb -in ""` is refused where makeblastdb crashes on it, after the head
  of the log, the removal of an existing database of the output name and the `-max_file_sz` and
  `-taxid_map` checks.
