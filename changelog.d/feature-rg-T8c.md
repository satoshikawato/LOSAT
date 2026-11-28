- `losat makeprofiledb` that stops at a closed standard output (a reader such as `head` gone before the end) leaves
  the files NCBI's makeprofiledb 2.17.0 leaves there, byte for byte: the `.rps` and `.freq` files hold the profiles
  added before that point, and `.aux`, `.loo` and the volume's files hold only what NCBI's had written (before,
  LOSAT wrote out everything it had, so a stopped build looked more complete than NCBI's). LMDB's lock files of the
  volume being built stay, as NCBI's do; `losat makeblastdb` stopped by a closed pipe of its log now leaves them too.
- `losat blastdbcmd`, `dustmasker`, `segmasker`, `windowmasker`, `convert2blastmask`, `makeprofiledb` and
  `makembindex` with `-logfile /dev/stdout` or `-logfile /dev/stderr` on a pipe whose reader has gone now stop at the
  first message they write to that log, with exit status 0 and nothing more written, where NCBI's tools are ended by
  the signal there (a message of the configuration too, such as a malformed `-conffile`). Before, LOSAT lost the
  message and ran on to its usual end and exit status. `losat windowmasker -infmt seqids` runs on to exit status 1 as
  NCBI's does, whose network loader ignores the signal: before, its `maskinfo_*` output on a closed standard output
  ended it with 0.
- `losat blastdbcmd -outfmt %d` and `-outfmt %b` now pass each object to the output as soon as it is written, and a
  FASTA retrieval that stops at an error first passes on the records before it, as NCBI's blastdbcmd does: with
  `-logfile /dev/stdout` the messages now come between the entries where NCBI's do. Before, LOSAT passed the objects
  on a block at a time and put the error before the FASTA records.
