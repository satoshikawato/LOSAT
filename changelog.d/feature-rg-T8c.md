- `losat makeprofiledb` that stops at a closed standard output (a reader such as `head` gone before the end) leaves
  the files NCBI's makeprofiledb 2.17.0 leaves there, byte for byte: the `.rps` and `.freq` files hold the profiles
  added before that point, and `.aux`, `.loo` and the volume's files hold only what NCBI's had written (before,
  LOSAT wrote out everything it had, so a stopped build looked more complete than NCBI's). LMDB's lock files of the
  volume being built stay, as NCBI's do; `losat makeblastdb` stopped by a closed pipe of its log now leaves them too.
