- `losat makembindex` stops at once with exit status 0 and no extra message when `-stat` (such as
  `-stat /dev/stdout` on a closed pipe or a pipe whose reader closes mid-run), `-show_filters`, or
  an index output meets a closed pipe outside Windows, matching where NCBI BLAST+ 2.17.0's
  `makembindex` is ended by SIGPIPE (141) without writing `done`, `.map`, `.shd`, or later volumes
  and leaving the incomplete `.NN.idx` with the same 44-byte header and complete 8 KiB blocks.
  Normal index builds (which write progress only to standard error) and non-pipe write failures
  (`/dev/full`) finish as NCBI's do.
