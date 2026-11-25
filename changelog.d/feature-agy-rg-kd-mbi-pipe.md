- `losat makembindex` stops at once with exit status 0 and no extra message when `-stat` (such as
  `-stat /dev/stdout` on a closed pipe or a pipe whose reader closes mid-run) or `-show_filters`
  meets a closed pipe outside Windows, matching where NCBI BLAST+ 2.17.0's `makembindex` is ended by
  SIGPIPE (141) without writing `done`, `.map`, `.shd`, or later volumes. Normal index builds (which
  write progress only to standard error) and non-pipe write failures (`/dev/full`) finish as NCBI's
  do.
