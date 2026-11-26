- CLI wording: when an NCBI-style exception has an empty message body (such as
  `GenCountsException::empty genome` in `losat windowmasker -mk_counts`), keep the error code text
  (`empty genome`) as the message instead of printing a blank `Error: ` line.
