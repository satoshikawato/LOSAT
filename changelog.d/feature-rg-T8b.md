- `losat rpsbproc` handles a report it cannot write as NCBI's rpsbproc does: on a full disk the rest of the
  report is dropped and the run goes on to its end with NCBI's exit status (it stopped with `rpsbproc failed to
  write output` and exit status 4); when the reader of a pipe closes it, the run stops at once with exit status 0,
  after the progress lines that NCBI writes before that point.
