- `losat rpsbproc` handles a report it cannot write as NCBI's rpsbproc does: on a full disk the rest of the
  report is dropped and the run goes on to its end with NCBI's exit status (it stopped with `rpsbproc failed to
  write output` and exit status 4); when the reader of a pipe closes it, the run stops at once with exit status 0,
  after the progress lines that NCBI writes before that point.
- `losat blast_formatter` writes its output through the stream that the search programs use, so a failed write ends as NCBI's
  `blast_formatter` 2.17.0 ends: a `-out` file that reaches its size limit in a `-outfmt 8` Seq-annot now posts NCBI's
  `Cannot close serializing output stream` and `BLAST failed to write output` with exit status 6 (it exited 0), and a standard
  output on a full device or a closed pipe after a standard error that could not be written exits with status 6 (it aborted).
  Through the library, a writer of yours that fails ends the call with `BLAST failed to write output` (exit status 6) instead of
  the abort text. The search programs' `-outfmt` 8, 9, 11, 12 and 16 on a full standard output after a failed standard error
  now exit with NCBI's status (0, or 6 for outfmt 16) instead of aborting.
- `losat makeblastdb` with its log on a pipe that closes before the sequences are added leaves the volume's files as NCBI's
  ended process leaves them (empty for a small input): the bytes its file buffers still held are dropped.
