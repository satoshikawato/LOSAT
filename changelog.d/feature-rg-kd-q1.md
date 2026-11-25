- `losat makeblastdb` whose log goes to a pipe that its reader has closed (`makeblastdb ... | head -0`)
  no longer builds a database nobody reads: it stops at the first line of the log that meets the
  closed pipe, as NCBI BLAST+ 2.17.0's `makeblastdb` does, and exits with status 0 and no message
  where NCBI's is ended by SIGPIPE (141). A reader that is gone when it starts leaves no new file and
  an existing database of the name as it was; a pipe that closes later stops the build at the next
  line of the log, and the files written before it stay. With `-logfile` naming a file the build
  runs to its end, as before.
