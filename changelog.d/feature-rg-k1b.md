- Library API and pylosat: `run_dustmasker` (`pylosat.dustmasker`) runs `dustmasker` with the
  command line's words, masking sequences given in memory as its standard input. `run_blastdbcmd`
  (`pylosat.blastdbcmd`) now takes the standard input too (`-entry_batch -`, `-ipg_batch -`,
  `-taxidlist -`): a Rust caller passes a reader as its new second argument, and neither function
  reads the process's standard input any more.
- `losat blastdbcmd`: after a failed write of the retrieved sequences (a full disk) the retrieval
  goes on, as NCBI's does, so an id that is not found later is still reported (`Skipped`, exit
  status 1); a failed write of `-metadata` ends with an error instead of exit status 0.
- `losat blastdbcmd` and `losat dustmasker` stop at once when the reader of their output closes
  the pipe (`| head`), with exit status 0 and no message, where NCBI's programs are ended by
  SIGPIPE. blastdbcmd writes its error lines after the output written before them, as NCBI does
  (with `2>&1` the order is NCBI's).
- `losat dustmasker -infmt blastdb` reads the alias files that only a search reads differently
  (`LENGTH 0`, a negative `NSEQ`, `FIRST_OID` of 2^31 or more, an OID range with another name), as
  NCBI's dustmasker does, instead of refusing them.
