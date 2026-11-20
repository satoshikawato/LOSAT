- LOSAT Web: `-logfile` with a name that has a directory (`logs/run.log`, `/tmp/run.log`) runs as on the
  command line, where it was refused ("the -logfile name logs/run.log has a directory, which LOSAT Web does
  not have"). The file comes with the run under that whole name, as the command line makes it where the
  directory exists, and is saved with each `/` as `_` (`…logfile.logs_run.log`). `-logfile /dev/stderr` sends
  the messages to the run's Warnings, as `-logfile -` does, and a name of a directory (`logs/.`) does too, as
  NCBI sends a log it cannot open to the standard error. `-logfile /dev/stdout` gives a file of that name,
  apart from the outputs, where the command line writes its lines into its standard output among the
  report's. `-logfile logs/` gives the file `logs`, and `/dev/fd/2` sends the messages nowhere, as on the
  command line.
- `-logfile` names are read as NCBI BLAST+ 2.17.0 reads them, in every program: a name that ends in `/`
  makes the file without it where no such directory exists (`-logfile logs/` writes to `logs`, where the
  messages went to the standard error), the spaces at the ends of a relative name are dropped (` run.log`
  is `run.log`), `a/../run.log` is the `run.log` beside `a` also where `a` is a link, and a log on a file
  system that reports no free space gets no message (tblastx with `-logfile /dev/fd/2` or
  `/proc/self/fd/2` no longer writes its messages to the standard error, nor with `/dev/fd/1` into its
  report). tblastx with `-logfile /dev/full` writes its report, where it stopped with "No space left on
  device". makeblastdb writes its log to the name as typed, as NCBI's does, also where its messages go to
  another file (` mk.log`).
