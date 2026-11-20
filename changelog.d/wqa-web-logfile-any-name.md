- LOSAT Web: `-logfile` with a name that has a directory (`logs/run.log`, `/tmp/run.log`) runs as on the
  command line, where it was refused ("the -logfile name logs/run.log has a directory, which LOSAT Web does
  not have"). The file comes with the run under that whole name, as the command line makes it where the
  directory exists, and is saved with each `/` as `_` (`…logfile.logs_run.log`). `-logfile /dev/stderr` sends
  the messages to the run's Warnings, as `-logfile -` does, and a name of a directory (`logs/`) does too, as
  NCBI sends a log it cannot open to the standard error. `-logfile /dev/stdout` gives a file of that name,
  apart from the outputs, where the command line writes its lines into its standard output among the
  report's.
