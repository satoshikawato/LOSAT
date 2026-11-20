- LOSAT Web: `-logfile` with a name that has a directory (`logs/run.log`, `/tmp/run.log`) runs as on the
  command line, where it was refused ("the -logfile name logs/run.log has a directory, which LOSAT Web does
  not have"). The file comes with the run under that whole name, as the command line makes it where the
  directory exists, and is saved with each `/` as `_` (`…logfile.logs_run.log`). `-logfile /dev/stderr` sends
  the messages to the run's Warnings, as `-logfile -` does, and a name of a directory (`logs/.`) does too, as
  NCBI sends a log it cannot open to the standard error. `-logfile /dev/stdout` gives a file of that name,
  apart from the outputs, where the command line writes its lines into its standard output among the
  report's. `-logfile logs/` gives the file `logs`, and `/dev/fd/2` sends the messages nowhere, as on the
  command line. A name where a user's command line cannot make or open a file (`/etc/run.log`, `/run.log`,
  `/data/run.log`, `/dev/run.log`, `/dev/kmsg`, a name of more than 255 bytes) sends the messages to the run's
  Warnings, and `/dev/ptmx` nowhere, as there; a user makes files in the working directory, in a home
  (`/home/<user>/...`) and where every user does (`/tmp`, `/var/tmp`, `/dev/shm`, `/run/lock`, also as
  `/var/lock`, `/var/crash`). A name may climb out of a directory of root's, as there
  (`/etc/default/../../tmp/run.log` is `/tmp/run.log`), but not out of one that a user cannot enter
  (`/root/../tmp/run.log`), and a name of 4 KB starts a run at once.
- `-logfile` names are read as NCBI BLAST+ 2.17.0 reads them, in every program and on Windows by its
  Windows rules (`\` and `/` both end a directory there): a name that ends in `/`
  makes the file without it where no such directory exists (`-logfile logs/` writes to `logs`, where the
  messages went to the standard error), the spaces at the ends of a relative name are dropped (` run.log`
  is `run.log`), `a/../run.log` is the `run.log` beside `a` also where `a` is a link, and a log on a file
  system that reports no free space gets no message (tblastx with `-logfile /dev/fd/2` or
  `/proc/self/fd/2` no longer writes its messages to the standard error, nor with `/dev/fd/1` into its
  report). tblastx with `-logfile /dev/full` writes its report, where it stopped with "No space left on
  device". `-logfile //run.log` keeps the messages on the standard error, as NCBI's does (it wrote
  `/run.log` where `/` is writable). In a working directory of 4095 bytes or more a relative name is the
  file in that directory, as NCBI's (it was a shortened name elsewhere). makeblastdb writes its log to the
  name as typed, as NCBI's does, also where its messages go to another file (` mk.log`).
