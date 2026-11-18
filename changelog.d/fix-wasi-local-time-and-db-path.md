- WASI command builds (`wasm32-wasip1`, `wasm32-wasip1-threads`): makeblastdb's log time
  (`Building a new DB, current time:`) and the database's date (blastdbcmd `-info`'s `Date:`, a
  report's `Posted date:`) are the local time, as natively, where they were UTC: the command reads
  `TZ`, `TZDIR` and `/etc/localtime` as glibc does, from the environment and the directories its
  host gives it. `New DB name:` and `Deleted existing ... named` print the absolute path, as
  natively, where they printed the `-out` name as given: the name is made absolute with the working
  directory that the command takes from the host's `PWD`. Native builds give the same bytes as
  before, except for a `-out` name or a working directory with C spaces (a tab, a newline) at its
  ends, which makeblastdb now cuts as NCBI's `ConcatPath` does. The WASI command checks
  (`LOSAT/tests/check_wasm_threading.py`) run makeblastdb, blastdbcmd `-info` and a blastn report in
  a working directory with `TZ` a zone file, a POSIX rule and unset through each WASI command.
