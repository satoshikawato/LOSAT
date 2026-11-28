- Web engine (npm package `losat`, ABI version 17): `makeblastdb` runs in the browser with any command
  line that `losat makeblastdb` accepts (nucleotide and protein, database versions 4 and 5,
  `-parse_seqids`, `-hash_index`, volumes, `-mask_data`, `-gi_mask`, `-taxid`, `-taxid_map`, every
  `-input_type`, gzip input, `-logfile`, help and version), on the user's files, with the same exit
  status, standard output, standard error and database files as the command line run in a folder of
  those files (`WasmEngine.runTool(['makeblastdb', ...], { directory })`); the database is ready for a
  `-db` search in the same page. A tool's run directory is now writable, as the command line's working
  directory is: the files that the run makes or changes are kept in the Engine worker's memory and come
  to the run's file sinks (or `files`) by name when it ends, and the user's files are never changed. The
  engine's local time is the browser's time zone (a POSIX `TZ` rule in its environment), read at the
  start of each run, so the dates in a database are the user's local time, as on the command line, also
  after the user's zone changed since an earlier run. `-logfile` in the browser writes where the command
  line writes: `/dev/stdout` and `/dev/stderr` among the run's output, `/dev/null` nowhere, a file
  outside the run's folder (`/tmp/mk.log`) as a file of the run, and a name that does not open (a
  missing folder, a folder) leaves the messages on the standard error, in the check of the words too.
  LOSAT Web's form does not offer makeblastdb yet.
- `makeblastdb` on WASI (the Web engine): two `-gi_mask_name` values that name one GI-mask file
  (`m` and `./m`, or a folder through a symbolic link) are refused as on the command line, instead of
  writing both masks into one file; names are resolved through symbolic links as `realpath(3)` does
  where the platform has no `canonicalize`.
- LOSAT Web's engine reads a file that a run reads in order (makeblastdb's input FASTA, a database's
  files) from the browser in larger pieces, up to 1 MiB ahead: makeblastdb of a 146 MB genome FASTA
  makes 143 reads of the file instead of 4,474 and runs 2 to 4 times faster in Chromium.
