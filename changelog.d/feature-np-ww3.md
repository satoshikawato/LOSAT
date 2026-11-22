- Web engine (npm package `losat`, ABI version 13): `windowmasker` runs in the browser with any command
  line that `losat windowmasker` accepts (`-mk_counts` in every `-sformat`, `-convert`, masking with
  `-ustat` in every `-outfmt`, with their options, help and version), on the user's files, with the
  same exit status, standard output, standard error and written files as the command line run in a
  folder of those files (`Engine.runTool(argv, { directory })`). Its output and every file it writes
  (`-out FILE`, `-logfile FILE`) stream out as they are written, so a genome's counts or masked FASTA
  need not fit in memory; `validate` checks a windowmasker command line as the command line does
  before it opens a file. LOSAT Web's form does not offer windowmasker yet.
