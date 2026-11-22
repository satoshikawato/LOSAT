- Web engine (npm package `losat`, ABI version 13): `windowmasker` runs in the browser with any command
  line that `losat windowmasker` accepts (`-mk_counts` in every `-sformat`, `-convert`, masking with
  `-ustat` in every `-outfmt`, with their options, help and version), on the user's files, with the
  same exit status, standard output, standard error and written files as the command line run in a
  folder of those files (`WasmEngine.runTool(argv, { directory })`), also when one file is named by
  both `-out` and `-logfile`. Its output, its messages and every file it writes (`-out FILE`,
  `-logfile FILE`) stream out as they are written, and the engine waits for a slow receiver on a
  cross-origin-isolated page, so a genome's counts or masked FASTA need not fit in memory;
  `validateTool(argv, { directory })` checks a windowmasker command line as the command line does
  before it reads its input. LOSAT Web's form does not offer windowmasker yet.
