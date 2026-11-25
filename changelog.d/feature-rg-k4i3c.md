- The native `losat` command line and the pylosat wheel no longer contain clap and the crates it
  brings (11 crates in all on Linux and macOS, 13 on Windows), and no build compiles clap's
  derive: the ten programs whose options NCBI's argument reader reads (`blastn`, `blastp`,
  `blastx`, `tblastn`, `tblastx`, `psiblast`, `rpsblast`, `rpstblastn`, `makeblastdb`,
  `rpsbproc`) now take their options straight from what the reader read, without a second
  conversion. What a run writes is unchanged (output, messages, exit statuses, defaults, help;
  the top level's help, `losat` or `losat -h`, still names the executable as it was called:
  `losat.exe` on Windows, or the name of a copy or a link). Only a diagnostic build's startup trace
  (`LOSAT_STARTUP_TRACE=1`) now says `[startup] after the options are read`.
  The `losat` binary for Linux x86-64 is 0.75 MB (2.4 %) smaller, a tiny `losat blastn` run
  executes 0.27 million (2.5 %) fewer instructions, and a program that checks two command lines
  with the library's `validate` executes half as many; a small gain, well under a millisecond per
  run. LOSAT Web is unchanged: its frozen ABI v1 keeps its own copy of clap's command for the
  messages of its BLASTX and TBLASTN.
