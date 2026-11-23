- `losat blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, `psiblast`, `rpsblast`, `rpstblastn`,
  `makeblastdb` and `rpsbproc` read their arguments as NCBI BLAST+ 2.17.0 does, as the other programs
  already did. An argument error (an unknown option, a missing value, a value that does not convert or
  lies outside NCBI's limits such as `-word_size 3`, two options that exclude each other, a missing
  required option) now prints NCBI's short usage, a line of `=`, and NCBI's two error lines
  (`Error: <message>` and `Error:  <message>`; LOSAT writes the second without NCBI's
  `(CArgException::eX) ` part), and exits with status 1, where LOSAT printed its own message and
  exited with status 2. With `-logfile`, the error lines go to the log file, as NCBI's do.
- `losat <program> -help`, `-help-full` and `-h` print NCBI BLAST+ 2.17.0's help texts, in LOSAT's words
  (`Protein-Protein BLAST (LOSAT 0.2.0)`), where they printed LOSAT's own help. `-help-full` is `-help`,
  except for `makeblastdb`, whose `-help-full` also lists the hidden `-randomize`. `losat <program> --help`
  (two dashes) is no longer an option, as in NCBI: it ends with `Unknown argument: "-help"` and status 1;
  use `-help`. `losat -h`, `losat --help`, `losat --version` and an unknown program name are unchanged.
- An NCBI option that LOSAT does not run is refused after NCBI's argument check, so NCBI's own errors
  come first, in one form with status 1 (it was status 2):
  `Error: the NCBI BLAST+ option -html is not supported by LOSAT's BLASTP`. With `-dryrun`, such a
  command line ends with status 0, as NCBI's does.
- The library API and `pylosat`: an argument error (`Error::Argv`, Python `ArgumentError`) carries exactly
  what the command line writes to its standard error for the same arguments, and the refusal of an NCBI
  option that LOSAT does not run carries the refusal's sentence. `describe` takes the choices of an option
  from NCBI's argument tables (a flag has none) and its default from NCBI's own definition; `-evalue` of
  `blastx` and `tblastx` has no default there (their help says `[default: 10]`), and `-template_type`
  and `-template_length` list their choices.
- LOSAT Web: the message of an invalid command line is the command line's standard error (the usage, the
  line of `=` and NCBI's error lines) or the refusal's sentence, and the form marks the fields that the
  error lines name. LOSAT Web's ABI v1, which is frozen, keeps its argument messages.
