- The library API and pylosat run the NCBI C++ Toolkit options `-logfile`, `-conffile` and `-dryrun` as the command line
  does: the `-logfile` file gets what `losat <program> ... -logfile FILE` writes there, `-conffile` is loaded in place of
  `<program>.ini`, and a `-dryrun` call checks its arguments and settings, then searches nothing and reads no record
  (pylosat: `Search.dry_run`). They were refused.
- The library's makeblastdb, blastdbcmd, dustmasker, segmasker and convert2blastmask run the same start: `-conffile` and
  `-dryrun` work there, and `run_makeblastdb` prints the versions to its output and writes its refusals to its
  diagnostics, as the other tools do.
- The library's database searches, and its searches with outfmt 18 or blastn's `-window_masker_taxid`, write NCBI's
  messages about a registry file (such as a malformed `.ncbirc`) and go on, or stop with NCBI's text, as the command line
  does. They refused the file.
