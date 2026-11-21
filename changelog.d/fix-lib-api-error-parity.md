- The library API and pylosat fail as the command line does in more cases. A report in a pairwise format (outfmt 0
  to 4) or the organism report (outfmt 18) that `run_to_file` (pylosat `Search.run_to_file`) cannot write, such as to
  a full disk, fails with `BLAST failed to write output` (the command line's exit status 6), where it gave the
  system's error (`No space left on device`); a short organism report that the command line writes without an error
  no longer fails. A directory given as the query or subject file (`read_input_file`, pylosat `Sequences.from_fasta`)
  is read as an empty input, where it was refused with `Is a directory`. `run_makeblastdb` writes the message of an
  argument that its parser refuses to its diagnostics too, as the command line writes it to its standard error.
- Every search call of the library API and pylosat reads `.ncbirc`, `<program>.ini` and the `NCBI_*` variables as the
  command line does at its start, and refuses what the command line refuses (such as
  `NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT=maybe`). Before, a search of sequences did so only with `-logfile`,
  `-conffile` or `-dryrun`, or where it needed NCBI's settings. Two such searches that worked now fail as the command
  line does: one with a `~/.ncbirc` entry that LOSAT refuses (such as `[BLAST] LONG_SEQID`), and one in a process
  that has no `HOME` and whose user ID has no entry in the user database (as in some containers).
- The documentation of the library's `Error` and of pylosat's errors says how each message is made from the command
  line's text.
