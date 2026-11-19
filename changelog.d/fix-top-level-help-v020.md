- `losat -help` names every program that `losat` runs: blastdbcmd, dustmasker, segmasker and
  convert2blastmask are listed under "Other programs" (derived from the entry names that
  `--install-shims` links). Stale notes are gone: makeblastdb says it writes format version 5 by
  default (4 with `-blastdb_version 4`), psiblast no longer says "its first round", blastx and
  tblastn no longer say "local" (both search `-db` too), and blastn names its four tasks.
- The `-outfmt` help of every search program says what it writes: blastn 0 to 18 and 20, blastp,
  blastx, tblastn and tblastx 0 to 16, 18 and 20, psiblast 0, 6, 7, 10 and 20, rpsblast 0, 6, 7 and
  10 (and what `-in_pssm` limits), that 6, 7, 10 and 20 take a field list, and the default fields.
  blastn no longer says custom fields are "not yet ported", tblastx no longer says "without custom
  fields", and blastp, blastx, tblastn, psiblast and rpsblast no longer show an empty line. One
  helper (`cli::outfmt_help`) writes all of them.
- LOSAT's own argument refusals and syntax errors (`-html`, `-remote`, an unknown option, a missing
  value) end with a line end, as clap's own errors do, so the shell prompt no longer lands on the
  message's line. NCBI's texts are not among them.
- No program's `-help` shows source references any more (blastn's `-perc_identity`, `-max_hsps` and
  `-subject_besthit` showed "NCBI reference: ncbi-blast/c++/... ```c ..."): the references are `//`
  comments next to the fields in every program's argument structs. The hint at the foot of blastn's
  help is true: `-h` ("see more with '-help'") prints the summary and `-help` ("see a summary with
  '-h'") the full text.
