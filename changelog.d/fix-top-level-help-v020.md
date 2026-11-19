- `losat -help` names every program that `losat` runs: blastdbcmd, dustmasker, segmasker and
  convert2blastmask are listed under "Other programs" (derived from the entry names that
  `--install-shims` links). Stale notes are gone: makeblastdb says it writes format version 5 by
  default (4 with `-blastdb_version 4`), psiblast no longer says "its first round", blastx and
  tblastn no longer say "local" (both search `-db` too), and blastn names its four tasks.
