- LOSAT Web: when `-out_pssm` names a file that cannot be opened (a folder that does not exist, a file used as a
  folder, a directory or a name ending in `/`, a name that is too long, a place the user cannot write), the error
  now gives the command line's reason, for example `No such file or directory (os error 2)` for `sub/x.asn`, where
  it said `cannot open file`. A PSSM file written to `/dev/full` now says `No space left on device (os error 28)`.
  The exit code and the kept output are unchanged.

- `losat psiblast -save_pssm_after_last_round` with an `-out_pssm` name that ends in `/`, `.` or `..` (`nosuch/`,
  `nodir/..`) now stops after the query's report, where NCBI BLAST+ stops, and names the system's reason for that
  name. Before, it wrote every report and then refused, and said `not a file name` or `Not a directory` where the
  system says `No such file or directory` or `Is a directory`.
