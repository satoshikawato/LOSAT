- LOSAT Web: when `-out_pssm` or `-out_ascii_pssm` names a file that cannot be opened (a folder that does not
  exist, a file used as a folder, a directory, a name that is too long), the error now gives the command line's
  reason, for example `No such file or directory (os error 2)` for `sub/x.asn`, where it said `cannot open file`.
  A PSSM file written to `/dev/full` now says `No space left on device (os error 28)`. The refusal, the exit code and
  the kept output are unchanged.
