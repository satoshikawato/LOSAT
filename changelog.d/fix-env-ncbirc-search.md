- The configuration files `.ncbirc` and `<program>.ini` are looked for where NCBI BLAST+ 2.17.0 looks: an
  empty `NCBI_CONFIG_PATH` means no directory (LOSAT searched the default path), the home directory is
  `APPDATA` before `USERPROFILE` on Windows and the account's when `HOME` is not set on Unix, and the
  program's directories are those NCBI has when it makes the path, before the program's arguments are set:
  the executable's directory on Linux, none on macOS and Windows (LOSAT added `losat`'s directory there).
  With `BLAST_USAGE_REPORT` set to a false value, or `.ncbirc` turned off in the environment
  (`NCBI_DONT_USE_NCBIRC`, `NCBI_CONFIG__NCBI__DONT_USE_NCBIRC`), NCBI makes the path later, and the directory
  of the name the program was run by comes first, so a `blastn` link's directory is searched. Without `HOME` and without an
  account entry for the user, NCBI looks the account up by the login name, which LOSAT cannot do: it stops.
- tblastx (before it writes anything), and blastx with `-ungapped -comp_based_stats 0` (at its first query batch,
  in the frame of the batch's error), stop with a message when NCBI BLAST+ would read a
  file for a `-matrix` that it does not have built in (`NCBI_DATA_PATH`, `[NCBI] Data`, `BLASTMAT`,
  `BLASTMAT/nt` or `./data`): NCBI BLAST+ 2.17.0 aborts after a matrix it reads, and fails with another code on
  a file that is not a matrix. Without such a file, or with one it cannot open, the output is NCBI's error, as
  before.
