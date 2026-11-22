- `losat blastdbcmd` reads a hexadecimal `-line_length` or `-ipg` as NCBI BLAST+ 2.17.0's blastdbcmd
  does, where it refused one before: `-line_length 0x5` writes lines of 5 letters and `-ipg 0x3E9`
  retrieves PIG 1001; `0x`, `0x0`, `0x1p3` and a value beyond 32 bits get NCBI's error message
  and exit status.
