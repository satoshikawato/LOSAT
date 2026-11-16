- The NCBI C++ Toolkit option `-dryrun` runs as in NCBI BLAST+ 2.17.0 in blastn, blastp, blastx, tblastn, tblastx,
  psiblast, rpsblast, makeblastdb and blastdbcmd: the arguments are read and checked, then nothing runs, exit status 0.
- blastn: `-window_masker_taxid` with `-no_taxid_expansion` is an argument error, as in NCBI (it ran into a BLAST engine
  error). makeblastdb: `-taxid` reads NCBI's integer forms (`0x10` is 16, as in NCBI; it was refused).
