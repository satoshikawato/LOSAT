- The NCBI C++ Toolkit option `-help-full` runs in every program: it writes the help of `-help`, as NCBI BLAST+ 2.17.0
  does (LOSAT's programs hide no option that they run). It was refused in blastn, blastp, blastx, tblastn, tblastx,
  psiblast, rpsblast, rpstblastn, rpsbproc and makeblastdb.
- As in NCBI, the arguments after `-h`, `-help` or `-help-full` are not read: `-help -bogus` and `-help -out` print the
  help (they were argument errors in those programs but rpstblastn). makeblastdb refuses `-xmlhelp` by name, as the others do.
