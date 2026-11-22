- `losat blastdbcmd`, `dustmasker`, `segmasker`, `convert2blastmask`, `makembindex` and `windowmasker`
  read a command-line word that is not UTF-8 as NCBI BLAST+ 2.17.0 does. Before, they refused the whole
  command line before reading it. Now NCBI's errors in the other words (an unknown argument, a missing
  value, a value outside its list) and the help words (`-h`, `-help`, `-xmlhelp`) come first, and an
  error quotes the word's bytes as NCBI's does. `-out` and `-logfile` with a name that is not UTF-8
  write that file, and `-dryrun` ends with exit status 0. Any other value that is not UTF-8 is still
  refused, now where the program would start its work.
