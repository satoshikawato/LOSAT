- `losat blastdbcmd`, `dustmasker`, `segmasker`, `convert2blastmask`, `makembindex` and `windowmasker`
  read a command-line word that is not UTF-8 as NCBI BLAST+ 2.17.0 does, where they refused the whole
  command line before reading it. NCBI's errors of the other words (an unknown argument, a missing
  value, a value outside its list) and the help words (`-h`, `-help`, `-xmlhelp`) come first, and an error
  quotes the word's bytes as NCBI's does. `-logfile` with a name that is not UTF-8 runs and writes that
  file. Any other value that is not UTF-8 is still refused, now after the arguments are read.
