- A `-conffile` name that is not UTF-8 is read as NCBI BLAST+ 2.17.0 reads it, in all 16 programs
  (the search programs, `makeblastdb`, `blastdbcmd`, `dustmasker`, `segmasker`, `convert2blastmask`,
  `makembindex`, `windowmasker` and `rpsbproc`). A name that is no file ends the program with NCBI's
  message, which holds the name's bytes (`Critical: (106.23) Registry file "conf\xff" not found,
  aborting.`), exit status 255; a name with a directory whose file cannot be opened gives NCBI's
  report with the bytes (exit status 2). Before, LOSAT refused such a name with its own error and
  exit status 1.
