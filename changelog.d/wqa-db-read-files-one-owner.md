- `-logfile` and the files of a BLAST database: the command line and LOSAT Web now refuse the same
  files, the ones that reading the database of `-db` / `-filtering_db` checks or opens. The command
  line now also refuses the names that it looks up on the way (a missing alias `.nal` / `.pal`,
  `index.alx`, the files tried for a database name that does not resolve), which it used to create
  and then read as a database, and every volume of an `-index_name base,,first,last` range; it no
  longer refuses the MegaBLAST index files of a search without `-use_index true`, the copy of an
  `OIDLIST` file beside the first volume when the file the alias names exists, or GI-mask data files
  that the run does not read. A database's files that a run cannot open are no longer refused on the
  command line either: the GI-mask files, the GI index and the mask columns of a search without
  `-db_soft_mask` / `-db_hard_mask`, the id-lookup files of a search without an id list or taxonomy
  option, the RPS files of every search but `rpsblast` and `rpstblastn`; for `dustmasker`,
  `segmasker`, `windowmasker`, `makembindex` and a `makeblastdb -input_type blastdb` copy, which
  read the sequences in order, a volume's string index, id-lookup and metadata files, and the mask
  files except for `makeblastdb` and `makembindex -db_mask`, which read them. LOSAT Web now refuses
  the same files as the command line, from the same function, and a volume's component files.
