- `-logfile` and the files of a BLAST database: the command line and LOSAT Web now refuse the same
  files, the ones that reading the database of `-db` / `-filtering_db` checks or opens. The command
  line now also refuses the names that it looks up on the way (a missing alias `.nal` / `.pal`,
  `index.alx`, the files tried for a database name that does not resolve), which it used to create
  and then read as a database, and every volume of an `-index_name base,,first,last` range; it no
  longer refuses the MegaBLAST index files of a search without `-use_index true`, the copy of an
  `OIDLIST` file beside the first volume when the file the alias names exists, or GI-mask data files
  that the run does not read. LOSAT Web now also refuses a volume's component files and the GI-mask
  files of an alias's `MASKLIST`, as the command line does.
