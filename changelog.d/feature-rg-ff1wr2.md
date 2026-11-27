- LOSAT Web reads the name of `-export_search_strategy` against the run's files, as it reads `-logfile`: a name in a
  folder that the run's files do not have (`out/s.asn` without a folder `out`) or in one of their files, and the name
  of one of their folders, is NCBI's `File is not accessible` error, as on the command line run in that folder; before,
  the Web GUI wrote the file. A file of the run that no option names (such as `.ncbirc`) and a file that the search
  reads (such as the alias file of a `-db` database, even where it is missing) are refused in the Web GUI, as for
  `-logfile`; before, the Web GUI wrote a file of that name, where the command line's file takes the place of the run's.
