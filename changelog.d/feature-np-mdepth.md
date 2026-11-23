- makeblastdb: `-mask_data` reads mask files whose locations nest inside each other (`mix`,
  `equiv`) to any depth, where LOSAT stopped at 97 levels and more. It writes the database that NCBI
  BLAST+ writes with no stack limit; NCBI's makeblastdb crashes from about 6,500 to 6,900 levels
  with its default 8 MiB stack. A file damaged deep inside such a nest gets NCBI's report.
