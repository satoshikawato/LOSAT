- LOSAT Web's page uses `-conffile`, `-logfile` and `-dryrun` typed in Additional options: it asks for the file that
  `-conffile` names and the run reads it under that name; the `-logfile` file is shown after the output formats, saved
  with Download All, and kept under a failed run's status (TBLASTX writes its run messages and NCBI's error there); a
  dry run says that it searched nothing. Before, every `-conffile` run failed with NCBI's `Registry file "NAME" not
  found`, the `-logfile` file was dropped, and a dry run looked like a search with no hits.
