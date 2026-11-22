- LOSAT Web writes the search strategy file of `-export_search_strategy NAME` as the command line does: the file, with
  the command line's bytes for the argv, comes among the run's files (as the `-logfile` file does) and can be saved from
  the run's results. Before, the Web GUI took the option, searched and wrote no file. The name is read as the command
  line on Linux opens it: a name it cannot open is NCBI's `File is not accessible` error, `/dev/null` gives no file, and
  `-` and the standard output are refused as on the command line. The standard error stream, the system's devices and
  files, and the file of `-query`, `-subject`, `-logfile` or `-window_masker_db` are refused in the Web GUI, where the
  command line writes the requests into that stream or file. The library and pylosat still refuse the option.
