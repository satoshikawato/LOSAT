- LOSAT Web writes the search strategy file of `-export_search_strategy NAME` as the command line does: the file, with
  the command line's bytes for the argv, comes among the run's files (as the `-logfile` file does) and can be saved from
  the run's results. Before, the Web GUI took the option, searched and wrote no file. The name is read as the command
  line on Linux opens it: a name it cannot open is NCBI's `File is not accessible` error, `/dev/null` gives no file, and
  `-` and the standard output, where the command line writes the requests, are refused. The standard error stream, the
  system's devices and files, the file of `-query`, `-subject`, `-logfile` or `-window_masker_db`, and a file that the
  search reads from the run's folder (a file of the database of `-db` or `-filtering_db`, the list of a filter) are
  refused in the Web GUI, where the command line writes the requests into that stream or file. A search of `-db` writes
  its file too, with the database in its requests. The library and pylosat still refuse the option. With the option, a
  `-outfmt` value that is not a format is reported where the command line reports it, after the errors of the other
  arguments (such as `-query_loc 10-5`), in the Web GUI as on the command line; before, the Web GUI reported the
  `-outfmt` value first.
- `-html` with a malformed `./.ncbirc` and `-export_search_strategy -`: BLASTX, BLASTN and TBLASTX write the start of
  the HTML report and exit with status 255, as NCBI BLAST+ and the run without the option do; before, they wrote
  nothing to the standard output.
