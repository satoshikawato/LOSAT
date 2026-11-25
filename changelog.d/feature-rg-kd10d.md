- `losat blastdb_aliastool` now reads every option value as its bytes, as NCBI BLAST+ 2.17.0 does: a
  database, list, input or output name with a byte that is not UTF-8 (a Latin-1 name, say) opens, is
  written into the alias file, the binary list and the messages as those bytes, and a missing one gives
  NCBI's message with that byte; a title and the lines of `-dblist_file` keep their bytes too. Before,
  it refused such a value and such a `-dblist_file` line.
