- With a taxonomy database (`taxdb.bti`) on the search path, `-db` searches of every program write NCBI's
  messages about it on stderr in outfmt 0, 5, 6 and 7, and SAM for blastn, as NCBI BLAST+ 2.17.0 does (they
  were written only in 13 to 16): the texts of opening it (for example `Tax database file has wrong magic
  number.`) and of a volume's lookup of a taxid (`Offset error at end of taxdb file.`), after the reader's
  warnings of each query batch and before its first query's own warnings. Without `taxdb.bti` nothing changes.
