- `blastp`, `blastn`, `blastx`, `tblastn` and `tblastx` take NCBI's taxid filters of a `-db` search: `-taxids`,
  `-taxidlist`, `-negative_taxids`, `-negative_taxidlist` and `-no_taxid_expansion`, as BLAST+ 2.17.0 reads them
  (the taxid lookup files of a version 5 database, `taxonomy4blast.sqlite3` on the database search path for the
  descendants, the list files with any line end). The search and the reports see only the kept sequences: their count
  and length in the statistics, compact ordinals (`gnl|BL_ORD_ID|N`), and only the definition lines of kept taxids.
  A version 4 database, a positive taxid that no sequence has, and a list file that cannot be read stop with NCBI's
  errors; without `taxonomy4blast.sqlite3` NCBI's warning is printed and the database is searched whole. `-taxids 1`
  (whose root is its own parent, where NCBI's expansion does not end) keeps every sequence of the taxonomy. A negative
  list that excludes every sequence is refused (NCBI searches an empty database). These options were refused before;
  `-seqidlist` and `-negative_seqidlist` still are.
