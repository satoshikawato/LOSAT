- `blastp`, `blastn`, `blastx`, `tblastn` and `tblastx` take NCBI's `-seqidlist` and `-negative_seqidlist` of a `-db`
  search, as BLAST+ 2.17.0 reads them: a text list (white-space separated ids, a leading `>` skipped, `#` comments) or
  a binary list of `blastdb_aliastool -seqid_file_out`, found on the database search path; on a version 5 database the
  ids are rewritten as NCBI rewrites them (`gi|N` dropped, upper-cased accessions with their versions when given) and
  looked up in the accession table (the first OID of an accession, a versionless id for every version; a negative list
  excludes a sequence only when every id the database stores for it is covered); on a version 4 database through the
  string ISAM. The search and the reports see only the kept sequences, as the taxid filters do: their count and length
  in the statistics, compact ordinals, and only the definition lines whose Seq-id the list names (or, for a negative
  list, does not name). NCBI's warnings are printed as it prints them (a text list on a version 5 database, a list
  whose ids are all gis, which leaves the database whole, a binary list recorded for another database); a binary list
  on a version 4 database, an empty list file, a database without the accession table and a list file off the search
  path stop with NCBI's errors. Two NCBI results are not reproduced: an accession of more than about 500 sequences
  (NCBI fails its lookup; LOSAT takes the first), and a sequence that a positive list keeps while matching none of its
  definition lines (NCBI leaves its hits out of the report without notice; LOSAT refuses the search). These options
  were refused before.
