- `-seqidlist` with a text list whose ids are all gi ids (`gi|N`), on a version 5 database of two or more volumes (a
  multi-volume database, several databases in `-db`, an alias file over several), or on one volume inside an alias file's
  `FIRST_OID`-`LAST_OID` range: NCBI BLAST+ drops the gi ids, warns `Empty seqid list`, and searches no sequence (each
  query's `Filtering resulted in an empty database.`); LOSAT searched the whole database. LOSAT now does the same as NCBI
  BLAST+, in every search program. On one volume that the alias files leave whole, both still search the whole database.
