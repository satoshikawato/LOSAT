- `-seqidlist` with a text list whose ids are all gi ids (`gi|N`), on a version 5 database: NCBI BLAST+ drops the gi ids
  and warns `Empty seqid list`. It then searches the whole database only when the `-db` names reach one volume and at
  least one path of alias files to it crosses no `FIRST_OID`-`LAST_OID` range that cuts the volume; otherwise (two or
  more volumes, such as a multi-volume database, several databases in `-db` or an alias file over several, or a volume
  that every path reaches through such a range) it searches no sequence, with each query's `Filtering resulted in an
  empty database.`. LOSAT searched the whole database in every case. It now does what NCBI BLAST+ does, in every search
  program. A negative list of gi ids only still excludes nothing. LOSAT also writes NCBI's warning `Ignoring extra parts
  (synonyms?) in FASTA-style ID` for a list id of several parts (such as `gi|1|ref|NP_1.1|`) on a version 5 database.
