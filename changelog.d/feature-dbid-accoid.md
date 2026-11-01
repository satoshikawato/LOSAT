- `-db` with a database built with `-parse_seqids`: every search program now looks the Seq-align id of a
  reported sequence up in the database as NCBI BLAST+'s formatters do (the LMDB keys of version 5, the ISAM
  files of version 4, under the search's id or taxonomy list) and refuses a sequence for which that lookup
  gives another sequence first, in the formats that read the subject's sequence (outfmt 0, 5, 6, 7, 10, 13
  to 17). NCBI BLAST+ 2.17.0 writes the names of the other sequence there, in BLASTP also the identities
  counted against its letters (`gb|AAA00001|` after `gb|AAA00001.1|`, `emb|X.1|` after `gb|X.1|`), or stops
  (`lcl|X` after `gb|X.1|`); LOSAT wrote the searched sequence's names. The same record in two volumes or in
  two databases of an alias is not refused (NCBI shows the copy, which is the same), and the Seq-annots and
  the archive (8, 9, 11, 12) are unchanged. A version 4 database without the ISAM files of `-parse_seqids`,
  in which NCBI finds no sequence, is now refused in those formats too.
