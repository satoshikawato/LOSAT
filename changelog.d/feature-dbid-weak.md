- `-db` with a database built with `-parse_seqids`: LOSAT now refuses a reported sequence whose Seq-align id NCBI
  BLAST+'s scope matches weakly to the id of a query read with `-parse_deflines` (`sp|P00001.1|NAME` and a query
  `sp||NAME` or the other way, `gb|X|` and a query `gb|X.1|`, `pir|ABCD|` and a query `pir||ABCD`). NCBI BLAST+
  2.17.0 writes the query's names for such a sequence (BLASTX can stop instead); LOSAT wrote its own, with exit
  status 0.
- `-db` with a database built with `-parse_seqids` and the environment variable `PRE_FETCH_SEQS_LIMIT` set to an
  integer other than 2147483647: LOSAT now refuses a reported sequence whose Seq-align id an id of another
  sequence of the database matches weakly (`gb|X|` and `gb|X.1|`, `sp|P00001.1|NAME` and `sp||name`), in a version 5
  database any whose Seq-align id has a name or is a `pir` or `prf` id, in outfmt 0, 5, 6, 7, 10, 13 to 17 and 20. With 0 or a value below the number of
  Seq-aligns of a batch, NCBI BLAST+ reads the reported sequences one by one and can write the other sequence's
  names for it; LOSAT wrote its own.
