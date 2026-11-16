- A `-db` list of version 5 databases made with and without `makeblastdb -parse_seqids`, where a query batch stops
  at a parsed Seq-id that a database of the list cannot look up and another reported sequence's load writes a
  message of the taxonomy database (a damaged `taxdb`): every search program now writes the messages of the
  sequences that NCBI BLAST+ loads before the stop (in the order of its Seq-id handles: gis, accessions, and
  `gnl|db|n` such as `gnl|BL_ORD_ID|n`, by their number, other ids last) and `BLAST Database error: DB contains no
  accession info.` (exit status 2), as NCBI BLAST+ does; LOSAT refused such a run. A load whose order NCBI takes from
  memory addresses (the same number, ids that are not packed, more than 200 ids in the batch) that writes a message
  stays refused. Newly refused, where LOSAT wrote the stop and NCBI BLAST+ goes on: a stopping Seq-id that another
  sequence of the database also holds, and one of the form `lcl|Query_N`.
