- `-db` with a database built with `-parse_seqids`: a reported sequence whose Seq-align id is `bbs|`, `nat|` or
  `lcl|<integer>` in a version 5 database, or `gnl|BL_ORD_ID|n` of the sequence's own number beside another id
  (version 4 and 5), is now reported as NCBI BLAST+ 2.17.0 reports it. LOSAT refused these by the kind of id;
  NCBI BLAST+'s formatters look every Seq-align id up in the database, and find these sequences (version 5 by
  their LMDB keys). LOSAT now refuses only where that lookup finds no sequence or another one first: `gnl|PIG|`,
  `bbs|` and `nat|` of version 4, `lcl|300` after `bbs|300` (one LMDB key), `gnl|BL_ORD_ID|n` of another number.
  The refusal of such a `gnl|BL_ORD_ID|n` that names another sequence now says that the lookup gives that
  sequence first.
