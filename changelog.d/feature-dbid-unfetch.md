- `-db` with a database built with `-parse_seqids`: a reported sequence whose Seq-align id is `bbs|`, `nat|` or
  `lcl|<integer>` in a version 5 database, or `gnl|BL_ORD_ID|n` of the sequence's own number beside another id
  (version 4 and 5), is now reported as NCBI BLAST+ 2.17.0 reports it, unless another sequence holds the same id.
  LOSAT refused these by the kind of id; NCBI BLAST+'s formatters look every Seq-align id up in the database and
  find these sequences (version 5 by their LMDB keys). Where that lookup finds no sequence or another one first
  (`gnl|PIG|`, `bbs|` and `nat|` of version 4, `lcl|300` after `bbs|300`, `gnl|BL_ORD_ID|n` of another number),
  LOSAT refuses as before; the refusal of a `gnl|BL_ORD_ID|n` that names another sequence now says that the
  lookup gives that sequence first.
- `-db` with a database built with `-parse_seqids`: LOSAT now refuses a reported sequence whose Seq-align id
  another sequence of the database also holds, in any letter case (`lcl|555`, and `gb|X.1|` with `lcl|555` in its
  second definition line; `lcl|abc` and `lcl|ABC`; one id in two databases of an alias), unless the two are the
  same record. Where both are reported, NCBI BLAST+'s scope can write the other sequence's names and identities
  for it; LOSAT wrote its own, with exit status 0.
