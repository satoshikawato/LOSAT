- `losat makeblastdb -mask_data` reads a mask string with a byte outside the visible ASCII range as makeblastdb does: the
  byte is stored as `#` and NCBI's error for it is written (`Error: (803.7) [makeblastdb]` with the member's path, the line
  or byte and the part of the string read), and a UTF8String tag for such a string gets NCBI's warning once per run. These
  files were refused before.
- `losat psiblast -in_pssm` and `losat tblastn -in_pssm` read a PSSM whose strings hold a byte outside the visible ASCII
  range as NCBI BLAST+ 2.17.0 does: the byte is stored as `#` and NCBI's error is written per byte
  (`Error: (803.7) [psiblast]` or `[tblastn]` with the member's path such as `PssmWithParameters.pssm.query.seq.descr.E.title`,
  the line and the part of the string read), after the formatting warnings and before the search's own messages. These
  PSSMs were refused before. The query's descriptors `name`, `comment`, `region`, `het`, `maploc`, `dbxref`,
  `create-date`, `update-date` and `user` are read as NCBI reads them (a byte outside the visible range is stored as `#` with
  NCBI's error; a `user` field's UTF8String takes any byte without one). Still refused: such a byte in a value that LOSAT's
  reader skips (the other descriptors, such as `org`, `pub`, `source` and `molinfo`; the query's `annot`; the `Seq-inst`
  members other than `repr`, `mol`, `length` and `seq-data`; the `constraints` and `annotatedSites` of the parameters).
- `losat rpsbproc` refuses an archive with a byte outside the visible ASCII range in any VisibleString (a Seq-id's, a Dbtag's,
  a date's as well as a title's), naming its line and the part of the string read.
