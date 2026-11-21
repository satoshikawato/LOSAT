- Alias files with `OID_MASK_TYPE` are read as NCBI BLAST+ 2.17.0 reads them (they were refused), such as the
  aliases `blastdb_aliastool -oid_masks 1` writes over a `makeblastdb -oid_masks 1` database: with the
  value's bit 0 set, the definition lines that hold a predicted RefSeq accession (`XP_`, `XM_`, `XR_`, ...)
  are hidden. A sequence left with no line is still searched, as in NCBI, but none of its hits is reported
  (after the `-max_target_seqs` cut, as NCBI's), and `blastdbcmd`, `segmasker` and the `windowmasker`
  masking step skip it where NCBI does; `dustmasker` stops at it with NCBI's error. The value is read as
  NCBI reads it (an invalid value gives NCBI's error), the last OID list above a volume decides, and
  `makembindex`'s new style index, as NCBI's, reads no alias key. A sequence with an id that NCBI warns
  about when it reads the id's accession is refused under the key, and so is an old style `makembindex`
  index of a sequence left with no line, where NCBI stops with a null-pointer report. A definition line with
  bytes NCBI reads as `#` is refused also where the key or a `-seqidlist`, `-negative_seqidlist` or taxid list
  removes it, as NCBI writes its error for such a line before it removes it. In an RPS database of several volumes, which NCBI
  searches without the alias files, a reported profile left with no line is refused, as NCBI's report does not
  find it (it writes `Unknown` or stops with a null-pointer report); with one volume its hits are dropped, as NCBI's.
  `makeblastdb -input_type blastdb` copies such an alias with the lines it keeps, as NCBI does, and refuses a sequence
  left with no line, where NCBI's `makeblastdb` stops with a null-pointer report.
