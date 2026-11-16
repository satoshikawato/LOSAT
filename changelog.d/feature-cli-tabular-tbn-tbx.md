- `tblastn` and `tblastx` write NCBI's tabular reports with a field list: `-outfmt "6|7|10|20 [delim=X] <fields>"`
  with all 50 fields of BLAST+ 2.17.0 (`qseqid` ... `qcovus`, `std`), outfmt 10 and 20, and `delim=` (NCBI's
  delimiter errors included), on `-subject` and `-db`. Field lists and `delim=` were refused before. The default
  outfmt 6 and 7 keep their writers.
