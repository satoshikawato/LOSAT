- `blastp`, `psiblast` and `blastx` write NCBI's tabular reports with a field list through the shared writer:
  `-outfmt "6|7|10|20 [delim=X] <fields>"` with all 50 fields of BLAST+ 2.17.0 (`staxids`, `sstrand`, `qcovs`,
  `qcovhsp`, `qcovus`, `sall*` and the others that were refused before), outfmt 10 and 20, and `delim=`. The
  default outfmt 6 and 7 keep their writers.
- The taxonomy database's names (`ssciname`, `scomname`, `sblastname`, `sskingdom` and their plural forms) are
  looked up with `taxdb.bti` on the search path in all five programs, with NCBI's messages on standard error;
  they were refused before.
