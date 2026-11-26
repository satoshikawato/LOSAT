- BLASTN: for `-db` searches, an alignment that NCBI BLAST+ keeps outside its query context whose
  Seq-align query start lies below 0 in record coordinates (for example with `-xdrop_ungap 1e10`,
  on either strand, gapped or ungapped, with or without `-query_loc`) is no longer refused before
  reporting: LOSAT now reproduces NCBI BLAST+ 2.17.0's exact per-report behavior across all output
  formats (outfmt 0 skips the alignment with `alignment skipped` or stops with `Row is empty`
  before a gapped description table or query-range sort; outfmt 1 to 4 stop in the anchored merge
  or before a gapped description table; outfmt 5, 15 and 16 post NCBI's error and finish while 13
  and 14 stop; outfmt 6, 7, 10 and 20 write rows for fields that do not read the query row and stop
  at the alignment or at gapped `qcovs`/`qcovus` otherwise; outfmt 8, 9, 11, 12 and 18 write the
  report with NCBI's wrapped unsigned `TSeqPos` starts; outfmt 17 writes SAM records or stops on
  query-reference sort after flushing earlier queries). Batch prefetch still runs before report
  formatting, and `-subject` searches keep their batch refusal after earlier batch output.
