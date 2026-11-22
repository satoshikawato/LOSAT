- tblastx against a single subject sequence (or with other subjects of fewer than 3 bases in all), when one query
  batch holds so many very short queries that its average translated length is 0 (for example a 1,200-nt query with
  1,300 queries of 3 nt, or 300 with `-strand plus` or `minus`), keeps writing linked hits where NCBI BLAST+ 2.17.0
  crashes; this is now a documented difference (the manual's "Differences from NCBI BLAST+" page). When the other
  subjects add at least one codon, the results are NCBI's.
