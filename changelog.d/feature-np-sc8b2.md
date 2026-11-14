- `blastp` runs the first stage of a search whose report the query blocks do not write (outfmt 10, 18 and 20, `-query_loc`,
  a custom outfmt 6 or 7 that reads more of a `-db` sequence's titles, and others) in blocks of whole query batches too,
  each block with its own lookup table, and reports all queries together as before: the output does not change, and the
  lookup table of all queries is no longer held at once. 3,000 proteins against a six-proteome database with
  `-word_size 4 -outfmt 10` at 8 threads peak at 248 MB instead of 6.1 GB (NCBI 1.0 GB) in 15.2 s instead of 14.2 s
  (NCBI 40.0 s); with the default word size 3, 43.0 s and 48 MB instead of 58.5 s and 215 MB at one thread (NCBI 132.6 s
  and 77 MB). `LOSAT_X_BPQBLOCK=0` keeps one table for all queries ([`docs/agents/diagnostics.md`](docs/agents/diagnostics.md)).
