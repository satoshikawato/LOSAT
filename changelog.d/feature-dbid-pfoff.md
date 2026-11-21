- A `-db` list of version 5 databases made with and without `makeblastdb -parse_seqids`, where `PRE_FETCH_SEQS_LIMIT`
  turns a query batch's prefetch off and a query's report reads a parsed Seq-id that a database of the list cannot
  look up: the JSON and XML2 formats (`-outfmt 13` to `16`) are now written as NCBI BLAST+ writes them (apart from the version line, which names LOSAT). `-outfmt 15`
  and `16` write `Error: [<program>] DB contains no accession info.` and NCBI BLAST+'s exception report, leave that
  query's report out and go on with the next query (exit status 0); `-outfmt 13` and `14` leave that query's file
  empty and stop with `LOSAT database error: DB contains no accession info.` (exit status 2). The messages of the
  taxonomy database that the sequences shown before that hit write come first, and a message that NCBI BLAST+ writes
  once per run (opening the taxonomy database, a volume's first lookup of a tax id) comes with the next sequence that
  NCBI BLAST+ reads (a later query's hit, or the prefetch of a later query batch, before that batch's error), also
  when the first such read comes after the hit that fails; LOSAT refused such a run. The other formats keep the
  refusal.
