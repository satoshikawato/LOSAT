- `tblastx`: LOSAT's refusals that a later query batch meets now come after the reports of the query
  batches before it (before, nothing was written before them): a query batch whose concatenated query
  and `-window_size` pass 2^30 letters and a database whose sequences cannot be read (in every output
  format, after the batch's query messages; in the first batch after the outfmt 0 header), and a
  reported sequence whose names, Seq-ids, taxonomy record or taxonomy-cache slot LOSAT does not write
  as NCBI BLAST+ does, in a `-db` or `-subject` set (with `-outfmt 0` to `4`, `6` and `7`). `-db` runs
  no longer read every definition line before the search.
