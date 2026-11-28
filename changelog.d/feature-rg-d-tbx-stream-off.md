- `tblastx` with `-db`: a run on a database with parsed Seq-ids of more than 2^20 sequences (such as
  NCBI's `nt` and `core_nt`), or one with a sequence whose report might be refused, now writes each
  query batch's report as soon as the batch is searched, as other runs do; before, it held every
  report until the end of the run. Runs on smaller databases with parsed Seq-ids start the search
  without first decoding every definition line (2.3 to 5.6 µs per sequence at one thread, up to about
  6 s for 2^20 sequences).
- `tblastx`: a refusal that a later query batch meets now comes after the reports of the query
  batches before it, as NCBI BLAST+ writes its errors of a batch; before, nothing was written before
  it. This covers a query batch whose concatenated query and `-window_size` pass 2^30 letters and a
  database whose sequences cannot be read (in every output format, after the batch's query messages),
  and a reported sequence whose names, Seq-ids, taxonomy record or taxonomy-cache slot LOSAT does not
  write as NCBI BLAST+ does, in a `-db` or `-subject` set (with `-outfmt 0` to `4`, `6` and `7`).
- `tblastx`: in the first query batch too, including a run of one batch, the window and read
  refusals now come after what NCBI BLAST+ writes before its error there: the `-outfmt 0` header, the
  XML2 and JSON prologs (`-outfmt 13` to `16`), the `-outfmt 20` header line and the batch's query
  messages. Before, LOSAT wrote nothing before them.
