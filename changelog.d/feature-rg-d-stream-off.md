- `blastn`: LOSAT's refusals that a later query batch meets now come after the reports of the query
  batches before it on every output path (as LOSAT's other refusals of a query batch do): a query of
  2^30 letters or more with `-ungapped` or matrix-only scoring (`-reward 0 -penalty 0`; before, nothing
  was written before it), a query batch whose concatenated query and `-window_size` pass 2^30 letters,
  and a query batch with an invalid query under scoring options without Karlin-Altschul values (before,
  the reports of the batches before were written only with some output formats). In the first batch
  these refusals come after the outfmt 0 header and the batch's query messages.
