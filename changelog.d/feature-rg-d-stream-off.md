- `blastn`: with `-ungapped` or matrix-only scoring (`-reward 0 -penalty 0`), the refusal of a query of
  2^30 letters or more now comes after the reports of the query batches before it (in the first batch,
  after the outfmt 0 header), as LOSAT's other refusals of a query batch do; before, nothing was written
  before it.
