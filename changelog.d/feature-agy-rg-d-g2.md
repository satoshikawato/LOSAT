- Treat 64-bit effective search space overflow (`SearchSpaceOverflow`) as a
  per-batch failure after earlier batches' reports (and after the prolog on a
  single-batch refusal) in `tblastn`, `blastn`, and `tblastx`, and remove the
  now-redundant `Off::SearchSpace` report-stream guard in `blastn` and
  `tblastx`.
- In `blastn`, avoid writing report-time prologs (`outfmt 15`, `16`, and `20`)
  twice when a refused first batch (`-parse_deflines` duplicate query IDs) hits
  a non-`NativeError` set-up refusal (count 0, greedy set-up, or search space
  overflow).
- In `blastn` split query batches where the batch has a valid context but every
  chunk has none, run the chunk OID-0 read check for empty databases and the
  batch-level traceback `BLAST_GapAlignSetUp` checks (count 0, average subject
  length reporting `Traceback failed`, then greedy set-up refusals) before
  traceback.
