- Treat 64-bit effective search space overflow (`SearchSpaceOverflow`) as a
  per-batch failure after earlier batches' reports in `tblastn`, `blastn`, and
  `tblastx` (matching NCBI BLAST+ partial output before the failing batch).
- In `blastn` split query batches where the batch has a valid context but every
  chunk has none, run the chunk OID-0 read check for empty databases and the
  batch-level traceback `BLAST_GapAlignSetUp` checks (count 0, average subject
  length reporting `Traceback failed`, then greedy set-up refusals) before
  traceback.
