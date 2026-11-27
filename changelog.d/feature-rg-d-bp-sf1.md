- `blastp` and `rpsblast` now write the results of the query batches before the batch where they
  stop, as NCBI BLAST+ writes them, for two refusals: an alias file's `NSEQ 0` with `-searchsp` or
  `OLD_FSC`, after batches whose queries have no letter to search, and an effective search space
  that does not fit a 64-bit integer (a long `-dbsize` or alias `LENGTH`), after batches with hits.
  Before, both stopped before any output. Queries with no letter to search no longer meet the
  second refusal: their reports are written as in NCBI. An ungapped `blastp` search over a subject
  longer than 5,000,000 residues also writes those batches before its refusal.
