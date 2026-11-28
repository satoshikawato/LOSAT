- `blastp` and `rpsblast` now write the results of the query batches before the batch where they
  stop, as NCBI BLAST+ writes them, for two refusals: an alias file's `NSEQ 0` with `-searchsp` or
  `OLD_FSC`, after batches whose queries have no letter to search, and an effective search space
  that does not fit a 64-bit integer (a long `-dbsize` or alias `LENGTH`), after batches with hits.
  Before, both stopped before any output. Queries with no letter to search no longer meet the
  second refusal: their reports are written as in NCBI. An ungapped `blastp` search over a subject
  longer than 5,000,000 residues also writes those batches before its refusal.
  With `-outfmt` 0 to 4 the report's opening lines (the program, the references and the database)
  now come before these refusals also where the first batch meets them, as NCBI writes them.
- `psiblast -phi_pattern` no longer refuses a search over a long `-dbsize` or alias `LENGTH` where
  a query's `blastp` search space would not fit a 64-bit integer: a pattern search's space is the
  database length less the pattern's length per sequence, and LOSAT now reports it as NCBI BLAST+ does.
