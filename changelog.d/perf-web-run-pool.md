- LOSAT Web: a run's searches (one per hit list size of its formats with `-num_descriptions` / `-num_alignments`, and
  BLASTP's block search and its fallback) share one thread pool, and the browser holds one thread worker per thread
  instead of two (a spare set for a second pool). At 16 threads the page holds 15 fewer workers (~2 MiB each).
