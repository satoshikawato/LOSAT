- Web: a run with a structured `-outfmt` (5, 8, 9, 11 to 17) also shows the queries' warnings, such
  as "Sequence contains no data", that the command line writes for outfmt 0, 6 and 7, which the run
  also writes, each once. A BLASTX run whose tabular results stop after a query's warning
  while the pairwise results complete (a run with `-num_descriptions`) shows that warning once,
  not twice. The command line's own output and messages do not change.
