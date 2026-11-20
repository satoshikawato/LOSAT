- Web: a run with a structured `-outfmt` (5, 8, 9, 11 to 17) also shows the queries' warnings, such
  as "Sequence contains no data", that the command line writes for outfmt 0, 6 and 7, which the run
  also writes. They were missing; the command line's own output and messages do not change.
