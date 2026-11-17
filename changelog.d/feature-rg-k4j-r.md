- Tabular reports (outfmt 6, 7, 10 and 20) take fewer instructions and write the same bytes: a row compares its letters
  with the scoring matrix only when the report shows `positive` or `ppos`, and BLASTX writes its default outfmt 6 and 7
  with the same writer as any field list. At one thread, blastx `-db` outfmt 6 and 7 run about 0.5 % fewer instructions,
  a blastx field list with `qseq` and `sseq` 0.9 %, blastx `-subject` outfmt 6 0.3 % and blastp `-outfmt "6 std qseq sseq"`
  0.5 %. New fixture rows cover BLASTX query and `-subject` titles of every byte class, also bytes that are not UTF-8.
