- BLASTX's tabular reports (outfmt 6, 7, 10 and 20, the default fields and field lists alike) and the field lists of
  blastp, psiblast, tblastn and tblastx that do not show `positive` or `ppos` take slightly fewer instructions, with the
  same output: at one thread about 0.2 to 0.9 % fewer for the blastx and blastp searches measured, depending on the
  search. The default outfmt 6 and 7 of blastp, blastn, tblastn and tblastx are not affected.
