- BLASTX's tabular reports (outfmt 6, 7, 10 and 20, the default fields and field lists alike) and the field lists of
  blastp, psiblast, tblastn and tblastx that do not show `positive` or `ppos` take slightly fewer instructions, with the
  same output: at one thread about 0.2 to 0.9 % fewer for the blastx, blastp, tblastn and tblastx searches measured,
  depending on the search. The default outfmt 6 and 7 of blastp, blastn, tblastn and tblastx are not affected.
- The tabular reports of blastx and rpstblastn, and the field lists of blastp, psiblast, rpsblast, tblastn and tblastx,
  now write each row as soon as it is made instead of holding all of a query's rows first, so a query with many hits no
  longer adds their text to the memory in use. A library caller's writer receives one write per row, as it already did
  from the default outfmt 6 and 7 of blastp; the writers the library is given are not buffered by it.
