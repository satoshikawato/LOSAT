- The variables `LOSAT_X_BNTABPAR` and `LOSAT_X_TBXRPTPAR` are gone: `LOSAT_X_PAIRWISEPAR=0` keeps the whole report
  stage on one thread (the outfmt 0 text of `blastn`, `blastp` and `tblastx`, the outfmt 6 / 7 rows of `blastn` and
  `tblastx`, and the displayed identities, engine counts and outfmt 0 rows of `tblastx`). The outfmt 6 / 7 rows of
  `tblastx` are now formatted on the search's threads for every query, in the pieces and windows `blastn` uses; before,
  a query of 1,024 rows or fewer was written on one thread. Reports do not change. <<NUMBERS>>
