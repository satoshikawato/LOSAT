- BLASTN drops a query batch's letters once the batch has encoded them, where no report shows them (outfmt 6, 7, 8,
  9, 12 and SAM without `SQ SR`) and the records are read before the search (`LOSAT_X_BNQDROP`, in the strict set;
  `=0` keeps them to the end of the search as before). Outputs do not change.
