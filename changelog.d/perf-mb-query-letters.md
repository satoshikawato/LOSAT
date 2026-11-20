- BLASTN drops a query batch's letters once the batch has encoded them, where no report shows them (outfmt 6, 7, 8,
  9, 12 and SAM without `SQ SR`) and the records are read before the search (`LOSAT_X_BNQDROP`, in the strict set;
  `=0` keeps them to the end of the search as before). Peak anonymous memory of megablast of a 5.5 Mb genome falls
  from 131.6 to 126.4 MiB at one thread and from 135.5 to 129.1 MiB with `-db` at 16 threads; megablast of 31,168
  cDNAs at one thread from 219.7 to 215.3 MiB. Time and outputs do not change.
