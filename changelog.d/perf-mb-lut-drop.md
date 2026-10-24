- BLASTN drops a query batch's lookup tables when its preliminary stage ends, before its traceback, where
  NCBI's traceback stage drops its table (`LOSAT_X_BNLUTDROP`, in the strict set; `=0` holds them to the end of
  the batch as before). Peak anonymous memory of megablast of a 5.5 Mb genome falls from 135.7 to 131.3 MiB at
  one thread and from 142.6 to 133.7 MiB with `-db` at 16 threads; blastn of 31,168 cDNAs at 16 threads from 741
  to 706 MiB. Time and outputs do not change.
