- BLASTN packs the letters of the query records it reads before the search as it reads them, where a batch drops
  its letters after encoding them (`LOSAT_X_BNQDROP`) and the query file holds two batches or more
  (`LOSAT_X_BNQPACK`, in the strict set; `=0` keeps the letters as before), and makes a batch's letters again when
  the batch starts. Peak anonymous memory of megablast of 31,168 cDNAs falls from 215.3 to 156.7 MiB at one thread
  (max RSS 221.9 to 163.5 MiB) and from 234.4 to 175.3 MiB with `-db` at 16 threads. Time and outputs do not change.
