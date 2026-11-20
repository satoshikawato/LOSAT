- BLASTN `-subject` searches whose query batches are searched ahead in slots (`LOSAT_X_QBATCHPAR`, 2+ threads with a
  query file NCBI splits) no longer hold the whole query file: the query store of `LOSAT_X_BNQDROP` / `LOSAT_X_BNQPACK`
  now serves the slots too, with packed letters (`LOSAT_X_BNQSLOTS`). dc-megablast of 31,168 cDNAs x 3 Mb at 32
  threads: peak anonymous memory 341 -> 286 MiB (NCBI 313 MiB max RSS at one thread), time x0.986.
