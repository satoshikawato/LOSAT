- blastn's query-batch slots (`LOSAT_X_QBATCHPAR`) also run in the threaded WebAssembly command: the loop's thread starts
  the threads of every slot's pool, as a WASI host starts threads only for the main thread. blastn and dc-megablast of
  31,168 D. melanogaster cDNAs x 3 Mb at 16 threads on the workstation (Node 26.8.2): 9.75 -> 7.88 s and 2.89 -> 2.67 s,
  max RSS 448 -> 922 MiB and 494 -> 789 MiB (native's slots: 106 -> 404 MiB and 136 -> 238 MiB). LOSAT Web starts no
  slots. Reports do not change.
