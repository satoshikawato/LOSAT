- `blastn`'s query-batch slots (`LOSAT_X_QBATCHPAR`) are off by default again: at 16 threads they
  raised the max RSS far more than they cut the time (31,168 cDNAs x 3 Mb: blastn 107 -> 400 MiB for
  6.0 -> 4.6 s, dc-megablast 136 -> 238 MiB for 3.2 -> 1.8 s). `LOSAT_X_QBATCHPAR=1` turns them on.
  Reports do not change.
