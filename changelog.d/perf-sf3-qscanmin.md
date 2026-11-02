- `blastn` searches that read their query file twice (`LOSAT_X_QSTREAM`) now read the query lengths first in parts at
  two threads or more on query files of any size, not only on files of 4 MiB or more; the tuning variable
  `LOSAT_X_QSCANPARMIN` is gone. Reports do not change.
