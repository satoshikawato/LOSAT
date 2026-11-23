- `-export_search_strategy` to a FIFO or a pipe: when LOSAT stops at an internal error of its own (a panic), the file now gets
  every request made before it, as a regular file already did; it lost the requests that LOSAT holds for such a file until the
  run reaches their batches.
