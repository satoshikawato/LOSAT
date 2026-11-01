- The threaded Wasm's shared-memory guard also covers atomic memory instructions and `memory.grow`. A rayon worker
  could stop a search with "memory access out of bounds" (Node exits 143 with empty output) when another thread grew
  the memory just before the worker used an address in the new pages; on ARM64 Linux and Android the same happened on
  plain stores, in 11-29 of 100 blastp runs at 16 threads. Each atomic instruction now refreshes a stale bound, and
  `memory.grow` records the new size so that the next atomic read of every thread sees it. The guard's record has
  version 3 and counts the atomics and grows. No output change.
