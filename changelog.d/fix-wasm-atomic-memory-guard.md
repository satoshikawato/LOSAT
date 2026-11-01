- The threaded Wasm's shared-memory guard also covers atomic memory instructions and `memory.grow`. A rayon worker
  could stop a search with "memory access out of bounds" (Node exits 143 with empty output) when another thread grew
  the memory just before the worker used an address in the new pages. On ARM64 Linux the same happened on plain
  stores, in 11-29 of 100 blastp runs at 16 threads with the atomics guarded (V8's source treats Android the same way;
  no device was run). Each atomic instruction now refreshes a stale bound, and `memory.grow` records the new size so
  that the next atomic read of every thread sees it. The guard's record has version 3 and counts the atomics and
  grows; a module guarded by an older version must be rebuilt. No output change. Blastx at 1 thread runs about 4 % more
  instructions (1-2 % more cycles); the other modes are within 0.7 %.
