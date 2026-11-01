- The threaded Wasm's shared-memory guard also covers atomic memory instructions. A rayon worker could stop a
  search with "memory access out of bounds" (Node exits 143 with empty output) when another thread grew the memory
  just before the worker used an address in the new pages. Each atomic instruction now refreshes a stale bound first,
  as `memory.fill` and `memory.copy` do. The guard's record has version 2 and counts the atomics. No output change.
