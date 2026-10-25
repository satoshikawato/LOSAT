- Wasm (threaded): a thread no longer traps "memory access out of bounds" at its start in about 0.3 % of
  16-thread runs: the Node thread host and the browser thread worker call `memory.grow(0)` before
  `wasi_thread_start`, so a worker prepared before the engine grew the memory sees its new size.
