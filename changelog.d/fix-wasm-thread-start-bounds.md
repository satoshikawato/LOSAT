- Wasm (threaded): the Node thread host and the browser thread worker call `memory.grow(0)` before
  `wasi_thread_start`, so a worker prepared before the engine grew the memory for its thread sees the
  new size. In the Node host about 0.3 % of runs of one 16-thread check trapped "memory access out of
  bounds" at the thread's start without it (not measured in a browser).
