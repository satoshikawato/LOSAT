- WASI command (`wasm32-wasip1`, serial and threaded): a `-query` or `-subject` named `.gz`, `.bz2` or `.zst`
  no longer stops the run with a panic ("no filesystem on wasm"). The inflated copy goes beside the input,
  as on the other platforms; WASI has no temporary directory to fall back to, so a directory that refuses
  the copy stops the run with LOSAT's message. Other platforms are unchanged.
