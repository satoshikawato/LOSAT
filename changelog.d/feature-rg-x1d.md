- Gapped extensions with an X-drop past what the 16-bit vector lanes hold (raw X-drop above about
  26,700, for example `-xdrop_gap_final 1e5` in blastn) run on 32-bit vector lanes on x86_64, aarch64 and
  wasm32 instead of the scalar loops, and an X-drop larger than any the extension can use is lowered to that
  bound first, with the same results: `blastn -task blastn -xdrop_gap_final 1e5` of 12 kb against 20 kb at one
  thread took 22.8 s of CPU instead of 72.5 s (NCBI BLAST+ 2.17.0: 82.3 s), and 265 MB instead of 302 MB.
  `LOSAT_X_DP32=0` runs the previous path.
