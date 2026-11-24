- The WebAssembly (WASI) `losat` command-line binaries are smaller: 0.59 MB (3.0 %) for
  `wasm32-wasip1` and 0.62 MB (2.5 %) for `wasm32-wasip1-threads`. They no longer carry LOSAT Web's
  programming interface (the `losat_web_*` functions) and the clap code behind it, which a command
  line never used: a WASI host now sees only the command's own entry. What the command writes is
  unchanged. LOSAT Web, LOSAT's WebAssembly library modules (`cargo build --lib`) and modules that
  others build from them keep that interface as before. To build the command line yourself as the
  release does, add `--features no-web-abi-v1` (see the installation chapter of the manual).
