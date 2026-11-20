- LOSAT Web and the WASI builds: the threaded Wasm modules declare a shared memory of 65,535 pages
  (4 GiB less one 64 KiB page) instead of Rust's default 1 GiB (`-C link-arg=--max-memory=4294901760`
  in `LOSAT/.cargo/config.toml` and `web/adapter/.cargo/config.toml`; the code and data sections do not
  change). The app and the Node hosts step the maximum down (49,152, 32,768, 24,576, 16,384 pages) where
  a device cannot reserve it, and take pointers and lengths as unsigned values, so that pointers of
  2 GiB or more work in the ABI binding and the WASI imports (`docs/web/abi_v2.md` §3, §6). Brought to
  `dev` from #136 and #137 (Owner decision of 2026-10-12, plan DW-28).
