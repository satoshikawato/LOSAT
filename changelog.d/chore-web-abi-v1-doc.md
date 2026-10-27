- Documentation only, no code change: new `docs/web/abi_v1.md` describes the v1 Web ABI
  (`losat_web_*`: return values, FASTA handles, result and error buffers, input checks,
  memory ownership, request format). It and `docs/web/abi_v2.md` have an "After a trap"
  section: with `panic = "abort"` a panic is a Wasm trap, so the host discards the
  instance and instantiates the module again.
