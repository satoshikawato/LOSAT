- CI: the nightly Wasm threading workflow (`wasm-threading.yml`) now runs
  `LOSAT/tests/test_benchmark_wasi_reuse.js`, the test of the Wasm benchmark runner's evidence rules
  (a run whose bytes do not match is `PARITY_FAIL`, a bad command is `FAIL`, a stale output file is
  refused), against the threaded artifact and, with serial compatibility, the serial one; no workflow
  ran it before. The unrun-tests check is now `check_tests_in_ci.py`
  and reads `test_*.js` files as well as `test_*.py`.
