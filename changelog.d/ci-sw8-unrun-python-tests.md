- CI: `rust lint` runs the tests of the Gate A certification tools (`compare_blastn_parity.py`,
  `certify_blastn_v010.py`, `audit_blastp_v010.py`, `audit_tblastx_v010.py`,
  `certify_integrated_runtime_v010.py`) and of `benchmark_wasm_threading.py`; five of these test
  files (35 tests) ran in no workflow, and the integrated-runtime one only in v0.1.0's dispatched
  release readiness. The new `check_python_tests_in_ci.py` fails `rust lint` when a `test_*.py` file
  of the repository is named in the `run:` commands of no workflow that runs on push, pull request or
  schedule (or of one they call); `python/tests/` runs as a directory in `python.yml`.
