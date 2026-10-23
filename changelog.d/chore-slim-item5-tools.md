- Removed developer tools that no workflow, script, skill or document runs: five TBLASTX hit-gap analysis scripts
  under LOSAT/tests/ (one of them also as a stray copy under LOSAT/LOSAT/tests/), tests/analyze_blast_results.py,
  tests/analyze_debug_output.py, and the seven BLASTX stage probes in LOSAT/examples/ (built by every
  `cargo test` and `cargo clippy --all-targets`). Outputs do not change.
