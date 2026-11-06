- Diagnostics: `LOSAT_X_REACH` takes only an absolute directory. A relative value such as `0` or `1` made a directory
  under the current one (12 stray files under `LOSAT/0/` and `LOSAT/1/`, removed by #715); the run now ignores it,
  writes nothing and says so once on stderr at its end. `.gitignore` covers `reach-*.tsv` and `LOSAT/[0-9]*/`. No report
  changes (one stderr line for a relative value).
- CI: a ratchet (`LOSAT/tests/check_lint_allow_ratchet.py`, run by `rust lint`) fails when a Rust file outside
  `LOSAT/tests/lint_allow_warnings.tsv` (lib.rs, main.rs and 23 test files at 2026-11-05) gets a blanket
  `allow(warnings)` or `allow(clippy::all)`, and when a listed file drops it without leaving the list (SWEBOK
  2026-11-04 audit, KA4).
