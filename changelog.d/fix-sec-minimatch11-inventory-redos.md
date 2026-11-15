- Security: the Web lint tooling uses minimatch 11.0.5 (npm `overrides` in `web/package.json`; dev-only, fixes
  GHSA-g3hv-75wg-m249) and `docs/evidence/losat_web_e2g/inventory_refs.py` matches function names with a regex that no
  longer backtracks exponentially (CodeQL `py/redos`; same matches). No change to LOSAT's output.
