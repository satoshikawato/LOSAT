#!/usr/bin/env bash
# Known advisories in every dependency lock of the repository: the three Cargo lock files (RustSec,
# cargo deny with deny.toml), the npm lock file of web/ (every dependency, the development ones too)
# and the pinned Python test tools (python/tests/requirements.txt with what they pull in; pylosat
# itself has no Python dependency). No NCBI counterpart.
#
# Dependabot alerts watch only the default branch, main, which lags dev by weeks: on 2026-11-13 the
# minimatch advisory GHSA-g3hv-75wg-m249 was raised on main's old lock file only. This script checks
# the commit it runs on, for nightly.yml (job `audit`, dev every night) and
# release-readiness-v020.yml (job `security`, the release candidate). The advisory databases are
# fetched at run time, so it stays off pull requests.
#
# Needs cargo-deny, npm and pip-audit on PATH. Run from the repository root:
#   scripts/dependency_advisories.sh
set -uo pipefail

status=0
for manifest in LOSAT/Cargo.toml web/adapter/Cargo.toml python/Cargo.toml; do
  echo "::group::cargo deny check advisories ($manifest)"
  cargo deny --manifest-path "$manifest" --config deny.toml check advisories || status=1
  echo "::endgroup::"
done

echo "::group::npm audit (web/package-lock.json)"
(cd web && npm audit --package-lock-only --audit-level=low) || status=1
echo "::endgroup::"

echo "::group::pip-audit (python/tests/requirements.txt)"
pip-audit --strict --progress-spinner off --requirement python/tests/requirements.txt || status=1
echo "::endgroup::"

exit "$status"
