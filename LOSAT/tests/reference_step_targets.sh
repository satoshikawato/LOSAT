#!/usr/bin/env bash
# The cargo test targets of the reference step (`LOSAT_X_DEFAULT=0`) of ci.yml and
# platform-checks.yml, as arguments for `cargo test` or `cargo nextest run` (CI): the library and
# binary unit tests and every integration test file without the marker line below.
#
# A file carries the marker when each LOSAT run it makes removes the caller's `LOSAT_*` variables
# (all but `LOSAT_BLASTX_WORKER_LOG`, which selects no code path) and sets its own, and it calls
# no library code in-process: its runs are then the same in both steps, and the default step
# already ran them. Other launches in such a file may keep the caller's environment only when the
# command reads no switch (`makedb`: `src/writedb` reads none). A file without the marker runs in
# both steps.
#
# The step runs no doc tests: the crate has none (2026-10-22); the default step runs them.
#
# Usage (from LOSAT/): targets=$(bash tests/reference_step_targets.sh)
#                      cargo test --all-features $targets
#                      (CI: cargo nextest run --all-features $targets)
set -euo pipefail
cd "$(dirname "$0")"

marker="// CI reference step: skipped"
args=(--lib --bins)
skipped=0
for file in *.rs; do
    if grep -qF "$marker" "$file"; then
        if ! grep -qF 'starts_with("LOSAT_")' "$file"; then
            echo "$file: marked, but it does not remove the caller's LOSAT_* variables" >&2
            exit 1
        fi
        if grep -Eq '(^|[^A-Za-z_])LOSAT::' "$file"; then
            echo "$file: marked, but it calls the library in-process" >&2
            exit 1
        fi
        skipped=$((skipped + 1))
        continue
    fi
    args+=(--test "${file%.rs}")
done
echo "reference step: $((${#args[@]} / 2 - 1)) integration test files, $skipped skipped (marked)" >&2
echo "${args[@]}"
