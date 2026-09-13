# Wasm performance optimization evidence

The current result is [run-20260913-02/REPORT.md](run-20260913-02/REPORT.md): a TBLASTX-only TurboFan profile for cold Node/WASI comparisons. [run-20260913-01/REPORT.md](run-20260913-01/REPORT.md) preserves the original decisions and known parity limitations.

## Git packaging

The commit includes implementation, required comparison/plotting prerequisites, reports, policy/decision/command/result metadata, timing and diagnostic logs, small plots, and the full evidence checksum inventories. Large raw outputs, V8 profiles, compiled artifacts, source/driver archives and frozen runner copies remain in the original local evidence directories; they are not deleted or rewritten. `commit-scope.json` lists the exact disposition and recorded checksums. A checkout of Git alone therefore does not contain the full raw evidence set, and each run's `SHA256SUMS` inventories that complete local set rather than only the Git subset. Evidence can be regenerated using the documented commands and fresh output paths, but that creates new measurements.

Source files are staged with the repository's LF endings; original working files and sealed evidence retain their validated bytes. This does not alter executable Python/Rust/JavaScript tokens. Existing unrelated changes, including line-ending-only differences elsewhere, remain outside the commit. The plotting helpers and Wasm-only wrapper are existing uncommitted prerequisites of the already validated comparison harness, included so a checkout has the tested dependency set.

The recorded benchmarks and 39-test gate were completed before this commit. Commit preparation confirms staged source bytes match the validated files after LF normalization. An unnecessary isolated test attempt stopped on fixture-containment checks; its failure log is retained. At the user’s direction, no further tests were run; commit-validation.json points to the existing successful gate. Native/Wasm search source is unchanged by this task. P4 remains PARTIAL for the documented megablast differences and incomplete release certification.
