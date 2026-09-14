# Commit handoff — provisionally adopted

Title: **Fix megablast X-drop parity and reduce Wasm search data movement**

Match NCBI's automatic megablast initial X-drop, keep BLASTN HSP payloads in
fixed storage during endpoint pruning, write BLASTP traceback rows through
reserved storage, and sort TBLASTX final linking results through references.
Align native x86_64 loops using the separately qualified Cargo configuration.
Use threaded Wasm for standard builds and comparisons, with serial available
as an explicit compatibility option.

The user requested provisional adoption and no further measurements. The report
records the verified cold gains, X reuse9/12 time failures and unexecuted final
groups. Full plan/release qualification is not claimed. BLASTP's measured Wasm improvement does not extend to native; TBLASTX has small condition-dependent improvements, including 2.04% at Wasm
n1 for AP027280 self. Shared runtime linear-memory growth remains a separate
user-approved follow-up, with actual memory nonregression evidence retained. Frozen
release authorities and the previously deferred certification scope are unchanged.

The user explicitly requested commit and push, including all files required
to build the tested implementation. The Git handoff includes the two earlier
local runtime/profile prerequisite commits, all current implementation/build
changes, and the concise evidence records. The delivered revision is identified
by Git history. `commit-build-input-binding.json` verifies the index against the
qualified source; unrelated README prose and existing line-ending changes remain
in the worktree. The worktree had extensive pre-existing changes;
`change-review/inventory.json` and its separate patches identify this task's
source, tests, tooling and documentation changes against the captured initial
bytes. Generated evidence and the local handoff archive are reviewed separately.
