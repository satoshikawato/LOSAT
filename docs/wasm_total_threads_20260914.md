# Total thread count and gbdraw Wasm handoff

The user requires `-num_threads N` to include the calling thread. For N > 1,
LOSAT uses the caller as Rayon worker zero and creates exactly N-1 child
threads. Thus n8 is one participating caller plus seven children. n1 remains
serial. This supersedes the September 13 dedicated-N-workers-plus-caller
contract; the historical reports retain their original results.

`LOSAT/src/utils/threading.rs` remains the single scheduler owner for native,
threaded command Wasm, and threaded reactor Wasm. A custom Rayon spawn handler
captures slot zero and starts only slots 1..N. The caller runs the ordinary
`ThreadBuilder::run` loop, whose startup hook executes the search. Dropping the
pool ends that loop and clears its TLS registration; scoped children are joined
before return. A user panic is caught inside the startup hook and resumed only
after cleanup. The leaked registration from `use_current_thread()` is not used.
Entry from an already registered Rayon worker is rejected before construction.

The private erased callback remains at its original address until the caller
loop and all child joins finish. Only slot zero dereferences it, exactly once;
other workers retain the hook without accessing the borrowed state. Partial
spawn failure never enters the caller hook and still joins existing children.

NCBI `c++/src/algo/blast/api/prelim_stage.cpp:145,172-180` is the reference for
the N-worker search context and completion ordering. Parent participation is
the user's scheduling requirement, not a claim about NCBI's caller thread.
Rayon-core 1.13.0 `registry.rs:282-315,496-512,687-693,910-938` supplies the
construction, same-pool execution, and TLS cleanup mechanics. Search candidates,
pruning, reductions, statistics, and output formatting are unchanged.

gbdraw retains its existing N-1 child preparation, input ownership, output
collection, AUTO selection, and serial fallback. The initial September 14
candidate and proposed consumer Worker changes were not adopted. Only the two
final Wasm files and their adjacent build notes are to be handed over.

This is `IMPLEMENT_EXISTING_AUTHORITY`: the user's explicit total-thread budget
and gbdraw's supported search behavior select the outcome. No new semantic
owner, execution path, compatibility branch, or product option is introduced.
There is no increase in owner/path/compatibility excess. Previous consumer
assets and notes are backed up for rollback.

Build inputs, binaries, and verification records are retained in
[`artifacts/gbdraw-wasm-total-threads-20260914`](../artifacts/gbdraw-wasm-total-threads-20260914/).
Verification passed with the built artifacts:

- Native scheduler: 3 tests, including caller participation, changing thread
  counts, errors, panic recovery, and TLS cleanup. Clippy (all targets/features,
  warnings denied), focused formatting, and 22 Python evidence tests passed.
- NCBI BLAST+ 2.17.0 comparison/command gate: 317 execution records including
  raw comparisons, expected failures, thread-budget checks, and artifact guards;
  zero format failures. Native and command Wasm n1/n2/n4/n8 are covered.
- The same threaded reactor instance: 144 calls, including changing thread
  counts, partial spawn failure and recovery, output formats, and 24 stress
  repetitions with stable guest memory in the final 12 runs.
- Chromium 149.0.7827.55 using unchanged gbdraw Worker sources: 23 checks across
  non-isolated serial fallback and isolated threading, repeat execution,
  cancellation/recovery, and BLASTP n8. All five BLASTN/megablast/TBLASTX/
  genetic-code-4/BLASTP fixtures matched native bytes in serial and threaded n2;
  BLASTP also matched at n8. External requests and page errors: zero.
- An actual n8 command reports pool_threads=8 and caller_participates=true;
  host spawn_attempt/spawned/ready/exited counts are each seven.

These local checks do not expand frozen release certification and make no speed
claim. The oracle is the installed NCBI 2.17.0 build; its exact binary hashes and
Node 26.8.2 runtime metadata are recorded by the gate.

The final serial and threaded command binaries were installed into
`gbdraw/gbdraw/web/wasm/losat/` with their notes. The adopted paths passed seven
additional BLASTP browser checks (serial/threaded repeat, n8, cancellation and
recovery), with no external requests or page errors. `installed.json` records
previous and installed hashes. Both existing gbdraw Worker source hashes are
unchanged.

The independent [read-only audit](../artifacts/gbdraw-wasm-total-threads-20260914/independent-audit.md)
verified source/artifact hashes, 189 raw candidate/oracle equalities, N-1 child
counts, reactor recovery, and both browser phases; it found no blocking defect.
