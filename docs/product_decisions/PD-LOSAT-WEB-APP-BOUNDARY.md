# Product Decision: LOSAT Web application boundary

- Decision ID: `PD-LOSAT-WEB-APP-BOUNDARY`
- Version: 1.0
- Date: 2026-09-29
- Status: Accepted by the maintainer on 2026-09-29 (plan decisions DW-1 and DW-12 in
  [`docs/losat_web_gui_plan.md`](../losat_web_gui_plan.md)).

## Scope

LOSAT Web is the in-browser application described in
[`docs/web/losat_web_design_v0.1.md`](../web/losat_web_design_v0.1.md) and planned in
[`docs/losat_web_gui_plan.md`](../losat_web_gui_plan.md). This decision governs:

- all code under `web/` (the Rust adapter crate `web/adapter/` and the TypeScript
  application `web/app/`);
- the engine entry points that exist to serve that code (the program-independent
  local-search entry in `LOSAT/src/api/local_blast.rs` and the fan-out of one search
  result to several output writers).

It does not change the rules in the root [`AGENTS.md`](../../AGENTS.md) for code under
`LOSAT/`. Any change under `LOSAT/`, including one made for LOSAT Web, still requires
NCBI source authority, NCBI reference comments and the parity gates.

## Authorities

- NCBI BLAST+ source remains the only authority for search behavior and for the
  compatibility outputs outfmt 0, 6 and 7 (`AGENTS.md`,
  [`PD-PURE-RUST-RUNTIME-AUTHORITY`](PD-PURE-RUST-RUNTIME-AUTHORITY.md)).
- The design document, as amended by the plan's maintainer decisions (plan §0.4), is
  the authority for application behavior: screens, jobs, storage, extraction and
  exports.

## Decision

1. Application features that have no NCBI equivalent (job queue, candidate tray,
   temporary files, sequence extraction, session files, CSV/JSON/report exports,
   offline caching and usage measurement) are permitted, but only under `web/` and
   only while the invariants below hold.
2. Invariants:
   1. **Engine input.** The engine receives exactly what an equivalent CLI invocation
      receives: an argv parsed by the same clap parser as the CLI, and FASTA bytes made
      only of whole original records, in their original order, as stored in the run
      snapshot (a single newline may be added between concatenated files that do not
      end with one). The application never changes sequences, options or statistics,
      and never splits one search into independent searches whose results are joined.
      The exact input bytes of every run can be exported, so the run can be reproduced
      with the CLI.
   2. **Compatibility outputs.** outfmt 0, 6 and 7 are produced only by the engine's
      formatters. The application stores and exports them byte for byte and never
      rewrites, filters or regenerates them.
   3. **BLAST-defined values.** Code under `web/` does not compute or format scores,
      E-values, identities, coverages, NCBI number formats or alignment text. Values
      shown to users come from engine output or from records the engine provides.
   4. **Application outputs.** CSV, JSON, report and session files are labelled as
      LOSAT Web formats and never claim identity with NCBI output.
   5. **Research data stays local.** Sequences, headers, file names, results and notes
      are never sent over the network.
3. **Subject retention.** For LOSAT Web, the engine may keep registered subjects and a
   warm Wasm instance between runs. This supersedes, for LOSAT Web only, the constraint
   "do not mix a persistent worker pool or a mutable cross-search cache" in
   [`docs/wasm_remaining_work_implementation_plan_20260915.md`](../wasm_remaining_work_implementation_plan_20260915.md)
   (line 86), which was written for that plan's performance experiments. The condition
   is that every output of consecutive runs is byte-identical to a fresh CLI run with
   the same argv and input.
4. **Engine plumbing authorized for LOSAT Web.** The following changes under `LOSAT/`
   are permitted although NCBI has no identical API, because they only route inputs and
   outputs and never change search behavior or output bytes. Each still follows the
   root `AGENTS.md` (NCBI reference comments for the behavior it routes, parity gates,
   independent audit):
   1. one program-independent core entry per program (`run_local` in
      `LOSAT/src/api/local_blast.rs`) through which the CLI, web ABI v1 and web ABI v2
      all run, available on every target (not only wasm32);
   2. handing one search result to several output writers, with format-dependent
      options resolved per format, and failing with an explicit error if the resolved
      search options differ between formats;
   3. a diagnostics writer that receives the warnings the CLI writes to stderr, each
      once;
   4. a callback that receives the final `PairwiseHit` list with stable HSP identifiers;
   5. a formatter observer that reports when the row or section of an HSP starts and
      ends, without changing the written bytes;
   6. keeping registered subjects and a warm instance between runs (item 3).
5. **Web ABI v1.** `losat_web_*` in `LOSAT/src/web_api.rs` stays frozen, except for
   fixes that make an unsupported request fail fast instead of returning a wrong
   result. Such a fix can turn a request that used to succeed into an error for gbdraw;
   it is recorded in the commit message and the gate record. Removing v1 is a separate
   decision, taken when gbdraw moves to v2.
6. Rules for writing code under `web/` are in [`web/AGENTS.md`](../../web/AGENTS.md).
   NCBI reference comments are not required for application code under `web/`,
   because invariant 2.3 keeps BLAST behavior out of it.

## Compatibility contract

- Identity gates (plan §6): the engine entry must reproduce the CLI bytes (V-NAT), the
  Wasm adapter must reproduce the frozen native bytes (V-ABI), and the browser must
  reproduce them through the real application path (V-BR).
- A verification badge may say that a run matches a certified profile only when the
  program, task, options and runtime path are inside a profile certified for that
  engine build and runtime path. Other runs are labelled "Engine-supported, outside
  certified profile". Application outputs never inherit an identity label.

## Change policy

Any change to the invariants requires a new version of this decision and maintainer
acceptance. An application feature that would need to break an invariant must instead
be implemented in the engine under the root `AGENTS.md` rules, or not at all.

## Non-goals

- It does not authorize new engine behavior or any NCBI runtime, build, FFI or
  subprocess dependency.
- It does not certify the browser runtime; certification comes only from the gates.
- It does not change the existing web ABI v1 (`losat_web_*` in `LOSAT/src/web_api.rs`),
  which gbdraw uses, beyond the fail-fast fixes in decision item 5.

## Revision history

- 0.1 (2026-09-29): proposed.
- 1.0 (2026-09-29): accepted. Defined the engine input as whole original records;
  replaced "one search per option set" with an explicit error.
