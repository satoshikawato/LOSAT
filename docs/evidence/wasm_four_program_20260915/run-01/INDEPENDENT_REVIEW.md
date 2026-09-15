# Independent read-only review

Reviewer: repository `ncbi_parity_auditor` agent, `/root/parity_audit`.
Status: source/correctness and the first three final pairs independently reviewed. **Integrated all-program nonregression is unsupported (six time-guard failures).**

Subsequent user decision: **ADOPTED**, accepting the documented tradeoffs; see
[ADOPTION.md](ADOPTION.md). The review below records the auditor's findings
before that decision and is preserved without relabelling any check.

## Source and correctness

- The cutoff correction follows the normal NCBI local-subject CLI database-mode
  path and uses the existing global query-context search space. The word7
  regression's raw output matches the fresh NCBI oracle.
- The C0-to-integrated composition and current owned Rust files were reviewed.
  No new source blocker was identified for the declared valid-query scope.
- The statistically invalid TBLASTX query's exit-status discrepancy occurs in
  an unchanged guard before X2 linking. It remains an open baseline defect;
  the review does not authorize an all-input or release parity claim.
- N2's native allocation counts and ordered batch evidence support fewer
  allocation calls. Peak requested bytes increase by 90 / 26 bytes; this is
  not evidence of lower memory use or a Wasm allocation measurement.

## Initial individual measurements

The reviewer independently recomputed all 36 paired conditions and checked
432 cold-process raw outputs, artifact identities, commands, warmup counts and
alternating order. All initial time/RSS guards pass. P1 remains mixed; X2's two
main default-Node gains are 12.47% and 3.59%, so the 5% target is not met on both.
P2's single-query n8 cold gain is 5.13%; multi-query observations do not measure
that branch's query reuse. No individual effects may be added together.

## Final measurement protocol

The final conditions and P2 fixture selection match the declaration. Body
timing must be extracted from exactly one finite positive
`[BODY_SCOPE_SECONDS]` value in each eligible diagnostic run, independently of
whole-process timing. The extractor implements that check.

Reuse uses two independent AB/BA process sessions with one warmup and five
invocations per condition in each. Within a reactor session, n1/n8 share the
same per-case instance. Memory and process-lifetime peak RSS therefore cannot
be attributed to isolated thread-count runs. Failed finite sessions and
malformed results are retained, independent declared sessions continue, and
the phase returns failure if any case fails. No deadline or observation window
is extended. Incomplete prefixes cannot establish a memory plateau.

The separately declared TurboFan-only **cold** comparisons satisfy the plan's
host-flag distinction. Both versions use the same flags; default Node remains
the primary runtime. No TurboFan body/reuse/browser claim is authorized.

## Final three-pair audit after user stop

The auditor independently selected exactly repeats 1–3: 270 timed processes,
45 conditions × three A/B pairs. All recorded output hashes match the NCBI
2.17.0 oracle. Time guards pass 39/45; max-of-three RSS guards pass 45/45.
The six failures and absolute medians are in `integrated-three-pair-results.md`.
All ten TBLASTX conditions pass both guards; threaded n8 cold gains are 16.32%
and 14.40%. Integrated all-program nonregression/adoption is unsupported.

The 54 completed repeat-4 records are retained but excluded under the user's
new three-repeat limit. Body, reuse, final control timings and TurboFan were
canceled before execution; the protocol review above is not execution evidence.

Initial individual builds lack contemporaneous source hashes; later snapshots
must not be represented as build-time source bindings. Their artifact hashes
and raw records remain valid. Final C0/integrated source/artifact bindings are
clean. No historical rebuild was run after the user's stop instruction.
