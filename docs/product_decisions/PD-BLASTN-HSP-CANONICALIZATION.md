# Product Decision: BLASTN HSP edit-script canonicalization

- Decision ID: `PD-BLASTN-HSP-CANONICALIZATION`
- Version: 1.0
- Date: 2026-08-29
- Status: Accepted; implementation pending

## Scope

This decision governs BLASTN HSPs that compare equal under the NCBI
common-endpoint comparator but have different normalized `gap_info` edit
scripts. It defines the contract for a later implementation change. It does not
change current LOSAT runtime behavior.

## Problem

Distinct score-equivalent traceback paths can have the same query context,
query and subject endpoints, and raw score while retaining different edit
scripts. Those scripts can produce different identity, alignment-length,
mismatch, and gap-open values even though the HSP coordinate and score keys are
unchanged.

Leaving the survivor dependent on an unspecified equal-element sort order makes
native, Wasm, and cross-platform results potentially disagree. This is also
user-visible downstream: gbdraw comparison filtering applies both identity and
alignment-length thresholds directly ([`gbdraw/io/comparisons.py`, lines
29-45](https://github.com/satoshikawato/gbdraw/blob/fd286c5e664e31f03e986f96c06c14d70580a963/gbdraw/io/comparisons.py#L29-L45)).

## NCBI source behavior

In the official NCBI BLAST+ 2.17.0 source,
`c++/src/algo/blast/core/blast_hits.c:2268-2387`,
`s_QueryOffsetCompareHSPs` and `s_QueryEndCompareHSPs` compare context,
endpoint coordinates, and raw score. Neither comparator examines `gap_info`,
and each returns equality after those fields tie.

`Blast_HSPListPurgeHSPsWithCommonEndpoints` at
`blast_hits.c:2455-2535` passes those comparators to C `qsort` and keeps the
first HSP in each common-start or common-end group. C `qsort` does not define a
stable relative order for comparator-equal elements. The release source
therefore does not select one canonical edit script for this tie class.

LOSAT currently mirrors that comparator field set. Native LOSAT routes the
sort through platform C `qsort`, while wasm32 uses Rust sorting. Consequently,
the current source does not define one cross-platform survivor when only the
edit script differs.

## Observed compatibility consequence

A fresh comparison against NCBI BLAST+ 2.17.0+ contains 14 cases: 13 are
byte-exact and one, `Sakai.MG1655.megablast`, is a value difference. That case
has the same row count, coordinate keys, E-values, and bit scores, but five rows
have different edit-script-derived values; two of those rows also differ in
percent identity. Static source analysis classifies the residual as
`CANONICALIZATION_DECISION_REQUIRED`.

Those five rows demonstrate the consequence, but they are not an implementation
oracle. A precompiled binary's arbitrary survivor among source-equal elements
does not add an ordering rule to the 2.17.0 source contract.

## Decision

LOSAT shall define a deterministic total order for normalized BLASTN edit
scripts when all fields used by the NCBI common-endpoint comparator are equal.
This rule is a **LOSAT deterministic compatibility extension**. It is not a
claim that NCBI BLAST+ 2.17.0 source defines the ordering.

The extension exists to provide native reproducibility, Wasm reproducibility,
cross-platform reproducibility, and stable downstream gbdraw behavior.

## Canonical owner

The normalized edit-script representation is the canonical secondary key.
Derived values such as percent identity, alignment length, mismatch count,
gap-open count, E-value, and bit score remain outputs of the selected script;
they must not select that script. Biological identifiers, case names,
accessions, and coordinates copied from observed residual rows are also
excluded from the key.

## Required ordering semantics

The future implementation shall use these semantics exactly:

1. Represent an edit script as the complete run-length-encoded sequence of
   `GapEditOp` values. For the comparison key, remove zero-length operations
   and coalesce adjacent operations of the same kind without integer wrap.
   Preserve all other operation boundaries and their stored order.
2. Include script presence in the key: `None` has presence rank `0` and
   `Some(...)` has presence rank `1`. This distinguishes an absent script from
   a present empty script.
3. Map each normalized operation to `(operation_kind_rank, operation_length)`.
   Operation-kind ranks are `Sub = 0`, `Del = 1`, and `Ins = 2`.
4. Compare operation lengths as unsigned integer counts in ascending order.
5. Compare the complete tuple vectors lexicographically in ascending order.
   At the first unequal tuple, the smaller tuple sorts first. If one vector is
   an exact prefix of the other, the shorter vector sorts first.
6. Use the stored normalized operation order. Do not sort the operations
   inside a script.
7. For reverse-strand HSPs, compare the internal stored script in the same
   order. Do not reverse the vector and do not exchange `Del` and `Ins`.
   Output-coordinate orientation is not part of this key.
8. Append this key only after every NCBI-defined field in the common-start or
   common-end comparator has tied. Apply it in both common-start and common-end
   sort passes, before their survivor-selection loops. The lower canonical key
   sorts first and is therefore the survivor considered first by the existing
   purge logic.
9. If the normalized representations, including presence, are structurally
   identical, the secondary comparison returns equality. Instance identity or
   allocation order is not a further key because either survivor has the same
   canonical script representation.
10. Native and wasm32 implementations must use these exact semantics. Platform
    sort stability must not decide between structurally different scripts.

The ascending ranks and direct tuple-vector comparison are intentionally
representation-derived. They do not encode preferences such as "fewer gaps is
better" and were not chosen to reproduce the five observed NCBI-binary rows.

## Native/Wasm requirement

The later implementation PR must wire one shared semantic comparator into both
native and wasm32 common-endpoint ordering. Tests must cover both purge passes,
reverse-strand stored order, prefix ordering, all operation-kind ranks, script
presence, identical normalized scripts, and native/Wasm output agreement.
Until that implementation lands, no deterministic cross-platform survivor
claim may be made for this tie class.

## Non-goals

- No production comparator, sort invocation, survivor selection, traceback,
  edit-script construction, metric derivation, or output behavior changes in
  this decision PR.
- No special case for `Sakai.MG1655.megablast` or any of its rows.
- No attempt to infer undocumented intent from one NCBI binary build.
- No change to source-defined BLASTN parity requirements.
- No gbdraw code or cross-repository runtime dependency.

## Release-gate impact

Source-defined BLASTN cases continue to require exact parity. For HSP ties that
the NCBI source leaves equal, release evidence shall require the deterministic
LOSAT canonicalization defined here and identical native/Wasm semantics. A
particular precompiled NCBI binary's arbitrary equal-element survivor is
diagnostic evidence, not the source-defined acceptance contract.

This document alone does not clear the release gate: the extension remains
implementation-pending, and the known 13-exact plus 1-value-difference result
is expected to remain unchanged in this PR.

## Reversibility / future revision process

Changing the canonical key is a compatibility change. A revision must update
this document's version, explain the reason, add representation-level component
tests and native/Wasm evidence, and identify downstream compatibility impact.
It must not silently rewrite expected binary rows or make derived statistics
the canonical owner.
