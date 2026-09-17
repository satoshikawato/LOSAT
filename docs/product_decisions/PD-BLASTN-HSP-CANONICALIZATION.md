# Product Decision: BLASTN common-endpoint equal-HSP policy

- Decision ID: `PD-BLASTN-HSP-CANONICALIZATION`
- Version: 1.3
- Date: 2026-09-17
- Status: Accepted by the maintainer during v0.1.0 distribution implementation

## Current decision

Keep the source-correct automatic megablast X-drop correction in
`f5955c5952998e50c1262186212d8cfc1527eb1d`. Retire the active five-row Sakai
output exception. All 14 declared Linux BLASTN cases now require exact raw
NCBI output; current LOSAT native/Wasm outputs must also satisfy Gate A version
2. This is an authority revision, not completed recertification.

Fresh reproduction at `89a85d6664a44c2dd41fe738158222c88c15b1e8` and independent
source review establish an upstream X-drop cause, not an unexplained survivor
choice. NCBI `api/blast_options_local_priv.cpp:49-54` zero-initializes options;
`api/blast_nucl_options.cpp:171-174` sets only the megablast window; and
`core/blast_parameters.c:380-383` selects the subject cutoff for zero X-drop.
See [the exact raw diff and review](../evidence/v010_distribution_20260917/AUTHORITY_REVIEW.md).

The general NCBI comparator semantics and component tests below remain valid.
They do not authorize the retired fixture discrepancy. No additional sort key,
normalization, runtime compatibility branch, or other output exception is added.
Gate B remains an independent exact platform fingerprint check. Old PR5/PR6
measurements remain historical and must not certify the corrected implementation.

## Historical Version 1.2 (superseded fixture classification)


- Decision ID: `PD-BLASTN-HSP-CANONICALIZATION`
- Version: 1.2
- Date: 2026-09-02
- Status: Accepted; Version 1.0 canonicalization requirement withdrawn for
  v0.1.0

## Scope

This decision governs BLASTN HSPs that compare equal under the NCBI
common-endpoint comparator but have different normalized `gap_info` edit
scripts. It revises the Version 1.0 product policy without changing LOSAT
runtime or test behavior.

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
therefore does not select one unique edit-script survivor for this tie class.

The static source audit also found LOSAT's greedy traceback, edit-script
construction, gap reduction, post-traceback re-evaluation, identity/length
counting, and tabular derivation semantically aligned with the relevant NCBI
2.17.0 source path for the supported case. The remaining observed discrepancy
is at this common-endpoint survivor boundary; the audit found no missing NCBI
secondary ordering rule to port.

## Decision

LOSAT v0.1.0 shall preserve the NCBI BLAST+ 2.17.0 source semantics for
common-endpoint HSP ties.

When NCBI's comparator treats distinct score-equivalent edit scripts as equal,
LOSAT does not introduce an additional canonical ordering solely to reproduce
a particular precompiled NCBI binary. Such ties are source-underdetermined and
may have more than one source-compatible edit-script representation.

A LOSAT-specific deterministic canonicalization may be reconsidered only if
supported native/Wasm execution is demonstrated to produce materially
inconsistent behavior that matters to the supported product contract. A
theoretical possibility of inconsistent ordering is not sufficient evidence
for adding that policy.

## Compatibility contract

### NCBI source compatibility

NCBI source permits multiple survivor representations for comparator-equal
edit scripts. LOSAT may therefore produce any representation that the source
permits without violating the source-level contract. This allowance is limited
to the demonstrated common-endpoint tie class; it does not authorize any other
algorithm or formatting difference.

### Precompiled-binary byte parity

Source-defined behavior continues to require exact parity with the official
NCBI BLAST+ 2.17.0 binary. For a source-underdetermined tie, one precompiled
binary's chosen edit-script representation is diagnostic evidence rather than
a normative source contract.

For the bounded Windows/macOS release gate,
[`PD-NCBI-PLATFORM-VARIANCE`](PD-NCBI-PLATFORM-VARIANCE.md) records exact
platform-local official-binary fingerprints independently from the frozen
LOSAT canonical bytes. This decision's LOSAT source-underdetermined authority
does not require those registered NCBI raw bytes to become LOSAT output and
does not allow a semantic diagnostic to rescue an unknown native fingerprint.

### Determinism

No cross-platform deterministic edit-script survivor behavior is claimed
without supporting evidence. LOSAT will monitor and test the supported native
and Wasm executions. If future evidence demonstrates a material inconsistency
between supported executions, this decision must be reopened before adopting a
LOSAT-specific ordering rule.

### Downstream gbdraw impact

Distinct source-compatible scripts can change percent identity and alignment
length. Those fields matter near the gbdraw comparison thresholds
([`gbdraw/io/comparisons.py`, lines 29-45](https://github.com/satoshikawato/gbdraw/blob/fd286c5e664e31f03e986f96c06c14d70580a963/gbdraw/io/comparisons.py#L29-L45)).
That downstream sensitivity is a reason to monitor and test supported
executions, but it does not by itself justify inventing an ordering rule before
an actual supported-platform incompatibility is observed.

## v0.1.0 release gate

For the current 14-case supported BLASTN manifest:

- The 13 source-defined cases require exact official NCBI BLAST+ 2.17.0 binary
  parity.
- The known source-underdetermined case may be accepted only when evidence
  shows the same HSP row count, the same HSP membership and coordinate keys,
  the same E-values, and the same bit scores, with differences limited to the
  source-underdetermined edit-script-derived fields currently observed:
  percent identity, alignment length, mismatch count, and gap-open count.

This exception does not permit coordinate changes, HSP-set changes, E-value
changes, bit-score changes, unrelated formatter differences, or unrelated
algorithm differences. It must not be generalized beyond the demonstrated
equal-HSP tie class.

## Existing component evidence

The component tests added with Version 1.0 remain useful and must not be
removed or weakened:

- `equivalent_scripts_keep_endpoints_and_score`
- `common_endpoint_comparator_ignores_edit_script`
- `source_compatible_first_survivor_accepts_either_equal_order`
- `source_compatible_purge_never_synthesizes_edit_script`

They establish that different edit scripts can be score/endpoint equivalent,
the NCBI-compatible comparator does not inspect the edit script, either
incoming equal order is source-compatible, and purge does not synthesize a
third edit script. They do not establish or require canonicalization.

## Non-goals

- No comparator, sort, purge, traceback, edit-script, metric, formatter, or
  other runtime behavior change.
- No Rust, Python, test, manifest, fixture, expected-output, workflow, BLASTP,
  or TBLASTX change.
- No special case for `Sakai.MG1655.megablast` or any of its rows.
- No weakening of exact parity for source-defined behavior.
- No claim that supported native/Wasm executions have already demonstrated one
  deterministic survivor.

## Revision history

Version 1.2 cross-references the bounded native-NCBI platform authority. It
does not change the Version 1.1 LOSAT source-underdetermined policy.

Version 1.0 correctly recorded the NCBI source-level ambiguity, then proposed a
LOSAT-specific representation-derived secondary total order using script
presence, operation-kind ranks, operation lengths, lexicographic comparison,
stored reverse-strand order, and application before both purge passes. That
canonicalization proposal is a withdrawn, non-normative historical record for
v0.1.0.

Version 1.1 preserves the source audit and withdraws the additional LOSAT
ordering requirement. Future canonicalization requires a new revision backed
by evidence of a material incompatibility in supported execution.
