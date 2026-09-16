# Round 03: TBLASTX n8 implementation and evidence

## Current result

I4 is integrated into production `linking.rs`. All 158 production build inputs
match the fully tested candidate. Root focused tests, clippy and format pass;
final independent evidence review found no material blocker.

- Ordinary Node cold n8 body: Mje **22.224 → 10.833 s (−51.25%)**;
  AP self **37.739 → 17.591 s (−53.39%)**.
- Candidate n8 is faster than candidate n1 on both primary inputs.
- All fixed Native/serial, compiled-module, same-instance reactor and browser
  controls pass their original time/memory guards. Browser Mje warmed n1 has
  a small 2.0–2.6% increase, which is retained in the report.
- 79 command conditions and four long-code4 conditions match reference bytes;
  802 scan cases and 352 complete group cases match across three targets, trace±.
- Rust suite: 627 passed / 3 ignored. Clippy, formatting, four ABI shapes,
  serial API recovery and browser component/lifecycle checks pass.
- The long-code4 reference is NCBI **-db with -query_gencode 4 -db_gencode 4**,
  using a DB made from AP027133; it is not NCBI local -subject output.
- The user explicitly removed 1.20 Native ratio as an adoption/completion
  condition. The decision centers on actual improvement over the baseline.

Full tables, methodology, known small regressions, scope and limits:
**[I4-RESULTS.md](I4-RESULTS.md)**.

## Reproducibility

- `I4-integration.json`: sole production edit and source binding.
- `I4-functional-index.json`: raw/state/API/browser/Rust records.
- `I4-integrated-checks.json`: focused checks after copying the tested source.
- `I4-*-desktop/`: all fixed samples and raw outputs, including warmups.
- Matching `*-environment/`: prospective background policy and observations.
- `work/I4/` and `work/I4body/`: builds, exact artifacts and source hashes.
- `I4-POLICY.md`: selection, measurement policy and latest user steering.

These are recorded normal-desktop comparisons, not isolated-machine proof.
Observed CPU limits cannot exclude I/O/cache/frequency interference or missed
short-lived processes. All limitations are stated in the result report.

## Preserved investigation history

- I1: cold n8 improvement, rejected for the fixed reuse regression.
- I2/I2b: corrected target cfg and prefix cursor; contaminated studies remain
  separately invalid, with numerical failures preserved.
- I3: rejected for the clean serial regression.
- I4: restores indexed Native/serial caller, preserving one macro body.
- `I4-serial-git-diagnostic` stays ineligible; it is never pooled.
- `DIAGNOSIS.md`, candidate policies and `PRE-INTEGRATION-STATUS.md` preserve
  prior findings. Lost round-02 measurements never become acceptance evidence.

NCBI authority: `c++/src/algo/blast/core/link_hsps.c:812–895`.
