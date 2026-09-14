# Remediation status — provisionally adopted

The user explicitly requested stopping further measurements and provisionally
adopting the current implementation on 2026-09-14. C-N4, C-P1, C-X2, megablast
parity and threaded-default changes are retained. The qualified C-NA64 native
configuration is now in root; all 589 guarded root files exactly match its
snapshot. See `user-provisional-adoption-20260914.json` and
`C-NA64-root-adoption.json`. The user subsequently authorized committing and
pushing this source together with its build prerequisites. See
`commit-build-input-binding.json` and `COMMIT_HANDOFF.md`.

| Scope | Recorded outcome |
|---|---|
| Megablast | Current three parity targets and 14 BLASTN task cases pass; automatic initial X-drop follows NCBI |
| BLASTN | Wasm n8 large input 46.37% shorter; short input effectively unchanged (+0.65 ms), with its explicit speed-floor waiver |
| BLASTP | Wasm n8 11.16%/11.40% shorter; n1 14.95%/14.64% shorter. Native 0.82–3.80% slower, within the fixed allowance |
| TBLASTX | Removes 4N final-sort payload clones. Native medians 0.03–2.55% shorter; AP self Wasm n1 2.04% shorter. Wasm n8 0.02–2.56% slower, within cold allowances |
| Threaded default | Standard distribution/checks use threaded artifacts. Serial remains explicit compatibility for current gbdraw AUTO/fallback consumers; gbdraw is unchanged |

All 37 completed cold conditions pass their recorded time/RSS allowances, with
actual raw output and independent reviews. The new native artifact has 14 fresh
conditions. Native format, 621 Rust tests (three pre-existing ignored), Clippy
and threaded builds pass. New command/reactor Wasm and four packaged JS files
are byte-identical to the previously tested integrated bundle. X n1 and n2/n4
raw/thread boundaries also pass.

**Provisional adoption is not a claim that every plan gate passed.** The original
X reuse cohort has 144 matching outputs and complete worker lifecycles, RSS and
budget 4/4 PASS, but time **9/12 FAIL**. Same-B1 time drift and overlapping bulk
read/hash work leave the attribution unresolved; no sample was removed or
corrected. The AP-only 12-job stage diagnostic completes, with both versions
about 25–27 seconds per call; it does not replace the failed three-input cohort.
The following 36-job interleaved diagnostic was stopped by user instruction
(exit130), retaining one complete parent process plus a partial next process.
The planned additional-five cohort was never run.

The six later groups remain unexecuted: final X thresholds, long code-4,
daily routes, route validation, capacity diagnostics and capacity summary.
The separate final native megablast recheck was also not run. Earlier protected
fixtures/build/contract evidence remains valid only for its recorded artifacts;
these pending groups are not marked PASS. No further measurement is authorized
by this handoff.

The user previously authorized keeping shared-runtime linear-memory growth as
a separate issue while retaining memory nonregression. Original N/P reuse and
fixed16 plateau failures remain. The allowed N/P additional-five time cohort
passes additional and pooled checks; pooled ten values are separate instance
cohorts. Observed N/P post-call live maxima/blocks are equal, with original RSS/
budget checks passing. One P132 fixed-window linear maximum is655,360 bytes
larger. X allocator live bytes are unmeasured/N/A; its new owners are local to
the search. Final capacity checks were stopped with the remaining scope.

Frozen release expected bytes, registered native platform fingerprints and the
46 deferred cases are unchanged. Corrected Sakai output differs from frozen
LOSAT bytes; formal recertification needs a separately reviewed authority
version. Current comparison success and user-directed provisional adoption do
not constitute full release certification.
