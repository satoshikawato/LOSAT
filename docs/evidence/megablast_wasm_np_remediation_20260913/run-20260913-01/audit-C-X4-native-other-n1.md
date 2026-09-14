# Independent audit: C-X4 remaining n1 controls

Reviewer: /root/audit_megablast (ncbi_parity_auditor), read-only.
Mela PASS / APself FAIL; C-X4 is not accepted.

Mela B median3.610508486s / C3.712206799s,2.82% slower, time/RSS PASS.
AP B median29.983601831s / C32.398969841s,8.06% slower. Delta2.415368s
exceeds allowance1.499180s by0.916188s. Five-sample ranges are disjoint:
B29.763869–30.283059s, C31.956876–33.052244s. CPU medians30.76/33.28s.
Both RSS comparisons pass.

All30 actual raw outputs (52,288,515 bytes) match their official oracles:
2 oracle,4 diagnostic,4 warmup,20 timed. All input/argv/artifact/source/runner/
usage/result identities, AB/BA order, n1 pool0 and monotonic non-overlap pass
independent review, with zero exclusions. Minimum process gap57.233ms.
Two realtime disagreements are retained; all wall times equal their monotonic
endpoint differences, BOOTTIME disagreement at most9.157 microseconds.
No correction or exclusion is warranted.

The remaining-steps ledger exit_status0 is the successful benchmark child exit,
not adoption PASS. The following engineering check failed and the wrapper
exited1. Native n8 and prepared X4 Wasm/capacity/route work remain unexecuted.
The alignment hypothesis is not established by this failed candidate cohort.
