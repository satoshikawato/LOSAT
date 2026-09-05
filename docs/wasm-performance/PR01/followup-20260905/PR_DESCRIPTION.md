Title: Revalidate PR01 baselines with preserved evidence

The original PR01 audit left 16 native/serial rows incomplete under its
120-second budget. Reuse verified artifacts and completed measurements,
and retry missing executions with separately declared 600-second manifests.
The reconciled 43-native/41-serial catalog now has 74 PASS and 10 TIMEOUT;
p06, p10 and d04 native/serial validation is complete. No completed-output
mismatch or new exception was found.

Preserve original manifests, raw evidence, canonical output, all pre-existing
source changes and baseline timing samples. Test-only scripts retain every
process and bind joined rows to exact commands and raw gates. A recorded
interruption switched the audit to at most three concurrent searches;
contended durations are excluded from performance comparisons. An additional
raw-gated serial BLASTP CPU profile supports the scoped Kappa/redo priority
without claiming allocation dominance or a measured optimization.

Validation: identity/preflight and retained-data revalidation passed; 17
focused negative tests passed; the read-only NCBI parity auditor verified the
new outputs, evidence join, profile attribution and scoped decision. Production
code, existing harnesses, canonical authority and old raw samples are unchanged.

PR01 remains NOT_PROVEN / NOT_ACCEPTED. Native/serial p11, d01, d02, d03 and
d05 still time out; p11's oracle and cold/warm baseline remain incomplete.
PR02 implementation cannot automatically depend on unaccepted PR01.
Native effective-thread observations constrain same-effective-N claims, while
lifecycle, browser and registered platform certification remain later gates;
the requirements are separated in REQUIREMENTS.md. No later PR was implemented
and no commit, push or publication was performed.
