# Adoption decision: I4

**Adopt I4.** The production change is integrated in linking.rs, with all 158
build inputs matching the tested candidate. Independent evidence review and
post-integration focused tests/clippy/format pass.

Ordinary Node cold n8 body improves **51.25% on Mje** and **53.39% on AP self**.
Candidate n8 is faster than candidate n1 on both fixtures. All specified
Native/serial/module/reactor/browser guards pass. Small measured regressions
remain visible, including browser Mje warmed n1 at approximately +2.0–2.6%.
The user explicitly excludes a 1.20 Native ratio from adoption/completion.

See [all results, scope and limitations](I4-RESULTS.md),
[independent review](I4-INDEPENDENT-REVIEW.md), and
[integration source binding](I4-integration.json).

Production: one Rust file; search semantics, thread counts and output remain
unchanged. I1/I3 failures and ineligible/contaminated studies are preserved,
not pooled or retroactively promoted. Prior H1 remains the baseline.
No deployment, publishing, tag, push or external application update occurred.
