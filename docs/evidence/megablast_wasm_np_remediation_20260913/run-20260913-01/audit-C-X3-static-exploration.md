# Independent static review: C-X3 exploration

Reviewer: /root/audit_megablast (ncbi_parity_auditor), read-only.
The extracted scan has no identified semantic or borrowing defect. This
supports static equivalence for exploration, not adoption or a speedup claim.

NCBI link_hsps.c827–861 visitation is retained: save the current helper index,
decrement, apply next_larger when the score cannot improve, test strict score
and trim bounds, copy the four selected values. The caller reaches this scan
with num/xsum/link still at their initial values; only h_sum may have changed.
The f64 value is copied, and the subsequent left-associative arithmetic,
helper updates and next_larger construction876–885 remain in the caller.
Immutable slice borrows introduce no pool mutation, allocation, payload clone
or unsafe code. Only linking.rs differs in the589-file frozen source tree.
Reviewed linking.rs SHA256:
5d02681fbd99d513a36bc0cefbec42eeb080b89f384f8d51849ed90dd84a735d.

The new1050-combination test checks empty scans, ties, trim equalities, initial
sums and the selected xsum bits. It currently uses MAX sentinels/fallback1,
so it does not cover the production zero-sentinel jump-to-0. This is a test
coverage gap, not an identified implementation bug. If promising, a new
revision should add that case with negative sums/initial values and describe
register/stack isolation as an intent, not a proven performance cause.
The current exploratory snapshot remains immutable.
