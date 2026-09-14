# Independent C-N4 source review

Reviewer: `ncbi_parity_auditor` (`/root/audit_megablast`), 2026-09-13. Read-only; no edits/builds/tests/measurements.

The reviewer supports the actual stable-arena implementation against NCBI blast_hits.c2455–2535. Payload storage never grows or reorders after construction. Initial handles are unique; stable sort preserves their permutation and incoming ties with the existing comparator and explicit NULL-last order. Only the active prefix is sorted. There is no new index tie-break or global sort API change.

Both trim paths mutate the referenced slot after immutable keeper/removed references end. Both deletion paths take/drop the payload and its script immediately before shifting. For decremented count c and removed index p, copy_within moves[p+1,c+1) to p, including the valid empty p=c case; tail assignment eliminates the temporary duplicate. This is the old shift's exact final state. Final materialization takes each live slot once in handle order and asserts missing/duplicate live slots, preserving NULL filtering, trimmed-tail order and extra_start.

Tests cover repeated trim/delete in both passes with NULL tails, all last-slot combinations, stable equal order, untouched prefix-external tail, payload addresses/content, and the updated post-sort injection seam. The previous last-slot test gaps are now covered. The added NCBI references/snippets were checked without errors. Other production files and Cargo/build settings match B1. Final patch/identity and actual source agree: SHA256 `4fa1b9a4510078c42f6bfe64fcce2555b168e408d8b023c29117f6ecc67e66f3`.

This supports source equivalence only. Execution raw parity, timing/RSS and adoption remain separate gates.
