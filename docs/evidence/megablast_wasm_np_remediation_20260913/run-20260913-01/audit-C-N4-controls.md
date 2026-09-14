# Independent C-N4 control review

Reviewer: `ncbi_parity_auditor` (`/root/audit_megablast`), 2026-09-13. Read-only; no edits, builds, tests or measurements.

The completed `C-N4-main-controls` group contains 86 raw outputs: two oracles, twelve diagnostic runs and 72 cold runs including warmups. The reviewer independently compared all actual output bytes with the oracle, recomputed all six five-pair medians, ranges, maximum RSS values and engineering decisions, and checked all 36 A/B pairs for identical search arguments, alternating order and non-overlapping monotonic intervals.

Source SHA256 `4fa1b9a4510078c42f6bfe64fcce2555b168e408d8b023c29117f6ecc67e66f3`, six artifacts, the common ext4 runners, Node and fixture hashes agree with metadata. All six conditions pass the specified time and RSS nonregression gates. LC native n1 is 40.61 ms slower (8.27%), within the explicit 50 ms allowance; this is not a speedup claim. Other conditions also pass.

At this review, `main-n8-controls` was running and `task-controls` had not started. Those groups and final adoption are outside this completed review scope.

## Completed main n8 review

The subsequent review independently checked all 58 actual output files from `C-N4-main-n8-controls` against the official oracle, all four five-pair medians/ranges/maximum RSS decisions, 24 alternating A/B pairs and their non-overlap, search arguments, and source/artifact/runner/Node/fixture identities. No exclusions or defects were found. All four nonregression and RSS decisions pass.

LC native n8 is 0.319462→0.321204 s; threaded n8 is 0.707173→0.680013 s (3.84% improvement, below the primary 5% target). AP native n8 is 6.303008→4.485849 s; threaded n8 is 10.691804→6.194712 s (42.06% improvement). Additional LC samples and final gates remain pending. Task controls were still running at this review.

## Completed task controls and additional LC review

The reviewer checked the remaining 105 actual raw outputs (90 task-control and 15 additional LC), all 42 A/B pairs, source and artifact identities, search arguments and common runners. All raw outputs match the official oracle and all seven time/RSS decisions pass. This completes independent review of 249 output files across the four control/confirmation groups.

Additional LC is 0.745443→0.749166 s, −0.4994% improvement (+3.72 ms), within nonregression. Both five-pair cohorts and their ten-pair pooled calculation are preserved. The pooled 7.0876% arithmetic is correct but does not establish a stable 7.09% improvement because the cohorts disagree. The user's explicit LC acceptance amendment permits candidate selection without a 5% short-input claim. Final integration, supported formats and reuse/linear-memory gates remain pending. No concrete blocker to this limited acceptance was found.
