# Independent additional N/P reactor timing confirmation

The required read-only auditor checked all four processes / 96 jobs and
recalculated the actual sample medians/ranges. **Additional 8/8, combined 8/8,
and process RSS 2/2 pass the unchanged gates.** The original reuse failure is
unchanged. All actual raw bytes (128546304 bytes read), argv/job templates,
source/artifact/Node/runner bindings, predeclaration times, AB/BA order, exact
process endpoints and every worker ID/event order/normal exit agree.

| Short LC reactor | Additional median difference | Combined ten-value difference | Allowance |
|---|---:|---:|---:|
| session 0 | +39.185704 ms | +32.375236 ms | 50 ms |
| session 1 | +6.339746 ms | +7.078536 ms | 50 ms |

The other six conditions pass both additional and combined evaluations.
Combined values pool two separate instance lifetimes; they are not ten
consecutive calls in one instance. No further timing confirmation is allowed.

Process peak RSS is 563924992 -> 558354432 bytes for session 0 and
578129920 -> 557273088 bytes for session 1. All comparisons pass.

Four complete memory series retain late growth within the last three calls:
session 0 baseline large N, session 0 candidate P078, session 1 candidate large
N, and session 1 candidate P078. The auditor checked every series. This timing
acceptance does not imply a memory pass. The original two-input memory proposal
was never executed; expanding the first fixed memory diagnostic to the full
four-input workload, 256 untimed jobs and a final-eight plateau was reviewed
before execution. Original and additional six-call observations remain.

Nineteen jobs (16 measured) retain realtime-clock disagreement. Maximum
difference +697.443064 ms is P132 repeat 1 in session1-baseline/process/output.txt.
No records were excluded or corrected.

The auditor supports a separate timing-confirmation disposition and the narrow
ownership note. Fixed memory observations, capacity evidence and final
independent adoption disposition are still required. The original stopped
pipeline and its failed N/P reuse step remain preserved.
