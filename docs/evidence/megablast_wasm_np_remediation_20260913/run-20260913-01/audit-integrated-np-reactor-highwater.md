# Independent audit: fixed reactor memory observation

The required read-only `ncbi_parity_auditor` independently confirmed
`FAIL_FIXED_OBSERVATION`: the predeclared final eight observations are constant
for 10 of 16 series and increase in six. All 256 jobs were untimed; exclusions: 0.

All 342,790,144 output bytes match the original official oracles. Inputs, argv,
artifacts, Node, runners, source binding, ABBA ordering, monotonic non-overlap,
worker IDs, event ordering and exit status were checked against actual records.

| Session | Version | Case | Late increase (zero-based repeat) |
|---|---|---|---|
| 0 | baseline | P078 | 15: 69,206,016 → 69,271,552 |
| 0 | candidate | P078 | 13: 69,074,944 → 69,140,480 |
| 0 | candidate | large BLASTN | 10: 132,055,040 → 137,297,920 |
| 1 | candidate | P078 | 15: 69,140,480 → 69,206,016 |
| 1 | baseline | short LC | 12: 46,530,560 → 46,858,240; 14: → 47,120,384 |
| 1 | baseline | large BLASTN | 15: 138,674,176 → 143,917,056 |

Every process remains within its RSS budget. All 51 job-level realtime clock
disagreements remain recorded and untimed; the largest is +0.840604444 seconds
(session 1 baseline P132 repeat 11). Outer monotonic/BOOTTIME disagreement is at
most 5.631 microseconds. Job-level BOOTTIME is not available.

The original reuse failure, separately accepted supplemental timing and this
fixed memory failure remain separate records. Baseline failures do not exempt
candidate failures. Source ownership and selected Vec capacities do not prove
that the observed growth consists exclusively of allocator high-water pages.
Direct live allocation accounting is the next diagnostic; it does not rewrite
the failed gate or justify extending the observation count until it passes.
