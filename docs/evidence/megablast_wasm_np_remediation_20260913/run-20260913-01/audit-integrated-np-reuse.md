# Independent integrated N/P reuse audit: original cohort fails

The read-only `ncbi_parity_auditor` independently recalculated all 20 parent
records and eight processes / 192 jobs. All actual raw bytes, oracle bindings,
worker lifecycles and process RSS checks pass. Adoption checks are **time 15/16,
process RSS 4/4, actual reactor plateau 14/16**. Command's 16 plateau booleans are
not-applicable sentinels, not observations of a persistent command instance.

| LC same-instance reactor | Baseline median | Candidate median | Difference | Result |
|---|---:|---:|---:|---|
| session 0 | 0.683540788 s | 0.735930843 s | +52.390055 ms | FAIL: exceeds unchanged 50 ms gate by 2.390055 ms |
| session 1 | 0.744808753 s | 0.741441646 s | −3.367107 ms | PASS |

Session 0 ranges are baseline 0.659560–0.708803 s and candidate
0.702348–0.807034 s. The short-input speed-floor waiver does not waive this gate.

Two candidate session 0 reactor memory series fail the original final-three
plateau check (repeat 0–5, bytes):

- Large N: `[121569280,126812160,126812160,126812160,132055040,132055040]`.
- P132: `[85131264,86704128,86704128,86704128,86835200,86835200]`.

Session 1 does not reproduce the late growth: large N reaches the same
132055040 bytes at repeat 2 and then remains constant; P132 is constant from
repeat 1 at 86966272 bytes. Neither this observation nor the original late
growth alone determines whether memory leaks.

| Session / mode | Baseline process peak RSS | Candidate process peak RSS |
|---|---:|---:|
| 0 / command | 1987371008 | 1726664704 |
| 0 / reactor | 562016256 | 555483136 |
| 1 / command | 1955295232 | 1651683328 |
| 1 / reactor | 573571072 | 555659264 |

The audit checked four inputs, outfmt 6, n8, explicit BLASTN task, NCBI 2.17.0+
oracle bytes, frozen source, Node 24.21.0 with no extra flags, unchanged common
host/imports, both command/reactor artifacts, and every job's argv. AB/BA order
and exact process monotonic endpoints have no overlap. All worker IDs, event
order, and exit codes are correct. No exclusions or corrections were made.

Thirty-nine jobs (32 measured) retain realtime-clock disagreement. The largest
difference is +647.244389 ms for session 0 candidate command P132 repeat 5.
All parent processes pass monotonic/boottime agreement, with maximum difference
5.134 microseconds. Individual jobs do not record boottime; they use
`performance.now()` with `Date.now()` as an adjustable-clock diagnostic.

The original pipeline stopped with exit 1. The original failure and all samples
remain immutable. Plan section 6 permits one additional five-sample cohort for
ambiguous results; the auditor supports a fixed confirmation preserving the
four-case workload, two AB/BA sessions and per-case instance scope. Additional
and combined ten-value medians must both pass the unchanged gate. The combined
values would come from separate instances, not ten consecutive calls in one
instance. No further confirmation is permitted if that remains inconclusive.

A separate fixed 16-call memory observation may test convergence over a longer
window, with the final eight values equal, two affected inputs, both versions
and two sessions (128 jobs). This cannot rewrite the original six-call failure
or prove unlimited repetition/scratch deallocation. Capacity and ownership
evidence remain separately required. Neither supplemental experiment had run
when this original-group audit was completed.
