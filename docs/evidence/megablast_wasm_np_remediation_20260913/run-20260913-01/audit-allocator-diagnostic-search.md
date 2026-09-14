# Independent audit: fixed96 untimed allocation observations

The required read-only `ncbi_parity_auditor` independently verified all96 actual
raw files (128,546,304 bytes) against the four original official oracles. All356
bound files, inputs/argv, Node, host/helpers, original/instrumented source and
artifacts, ABBA ordering, monotonic process non-overlap,768 worker IDs/3072 events,
per-worker ordering and exit0 are consistent. All jobs are untimed; exclusions0.

All16 series reach their observed maximum post-call live usable bytes/blocks
by the second call. Both versions and both sessions have the same final values:

| Case | Live usable bytes | Live blocks | Retained result bytes |
|---|---:|---:|---:|
| short LC BLASTN | 1,419,328 | 40 | 200,989 |
| large BLASTN | 6,334,528 | 40 | 4,238,551 |
| P078 | 1,616,940 | 41 | 385,994 |
| P132 | 2,141,132 | 41 | 530,562 |

Each first-to-second call increases17,600 bytes/16 blocks. Subsequent values
remain constant in13 series. Three temporarily decrease and return to the
same observed maximum: session0 baseline largeN at repeat3 and session0 candidate
largeN at repeat4 each decrease1164 bytes/one block; session1 baseline shortLC
at repeat2 decreases3492 bytes/three blocks. Repeat indexes are zero based.
No later call exceeds the second call's observed live maximum in this window.
Every before count equals the preceding after count for the same instance.
Accounting errors, failed allocations and null-zero-reallocs remain zero.

Linear memory can increase while post-call live counts remain constant. The
largest linear memory is143,917,056 bytes; maximum process RSS573,423,616 bytes.
Every instance stays within its1GiB configured linear-memory limit and every
process within4GiB RSS. Retained result bytes are included in live usage and must
not be added again. Instrumentation can change scheduling and allocation layout.

There are19 retained adjustable-clock disagreements, maximum+0.895623897 seconds.
Outer monotonic/BOOTTIME disagreement is at most5.123 microseconds. Individual
jobs have no BOOTTIME observation. No records are excluded or corrected.

This is a finite six-call observation on instrumented artifacts. It does not
identify every retained block by type or prove an unlimited repetition bound.
The original reuse failure (one timing/two plateau checks) and fixed16-call
failure (six series) remain unchanged. This diagnostic audit is not acceptance
of a replacement memory gate. The independently accepted additional timing
cohort remains separate from the unresolved memory disposition.
