# Independent review of the I2b reuse overlap

Recorded from the read-only `ncbi_parity_auditor` review. No performance run
was started by the auditor. All auditor tools closed at 2026-09-16 02:01:55 UTC.

The saved `/proc/20817/stat` field 22 is 112308 ticks at 100 Hz, establishing
external process start at boot time 1123.08 s. The saved simultaneous
CLOCK_MONOTONIC and CLOCK_BOOTTIME readings differ by 0.000000864 s. The
external argv identifies a gbdraw timing benchmark with 21 samples. Its CPU
affinity is CPU 3, overlapping LOSAT's 0–31 affinity.

| Session/process | Overlap with the external benchmark |
|---|---:|
| AB baseline | 0 s |
| AB I2b | 0 s |
| BA I2b | 38.447467 s |
| BA baseline | 52.645641 s |

The auditor independently recomputed BA's numerical failure: median difference
0.515018104 s (+5.41083%), above the allowed 0.475914165 s. The overlap proves
violation of the exclusive measurement condition; it does not prove that all
of the slowdown comes from that external load.

Preserving the complete original results and numerical FAIL while marking
acceptance invalid is appropriate. Once an exclusive window is confirmed, one
new complete fixed AB+BA set with identical artifacts, conditions, sample counts
and guards is justified in a separate directory. Do not accept only the old AB,
pool old/new samples, replace individual samples, or extend until passing.

The original process reports completed before the attempted stop; the wrapper
ended with its fixed guard failure. `external-overlap.json` records that the
attempted kill found no remaining process. No unrelated process was stopped.
