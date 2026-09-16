# I4: keep the indexed driver in the original caller

I3 is rejected: the fixed serial Mje screen is 13.540885 → 16.847026 s
(+24.416% body); all eight outputs match the oracle. Sampled pre/during/post
process observations showed no competing build/test/benchmark. There is no
sample extension. Actual timed-artifact group WAT differs despite inlining:
baseline 16,902 lines/131 locals, I3 16,961 lines/133 locals. Counts alone do
not explain the time delta; inlining did not restore the original code shape.

I4 restores the Native/plain-WASI indexed driver at its original source call
site. It preserves the original current_idx/helper bindings and -= decrement.
A module-scope macro owns the selection/trace statements once, with every
external local passed explicitly. Threaded WASI retains a real helper call
and the prefix driver. The helper is also available in Native unit tests;
those tests do not replace validation of the real indexed caller.

Before performance: check token/trace preservation; 802-case differential
must extract the real caller driver for Native/plain WASI and the helper for
threaded WASI; inspect serial emitted group code against baseline; run the
focused raw gate. First fixed performance screen remains Mje serial n1,
ordinary Node, warmup one plus three AB/BA/AB pairs, original time/RSS guards.
Stop on failure. No I3/I2b measurement is transferred. Only a passing screen
proceeds to remaining cold n8/n1, Native, reuse, browser and final gates.

NCBI authority: c++/src/algo/blast/core/link_hsps.c:812-895. Initial best setup,
scan order, strict ties, floating-point evaluation and caller updates remain
at their original semantic boundaries. No new runtime option or heuristic.

## Prospective desktop-environment protocol (2026-09-16, before new timings)

External editor Git indexing has continued/restarted for more than 30 minutes,
using roughly 0.1–0.2 CPU core in observations. The plan forbids simultaneous
builds, tests, benchmarks, browser QA and compression; these are all complete.
A new fixed study may run under recorded normal-desktop background activity.
It is not described as exclusive or free of background interference.

`study_guard.py` observes processes every second. The prospective validity bounds
are mean observed background CPU <= 0.5 core and every observed interval <= 1.0
core, with no observed external build/test/benchmark/browser/compression process.
Timing, raw-output, thread, RSS and linear-memory guards are unchanged. Every
sample is retained; no fixed set is extended. The previous diagnostic-only
serial set stays ineligible and is never pooled with this new study.

The observation cannot prove absence of subsecond/exited processes, I/O, cache
or CPU-frequency effects. A passing environment guard is a bounded desktop
comparison, not proof of hardware isolation. Noise or ambiguous results remain
inconclusive. Report the background scope alongside any accepted speed claim.
First condition remains serial Mje; only a passing condition proceeds to the
primary threaded n8/n1 conditions and remaining controls/reuse/browser gates.

## User steering (2026-09-16)

The user reiterated that the 1.20 Native ratio is not important: a real speed
improvement is sufficient motivation. The ratio is neither an adoption nor a
completion condition and will not lead the report. Preserve exact output and
required non-regression checks; judge I4 by its improvement over the baseline.
