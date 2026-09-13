# Independent timing audit — preliminary findings

Auditor: `ncbi_parity_auditor`, read-only. No benchmark/build/edit performed by the auditor. Audit cutoff: 2026-09-13 08:12:51 UTC while the controls run was still producing samples.

Verdict at cutoff: elapsed-time adoption claims remain provisional. Raw-output and worker-lifecycle findings are independent of this timing issue.

The frozen runner in `p1-controls/runner-source/wasm_performance.py` samples monotonic elapsed before reading GNU usage and hashing output, but samples second-resolution `ended_utc` afterward. UTC span therefore includes bookkeeping and quantization. It is not the same boundary. However, this cannot explain a 300-second watchdog expiring within 282.01989235001383 recorded monotonic seconds. That Pemo baseline row has no output or GNU usage and remains TIMEOUT; it must never become a 300-second success.

At the cutoff, 19 of 47 completed records had UTC span minus monotonic greater than max(1 second, 5% of monotonic): 10/32 cold PASS, 5/9 diagnostic PASS, 1/1 diagnostic TIMEOUT, 3/5 oracle PASS. Eighteen of 46 successful records exceeded this exploratory cutoff; it does not by itself prove a clock fault because end timestamps include bookkeeping.

Current available system evidence: WSL2 5.15.167.4-microsoft-standard-WSL2; `tsc` clocksource, `hyperv_clocksource_tsc_page` available; time-namespace offsets zero. Sequential monotonic and boottime readings differed by 867 ns. Suspend statistics were absent, and available filtered kernel/timesyncd logs provided no historical explanation. This cannot rule out host VM pauses or earlier clock events.

Node diagnostic last-event timestamps agreed with parent monotonic within about 30–80 ms. First-spawn to last-exit intervals were AP027280 67.141→23.989 s, MjeNMV 37.359→14.183 s, EDL933 megablast 2.214→2.843 s, WSSV BLASTP 1.062→1.812 s. These have narrower boundaries and use `performance.now()`, so they corroborate monotonic computation observations but not independent realtime calibration. GNU user+system CPU totals corroborate the direction of the changes; CPU is not wall time or worker utilization. The old `%U/%S/%M` files do not contain elapsed time.

Future runner records retain immediate completion timestamps, realtime and boottime durations, GNU elapsed, watchdog start/fire clock samples, and separate bookkeeping completion time. Reuse samples also compare Date.now() at the same runtime/close endpoints. Future clock disagreement excludes timing separately from parity. Old records remain immutable. Idle sleep/watchdog controls will be run after the current performance subprocess ends.

## Follow-up: mechanism and measurement policy

The idle `clock-probe/` reproduces the difference inside a one-second child sleep: realtime 3.308256 s, monotonic 1.000132 s, boottime 1.000134 s. Other short controls agree. This rules out evidence hashing as the cause of that discrepancy; it does not identify the host event or externally calibrate monotonic.

The independent auditor traced the installed CPython 3.13.3 executable, SHA-256 `56b34efef3ef5e625540e7b8f5df4b706e6f1d2d0ea5bad28b6e7435c55456db`. The actual call path is Timer → Event.wait → Condition.wait → lock_PyThread_acquire_lock → _PyMutex_LockTimed → _PyParkingLot_Park → _PySemaphore_Wait. This build has `HAVE_SEM_CLOCKWAIT=0` and `HAVE_SEM_TIMEDWAIT=1`. Its timed semaphore path calls PyTime_TimeRaw and sem_timedwait, and returns on timeout without validating the mutex's monotonic deadline. The local disassembly agrees with official [CPython parking_lot.c](https://github.com/python/cpython/blob/v3.13.3/Python/parking_lot.c) and [lock.c](https://github.com/python/cpython/blob/v3.13.3/Python/lock.c). Thus forward realtime changes can cause an early watchdog. Attribution of the historical Pemo timeout to a particular clock change remains an inference.

The declared elapsed metric remains CLOCK_MONOTONIC, consistently for baseline and candidate. Realtime/GNU elapsed offsets remain diagnostics; they do not independently invalidate monotonic results. Monotonic/boottime disagreement still rejects new timing, as do process, parity and lifecycle failures. The auditor accepts this bounded policy; it makes no physical-clock calibration claim. Old successful observations retain their original boundary, and old timeouts remain excluded.

The watchdog callback now rechecks its monotonic deadline and waits the remaining duration through its cancellation event if awakened early. This prevents premature termination from a forward realtime step and preserves cancellation. The initial realtime-sensitive wait may still return late under a backward realtime step, so this is not a strict upper-bound guarantee in such an environment. Tests inject an early timer wake to check both successful completion/cancellation and eventual timeout.


Verified source line ranges: CPython v3.13.3 `_threadmodule.c:777–791`, `lock.c:65–70,108–109,122–124`, `parking_lot.c:124–158,303–323`, and `timemodule.c:2259–2263`. Local `threading.py`: Condition.wait 327–373, Event.wait 641–660, Timer.run 1339–1343. GNU time's installed executable calls gettimeofday before fork and after completion (disassembly 0x2657/0x26c8).

Additional audited file identities:

| File | SHA-256 |
|---|---|
| `/home/kawato/micromamba/include/python3.13/pyconfig.h` | `0a5844461522c24b4d4e7322ecef319e6e4b30b7be1f7ef45aacec0920997b05` |
| `/home/kawato/micromamba/lib/python3.13/threading.py` | `dbc5f71103224bac1b650da8c3c306a44a1963315bcdbd9bee8e3421f42fff00` |
| `/usr/bin/time` | `3b11dec50514a8473e9f6efa7a34d584d0657538c09988f61b72d38ad4991a10` |

## Final consumer admission check

The independent audit found that the legacy `run()` and `run_warm()` consumers recorded clock disagreement without excluding it from medians. They now reject CLOCK_MONOTONIC/CLOCK_BOOTTIME disagreement from timing admission. The raw audit phase retains its independent byte gate; adjustable realtime disagreement alone remains admissible. A synthetic consumer regression exercises cold/warm rejection, realtime-only acceptance and raw-audit preservation. The detailed threading benchmark already had this check, so previously audited adoption measurements are unaffected.
