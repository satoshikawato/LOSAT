#!/usr/bin/env python3
"""perf S-I session 7: time one command with the monotonic clock.

usage: s7_time.py OUT_FILE COMMAND...
Writes "wall_s user_s sys_s max_rss_kb" to OUT_FILE and exits with the command's status. The wall time is
time.perf_counter() (CLOCK_MONOTONIC): this WSL host's wall clock steps forward 1.4-2.4 s about every 32 s
(memory wsl-wall-clock-steps), which /usr/bin/time's elapsed time includes.
"""
import os
import sys
import time

out, cmd = sys.argv[1], sys.argv[2:]
start = time.perf_counter()
pid = os.fork()
if pid == 0:
    try:
        os.execvp(cmd[0], cmd)
    finally:
        os._exit(127)
_, status, usage = os.wait4(pid, 0)
wall = time.perf_counter() - start
with open(out, "w") as f:
    f.write(f"{wall:.3f} {usage.ru_utime:.2f} {usage.ru_stime:.2f} {usage.ru_maxrss}\n")
sys.exit(os.waitstatus_to_exitcode(status) if os.waitstatus_to_exitcode(status) >= 0 else 128 - os.waitstatus_to_exitcode(status))
