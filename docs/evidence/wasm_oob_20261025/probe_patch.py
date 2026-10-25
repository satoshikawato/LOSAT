#!/usr/bin/env python3
"""Patch LOSAT/tests/wasi_thread_host.js (on the runner only) for the n16 OOB repro.

- After a trap in a thread worker's wasi_thread_start, print one PROBE line: the memory length
  the worker's isolate now sees, and the start arguments (stack, tls_base, start_func, arg).
  Nothing is read before the trap, so the patch changes no behaviour of a run that does not trap.
- OOB_REFRESH=1: call memory.grow(0) on the worker's memory right before wasi_thread_start
  (V8 refreshes the calling isolate's bounds of a shared memory on grow, also grow(0)).
"""
import sys
from pathlib import Path

p = Path(sys.argv[1]); s = p.read_text()
old = "        instance.exports.wasi_thread_start(tid, startArg);\n"
new = """        if (process.env.OOB_REFRESH === "1") memory.grow(0);
        try {
          instance.exports.wasi_thread_start(tid, startArg);
        } catch (error) {
          const buffer = memory.buffer, dv = new DataView(buffer), a = startArg >>> 0;
          const stack = dv.getUint32(a, true), tls = dv.getUint32(a + 4, true);
          process.stderr.write(`PROBE ${JSON.stringify({ tid, startArg: a, length: buffer.byteLength,
            stack, tls, tidStore: tls + 1384, inBounds: tls + 1388 <= buffer.byteLength,
            startFunc: dv.getUint32(a + 8, true), arg: dv.getUint32(a + 12, true),
            standby: standby ?? null, error: String(error) })}\\n`);
          throw error;
        }
"""
assert s.count(old) == 1, "anchor not found"
p.write_text(s.replace(old, new))
print(f"patched {p}")
