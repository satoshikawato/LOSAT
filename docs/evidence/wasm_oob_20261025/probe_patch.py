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
          // Print the error first; then what this isolate sees as the memory length before and
          // after memory.grow(0), and the start arguments where they are readable.
          process.stderr.write(`PROBE-ERR ${String(error)}\n`);
          const a = startArg >>> 0, read = buffer => {
            const length = buffer.byteLength, info = { length };
            if (a + 16 <= length) {
              const dv = new DataView(buffer);
              Object.assign(info, { stack: dv.getUint32(a, true), tls: dv.getUint32(a + 4, true) });
              info.tidStoreInBounds = info.tls + 1388 <= length;
            }
            return info;
          };
          const before = read(memory.buffer); memory.grow(0); const after = read(memory.buffer);
          process.stderr.write(`PROBE ${JSON.stringify({ tid, startArg: a, standby: standby ?? null, before, after })}\n`);
          throw error;
        }
"""
assert s.count(old) == 1, "anchor not found"
p.write_text(s.replace(old, new))
print(f"patched {p}")
