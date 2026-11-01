#!/usr/bin/env python3
"""Line WO: copies a tree's Node WASI host and adds a report after a worker's trap.

Nothing changes before the trap: the host's thread start only gets a try/catch. In the catch it
prints this isolate's view of the memory, the true size (after memory.grow(0)) and the value of
the i32 at WO_PTR (the address the trapping rmw.add used, when the trap is the Collector clone).

usage: mk_host.py <tree>/LOSAT/tests <dest dir>
"""
import pathlib
import shutil
import sys

src = pathlib.Path(sys.argv[1])
dest = pathlib.Path(sys.argv[2])
dest.mkdir(parents=True, exist_ok=True)
for path in src.glob("*.js"):
    if path.name.startswith("wasi_") or path.name == "run_losat_wasi_threads.js":
        shutil.copy(path, dest / path.name)

START = "instance.exports.wasi_thread_start(tid, startArg);"
REPORT = """try { instance.exports.wasi_thread_start(tid, startArg); } catch (trap) {
          try {
            const view = memory.buffer.byteLength;
            memory.grow(0);
            const size = memory.buffer.byteLength;
            const at = Number(process.env.WO_PTR || 0);
            const ptr = at ? new DataView(memory.buffer).getUint32(at, true) : 0;
            require("fs").writeSync(2, `[wo] trap tid=${tid} message=${trap.message} view_after=${view} true_size=${size} ptr=${ptr} ptr_in_true=${ptr + 4 <= size}\\n`);
          } catch (err) { require("fs").writeSync(2, `[wo] report failed: ${err.message}\\n`); }
          throw trap;
        }"""
patched = 0
for path in dest.glob("wasi_thread_*.js"):
    text = path.read_text()
    if START in text:
        path.write_text(text.replace(START, REPORT))
        patched += text.count(START)
if patched != 1:
    sys.exit(f"expected one thread start, patched {patched}")
print(f"host copy with the trap report in {dest}")
