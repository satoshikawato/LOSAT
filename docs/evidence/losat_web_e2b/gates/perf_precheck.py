#!/usr/bin/env python3
"""Before/after output equality of the S08 V-PERF cases (native, -num_threads 1 and 4)."""
import hashlib, importlib.util, subprocess, sys
from pathlib import Path
W = Path("/mnt/c/Users/genom/GitHub/LOSAT-web-gui")
spec = importlib.util.spec_from_file_location("pc", W / "docs/evidence/losat_web_e2b/perf_cases.py")
pc = importlib.util.module_from_spec(spec); spec.loader.exec_module(pc)
before, after, out = sys.argv[1], sys.argv[2], Path(sys.argv[3])
cases = "tblastx,tblastx-multi,tblastx-many,tblastn,tblastn-fmt0,blastp,blastp-fmt0,blastn-large-fmt0,blastn-large".split(",")
for case in cases:
    argv = pc.measure_perf.FIXTURES[case]
    for threads in (1, 4):
        h = {}
        for side, binary in (("before", before), ("after", after)):
            o = out / f"{case}.t{threads}.{side}.out"
            p = subprocess.run([binary] + argv + ["-num_threads", str(threads), "-out", str(o)], capture_output=True)
            h[side] = (p.returncode, hashlib.sha256(o.read_bytes()).hexdigest()[:16] if o.exists() else "-", len(p.stderr))
        same = h["before"] == h["after"]
        print(f"{case}\tthreads {threads}\t{'same' if same else 'DIFFERENT'}\tbefore {h['before']}\tafter {h['after']}", flush=True)
