#!/usr/bin/env python3
"""Run NCBI blastp for the capture cases whose LOSAT output changed since SD and compare bytes."""
import csv, json, subprocess, sys, tempfile
from pathlib import Path
REPO = Path('/mnt/c/Users/genom/GitHub/LOSAT-web-gui')
CAP = Path.home() / '.cache/losat-web-gui-target/capture-s08p'
NCBI = '/home/kawato/micromamba/bin/blastp'
cases = sys.argv[1:]
rows = {(r['program'], r['case_id']): r for r in csv.DictReader(open(CAP / 'hashes.tsv'), delimiter='\t')}
for case_id in cases:
    row = rows[('blastp', case_id)]
    command = json.loads(row['command'])
    assert command[0] == 'blastp', command
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / 'n.out'
        argv = [NCBI] + [str(out) if a == '{OUT}' else a for a in command[1:]]
        run = subprocess.run(argv, cwd=REPO, capture_output=True)
        ncbi = out.read_bytes() if '{OUT}' in command else run.stdout
    losat = (CAP / 'outputs/blastp' / f'{case_id}.out').read_bytes()
    print(f"{case_id}\tncbi_exit={run.returncode}\tlosat_exit={row['exit']}\tstdout={'same' if ncbi == losat else 'DIFF'}\tbytes={len(losat)}\tncbi_stderr={run.stderr.decode(errors='replace').strip()[:160]!r}")
