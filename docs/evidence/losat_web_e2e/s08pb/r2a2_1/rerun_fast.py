#!/usr/bin/env python3
"""Re-run the round-2 (a) second-pass rows that use -task blastp-fast (chaining) and every
row that differed, with NCBI blastp 2.17.0 and a LOSAT binary; compare stdout, stderr, exit."""
import concurrent.futures as cf
import os
import subprocess
import sys
from pathlib import Path

A = Path.home() / '.cache/losat-web-gui-target/s08pb-audit/a2'
F = Path.home() / '.cache/losat-web-gui-target/s08pb-audit/src2/LOSAT/tests/fasta/outfmt0'
NCBI = '/home/kawato/micromamba/bin/blastp'
LOSAT = sys.argv[1]
OUT = sys.argv[2]
ENV = {'PATH': '/usr/bin:/bin', 'HOME': os.environ['HOME']}

rows = []
for name in [] if os.environ.get('G2_ONLY') else ['w1_11.tsv', 'w1_12f.tsv', 'w2_21.tsv', 'w2_22.tsv', 'w3_31.tsv', 'w3_32.tsv']:
    for line in open(A / name):
        c = line.rstrip('\n').split('\t')
        if len(c) < 6:
            continue
        status, opts, q, s = c[0], c[3], c[4], c[5]
        if 'blastp-fast' in opts or not status.startswith(('SAME', 'LOSAT-REJECTS')):
            rows.append((name, status, opts.split(), q, s))
for line in open(A / 'g2.tsv'):
    c = line.rstrip('\n').split('\t')
    if len(c) < 7:
        continue
    status, files, opts = c[0], c[5].split(), c[6]
    if 'blastp-fast' in opts or not status.startswith('SAME'):
        rows.append(('g2.tsv', status, opts.split(), str(F / files[0]), str(F / files[1])))


def run(row):
    name, status, opts, q, s = row
    argv = ['-query', q, '-subject', s] + opts
    n = subprocess.run([NCBI] + argv, capture_output=True, env=ENV, timeout=1200)
    lo = subprocess.run([LOSAT, 'blastp'] + argv, capture_output=True, env=ENV, timeout=1200)
    if (n.returncode, n.stdout, n.stderr) == (lo.returncode, lo.stdout, lo.stderr):
        cls = 'SAME'
    elif lo.returncode == 1 and b"not supported by LOSAT" in lo.stderr:
        cls = 'LOSAT-REJECTS'
    elif n.returncode < 0 or n.returncode == 139:
        cls = 'NCBI-CRASH'
    elif n.stdout == lo.stdout and n.returncode == lo.returncode:
        cls = 'SAME-STDOUT'
    else:
        cls = 'DIFF'
    return f"{cls}\t{name}\t{status}\tn={n.returncode} l={lo.returncode}\tnb={len(n.stdout)} lb={len(lo.stdout)}\t{' '.join(opts)}\t{q}\t{s}"


with open(OUT, 'w') as out, cf.ThreadPoolExecutor(int(os.environ.get('JOBS', '4'))) as ex:
    counts = {}
    for line in ex.map(run, rows):
        out.write(line + '\n')
        cls = line.split('\t', 1)[0]
        counts[cls] = counts.get(cls, 0) + 1
    out.write(f"# rows={len(rows)} {counts}\n")
print(len(rows), counts)
