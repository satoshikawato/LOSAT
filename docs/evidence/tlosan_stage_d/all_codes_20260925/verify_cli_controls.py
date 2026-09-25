#!/usr/bin/env python3
"""Comparison-only pinned NCBI CLI code-1 controls for every code fixture."""
import hashlib
import os
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
saved = root / 'run_20260925'
out = Path(sys.argv[1])
out.mkdir()
cli, kappa, d_trace = map(Path, sys.argv[2:])
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
assert sha(cli) == 'e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0'
manifest = [f'CLI SHA256: {sha(cli)}', f'Kappa probe SHA256: {sha(kappa)}',
            f'D probe SHA256: {sha(d_trace)}', 'All 27 API code-1 controls equal pinned local -subject CLI code-1.']
for row in (root / 'fixtures.tsv').read_text().splitlines()[1:]:
    code = row.split('\t')[0]
    stem = f'code{code}.control'
    cmd = [str(cli), '-task', 'tblastn', '-query', str(root / 'fixtures' / f'code{code}.faa'),
           '-subject', str(root / 'fixtures' / f'code{code}.fna'), '-db_gencode', '1',
           '-num_threads', '1', '-outfmt', '6']
    for probe, prefix, suffix in ((None, None, 'report.tsv'),
                                  (kappa, b'K_TRACE_', 'kappa.trace'),
                                  (d_trace, b'D_', 'd.trace')):
        env = os.environ.copy()
        if probe: env['LD_PRELOAD'] = str(probe)
        result = subprocess.run(cmd, check=True, capture_output=True, env=env)
        data = result.stdout if prefix is None else b''.join(line for line in result.stderr.splitlines(keepends=True) if line.startswith(prefix))
        expected = saved / f'{stem}.{suffix}'
        assert data == expected.read_bytes(), (code, suffix)
        (out / f'{stem}.{suffix}').write_bytes(data)
    manifest.append(f'{stem}: output, Kappa trace, D trace byte-identical; command: {cmd!r}')
(out / 'manifest.txt').write_text('\n'.join(manifest) + '\n')
(out / 'outputs.sha256').write_text(''.join(f'{sha(path)}  {path.name}\n' for path in sorted(out.iterdir()) if path.is_file() and path.name != 'outputs.sha256'))
