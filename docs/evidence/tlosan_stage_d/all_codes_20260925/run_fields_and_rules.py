#!/usr/bin/env python3
"""Save custom numeric report fields and raw adjustment rules for all codes."""
import hashlib
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
STAGE_D = HERE.parent
OUT, ORACLE, RULE, D_TRACE, CLI = map(Path, sys.argv[1:])
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
manifest = [
    'NCBI source commit: 598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4',
    f'NCBI CLI SHA256: {sha(CLI)}',
    f'Derived fields oracle source SHA256: {sha(OUT / "fields_oracle.cpp")}',
    f'NCBI rule probe SHA256: {sha(STAGE_D / "ncbi_kappa_rule_trace.c")}',
    f'NCBI D probe SHA256: {sha(STAGE_D / "ncbi_d_call_trace.c")}',
    f'Oracle binary SHA256: {sha(ORACLE)}',
    'Custom formatter changes only CBlastFormat custom_output_format; NCBI search state remains the pinned CLI-calibrated API path.',
]
fields = 'qseqid sseqid score bitscore evalue nident positive length mismatch gaps gapopen qstart qend sstart send sframe'
def run(cmd, probe=None, prefix=None):
    env = os.environ.copy()
    if probe: env['LD_PRELOAD'] = str(probe)
    result = subprocess.run(cmd, capture_output=True, check=True, env=env)
    selected = b'' if prefix is None else b''.join(line for line in result.stderr.splitlines(keepends=True) if line.startswith(prefix))
    ordinary = result.stderr if prefix is None else b''.join(line for line in result.stderr.splitlines(keepends=True) if not line.startswith(prefix))
    return result.stdout, ordinary, selected
for row in (HERE / 'fixtures.tsv').read_text().splitlines()[1:]:
    code = int(row.split('\t')[0])
    query = HERE / 'fixtures' / f'code{code}.faa'
    subject = HERE / 'fixtures' / f'code{code}.fna'
    manifest += [f'code{code} query SHA256: {sha(query)}', f'code{code} subject SHA256: {sha(subject)}']
    for selected, label in ((1, 'control'), (code, 'selected')):
        stem = f'code{code}.{label}'
        cmd = [str(ORACLE), str(query), str(subject), str(selected), '6']
        plain, stderr, _ = run(cmd)
        probed, probed_stderr, rules = run(cmd, RULE, b'K_RULE_NEW\t')
        assert (plain, stderr) == (probed, probed_stderr) and rules
        (OUT / f'{stem}.fields.tsv').write_bytes(plain)
        (OUT / f'{stem}.rules.trace').write_bytes(rules)
        manifest.append(f'{stem}: command {(["$NCBI_FIELDS_ORACLE"] + cmd[1:])!r}; rule probe leaves stdout/stderr unchanged')
        if selected == 1:
            cli_cmd = [str(CLI), '-task', 'tblastn', '-query', str(query), '-subject', str(subject),
                       '-db_gencode', '1', '-num_threads', '1', '-outfmt', '6 ' + fields]
            cli_out, cli_err, _ = run(cli_cmd)
            assert cli_out == plain and cli_err == stderr, (code, 'custom fields CLI calibration')
            # Search input/call state must still match the pinned code-1 CLI.
            _, api_err, api_d = run(cmd, D_TRACE, b'D_')
            _, cli_d_err, cli_d = run(cli_cmd, D_TRACE, b'D_')
            assert api_d == cli_d and api_err == cli_d_err, (code, 'D call state')
            manifest.append(f'{stem}: custom report and D trace byte-identical to code-1 local CLI')
(OUT / 'fields_manifest.txt').write_text('\n'.join(manifest) + '\n')
(OUT / 'fields_outputs.sha256').write_text(''.join(
    f'{sha(path)}  {path.name}\n' for path in sorted(OUT.iterdir())
    if path.is_file() and path.name in {'fields_oracle.cpp', 'fields_manifest.txt'} or
    path.is_file() and path.suffix in {'.tsv', '.trace'} and path.name.startswith('code') and
    (path.name.endswith('.fields.tsv') or path.name.endswith('.rules.trace'))
))
