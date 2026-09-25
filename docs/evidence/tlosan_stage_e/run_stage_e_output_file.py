#!/usr/bin/env python3
"""Compare public TBLASTN -out files against pinned NCBI for 0/6/7."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from run_stage_e_cli_matrix import ROOT, NCBI, LOSAT, SOURCE, PIN, NCBI_SHA, digest


def main() -> int:
    out = Path(sys.argv[1]).resolve()
    out.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(['git', '-C', str(SOURCE), 'rev-parse', 'HEAD']).decode().strip() == PIN
    assert digest(NCBI.read_bytes()) == NCBI_SHA
    query = ROOT / 'docs/evidence/tlosan_stage_c/multi_hsp_20260924/query.faa'
    subject = ROOT / 'docs/evidence/tlosan_stage_c/multi_hsp_20260924/subjects.fna'
    rows = []
    for fmt in ('0', '6', '7'):
        expected_file = out / f'fmt{fmt}.ncbi.out'
        actual_file = out / f'fmt{fmt}.losat.out'
        common = ['-query', str(query), '-subject', str(subject), '-num_threads', '1',
                  '-outfmt', fmt]
        expected_cmd = [str(NCBI), '-task', 'tblastn', *common, '-out', str(expected_file)]
        actual_cmd = [str(LOSAT), 'tblastn', *common, '-out', str(actual_file)]
        expected = subprocess.run(expected_cmd, cwd=ROOT, capture_output=True)
        actual = subprocess.run(actual_cmd, cwd=ROOT, capture_output=True)
        expected_bytes = expected_file.read_bytes() if expected_file.exists() else b''
        actual_bytes = actual_file.read_bytes() if actual_file.exists() else b''
        passed = expected.returncode == actual.returncode == 0 and not expected.stdout and not actual.stdout and expected_bytes == actual_bytes
        rows.append({'outfmt': fmt, 'ncbi_command': expected_cmd, 'losat_command': actual_cmd,
                     'query_sha256': digest(query.read_bytes()), 'subject_sha256': digest(subject.read_bytes()),
                     'ncbi_exit': expected.returncode, 'losat_exit': actual.returncode,
                     'ncbi_file_sha256': digest(expected_bytes), 'losat_file_sha256': digest(actual_bytes),
                     'ncbi_stdout_sha256': digest(expected.stdout), 'losat_stdout_sha256': digest(actual.stdout),
                     'bytes': len(expected_bytes), 'equal': passed})
        print(f'{"PASS" if passed else "FAIL"} -out fmt{fmt}: {len(expected_bytes)} bytes', flush=True)
    (out / 'comparison.jsonl').write_text(''.join(json.dumps(row, sort_keys=True)+'\n' for row in rows))
    return 0 if all(row['equal'] for row in rows) else 1


if __name__ == '__main__':
    raise SystemExit(main())
