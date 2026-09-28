#!/usr/bin/env python3
"""Confirm unsupported public TBLASTN paths fail before producing a report."""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

from run_stage_e_cli_matrix import ROOT, LOSAT

QUERY = ROOT / 'docs/evidence/tlosan_stage_c/multi_hsp_20260924/query.faa'
SUBJECT = ROOT / 'docs/evidence/tlosan_stage_c/multi_hsp_20260924/subjects.fna'
CASES = {
    'database_search': ['-db', 'unported'],
    'fast_task': ['-task', 'tblastn-fast'],
    'checkpoint': ['-in_pssm', 'unported.asn'],
    'remote': ['-remote'],
    'subject_location': ['-subject_loc', '1-30'],
    'ungapped': ['-ungapped', '-comp_based_stats', 'F'],
    'threads': ['-num_threads', '2'],
    'matrix': ['-matrix', 'PAM30'],
    'intron': ['-max_intron_length', '1'],
    'composition': ['-comp_based_stats', '1'],
    'word_size': ['-word_size', '7'],
    'invalid_code': ['-db_gencode', '7'],
    'xml_format': ['-outfmt', '5'],
    'custom_fields': ['-outfmt', '6 qseq sseq'],
}


def sha(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def main() -> int:
    out = Path(sys.argv[1]).resolve()
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    for name, option in CASES.items():
        target = out / f'{name}.out'
        command = [str(LOSAT), 'tblastn', '-query', str(QUERY), '-subject', str(SUBJECT),
                   '-out', str(target), *option]
        result = subprocess.run(command, cwd=ROOT, capture_output=True)
        passed = result.returncode != 0 and not result.stdout and not target.exists()
        row = {'case': name, 'command': command, 'exit': result.returncode,
               'stdout_sha256': sha(result.stdout), 'stderr_sha256': sha(result.stderr),
               'stderr': result.stderr.decode(errors='replace'),
               'output_file_created': target.exists(), 'pass': passed}
        rows.append(row)
        print(f'{"PASS" if passed else "FAIL"} {name}: exit={result.returncode}', flush=True)
    (out / 'comparison.jsonl').write_text(''.join(json.dumps(row, sort_keys=True)+'\n' for row in rows))
    return 0 if all(row['pass'] for row in rows) else 1


if __name__ == '__main__':
    raise SystemExit(main())
