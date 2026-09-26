#!/usr/bin/env python3
"""Replay three v0.1.0 TBLASTX profile probes on the v0.2.0 serial module."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[3]
MODULE = ROOT / 'LOSAT/target/serial-command/wasm32-wasip1/release/LOSAT.wasm'
RUNNER = ROOT / 'LOSAT/tests/run_losat_wasi.js'
MANIFEST = ROOT / 'LOSAT/tests/tblastx_v010_parity_manifest.tsv'
CANONICAL = ROOT / 'LOSAT/tests/platform_native_v010_canonical.tsv'
MODULE_SHA = '912ea816dd5383b423011b02cacc32ade1d3c61cd02a391a4944f82b4b84eee5'
CASES = ('p03_mela_pemojnva', 'p12_lc738874_lc738875_default', 'd06_ap027131_ap027133_db4')


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    if not str(output_dir).startswith('/tmp/') or output_dir.exists():
        raise RuntimeError('output directory must be a new path under /tmp')
    if sha256(MODULE) != MODULE_SHA:
        raise RuntimeError('module is not the exact v0.2.0 candidate')
    rows = {row['case_id']: row for row in csv.DictReader(MANIFEST.open(), delimiter='\t')}
    frozen = {
        row['case_id']: row for row in csv.DictReader(
            (line for line in CANONICAL.open() if not line.startswith('#')), delimiter='\t'
        ) if row['program'] == 'tblastx'
    }
    output_dir.mkdir(parents=True)
    version = subprocess.run(['node', '--version'], capture_output=True, text=True, check=True).stdout.strip()
    results = []
    # NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
    # const string kArgQuery("query"), kArgSubject("subject");
    # const string kArgQueryGeneticCode("query_gencode"), kArgDbGeneticCode("db_gencode");
    # The host supplies the same local-query/local-subject options as the certified native profile.
    for case_id in CASES:
        row = rows[case_id]
        output = output_dir / f'{case_id}.tsv'
        command = [
            'node', '--no-warnings', '--experimental-wasi-unstable-preview1',
            str(RUNNER), str(MODULE), 'tblastx', '-query', str(ROOT / row['query']),
            '-subject', str(ROOT / row['subject']), '-outfmt', '6', '-num_threads', '1',
        ]
        if row['gencode_args'] == 'explicit':
            command += ['-query_gencode', row['query_gencode'], '-db_gencode', row['db_gencode']]
        command += ['-out', str(output)]
        run = subprocess.run(command, capture_output=True, check=False)
        actual = sha256(output) if output.is_file() else None
        expected = frozen[case_id]['losat_sha256']
        result = {
            'case_id': case_id, 'contract': row['contract'], 'module_sha256': MODULE_SHA,
            'node_version': version, 'command': command,
            'query_sha256': sha256(ROOT / row['query']),
            'subject_sha256': sha256(ROOT / row['subject']),
            'exit': run.returncode, 'stdout_sha256': hashlib.sha256(run.stdout).hexdigest(),
            'stderr_sha256': hashlib.sha256(run.stderr).hexdigest(),
            'output_bytes': output.stat().st_size if output.is_file() else None,
            'output_sha256': actual, 'frozen_native_sha256': expected,
            'equal': run.returncode == 0 and not run.stderr and actual == expected,
        }
        results.append(result)
        (output_dir / 'comparison.json').write_text(json.dumps(results, indent=2) + '\n')
        print(f"{case_id}: {'PASS' if result['equal'] else 'FAIL'}", flush=True)
    return 0 if all(row['equal'] for row in results) else 1


if __name__ == '__main__':
    raise SystemExit(main())
