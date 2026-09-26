#!/usr/bin/env python3
"""Run the four independent selected-subject-code TBLASTX audit cases."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'LOSAT/tests'))
import audit_tblastx_v010 as audit  # noqa: E402

CANDIDATE_SHA = '6c50c85d4f14379e3d71952c7cc7c9c3f95976e4bc09921ab7cc421390500b40'
EMPTY_SHA = hashlib.sha256(b'').hexdigest()


# NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
# const string kArgQuery("query"), kArgSubject("subject");
# const string kArgDbGeneticCode("db_gencode"), kArgNumThreads("num_threads");
# Reuse the fixed auditor's command builder and byte/HSP classifier.
def check(case: audit.Case, output_dir: Path, oracle: Path, binary: Path,
          canonical: dict[str, dict[str, str]]) -> dict[str, object]:
    directory = output_dir / case.case_id
    directory.mkdir()
    ncbi = directory / 'ncbi.tsv'
    losat = directory / 'losat.tsv'
    ncbi_command, losat_command = audit.build_commands(case, oracle, binary, ncbi, losat)
    ncbi_result = audit.run_command(ncbi_command, directory / 'ncbi.stderr')
    losat_result = audit.run_command(losat_command, directory / 'losat.stderr', native=True)
    classification, metrics = audit.classify_outputs(ncbi, losat)
    expected = canonical[case.case_id]['losat_sha256']
    record = {
        'case_id': case.case_id, 'contract': case.contract,
        'ncbi_command': ncbi_command, 'losat_command': losat_command,
        'query_sha256': audit.sha256_path(case.query),
        'subject_sha256': audit.sha256_path(case.subject),
        'ncbi_exit': ncbi_result.returncode, 'losat_exit': losat_result.returncode,
        'ncbi_stderr_sha256': hashlib.sha256(ncbi_result.stderr.encode()).hexdigest(),
        'losat_stderr_sha256': hashlib.sha256(losat_result.stderr.encode()).hexdigest(),
        'classification': classification, 'ncbi_sha256': audit.sha256_path(ncbi),
        'losat_sha256': audit.sha256_path(losat), 'frozen_losat_sha256': expected,
        'output_ids_valid': audit.output_ids_are_valid(case, losat), **metrics,
    }
    record['pass'] = (
        ncbi_result.returncode == losat_result.returncode == 0
        and record['ncbi_stderr_sha256'] == record['losat_stderr_sha256'] == EMPTY_SHA
        and classification == 'HSP_SET_DIFF'
        and record['losat_sha256'] == expected
        and record['output_ids_valid']
    )
    (directory / 'record.json').write_text(json.dumps(record, indent=2) + '\n')
    print(f"{case.case_id} {'PASS' if record['pass'] else 'FAIL'}", flush=True)
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(exist_ok=False)
    oracle = audit.DEFAULT_ORACLE
    binary = ROOT / 'LOSAT/target/release/LOSAT'
    if audit.sha256_path(oracle) != audit.EXPECTED_ORACLE_SHA256:
        raise RuntimeError('NCBI oracle SHA mismatch')
    if audit.sha256_path(binary) != CANDIDATE_SHA:
        raise RuntimeError('LOSAT candidate SHA mismatch')
    manifest = audit.load_manifest(ROOT / 'LOSAT/tests/tblastx_v010_parity_manifest.tsv', ROOT)
    chosen = [case for case in manifest if case.case_id.startswith(('d02_', 'd03_', 'd04_', 'd05_'))]
    with (ROOT / 'LOSAT/tests/platform_native_v010_canonical.tsv').open() as stream:
        canonical = {row['case_id']: row for row in csv.DictReader(
            (line for line in stream if not line.startswith('#')), delimiter='\t'
        ) if row['program'] == 'tblastx'}
    with ThreadPoolExecutor(max_workers=4) as pool:
        records = list(pool.map(lambda case: check(case, output_dir, oracle, binary, canonical), chosen))
    (output_dir / 'summary.json').write_text(json.dumps(records, indent=2) + '\n')
    if not all(row['pass'] for row in records):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
