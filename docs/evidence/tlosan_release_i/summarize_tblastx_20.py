#!/usr/bin/env python3
"""Verify the exact-candidate 20-case TBLASTX audit assembled from parallel runs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'LOSAT/tests'))
import audit_tblastx_v010 as audit  # noqa: E402

SOURCE_COMMIT = '005e3d4b6cba6b5808334088fe9595c89efe01f8'
CANDIDATE_SHA = '6c50c85d4f14379e3d71952c7cc7c9c3f95976e4bc09921ab7cc421390500b40'
EMPTY_SHA = hashlib.sha256(b'').hexdigest()


def records(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline='') as stream:
        return {row['case_id']: row for row in csv.DictReader((line for line in stream if not line.startswith('#')), delimiter='\t')}


def source_dir(case_id: str, args: argparse.Namespace) -> tuple[str, Path]:
    if case_id.startswith(('p01_', 'p02_', 'p03_', 'p04_', 'p05_', 'p06_')):
        return 'original_prefix', args.prefix / case_id
    if case_id == 'd01_nz_self_code4':
        return 'deviation_prefix', args.deviation / case_id
    if case_id.startswith(('d02_', 'd03_', 'd04_', 'd05_')):
        return 'direct_parallel', args.direct / case_id
    if case_id.startswith(('p07_', 'p08_', 'p09_', 'p10_', 'p11_')):
        return 'later_prefix', args.later / case_id
    return 'future_controls', args.future / case_id


# NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
# const string kArgQuery("query"), kArgSubject("subject");
# const string kArgQueryGeneticCode("query_gencode"), kArgDbGeneticCode("db_gencode");
# This summary reclassifies raw files made by the fixed oracle and candidate commands.
# It accepts only exact NCBI bytes or the six manifest-designated subject-code rows.

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', type=Path, required=True)
    parser.add_argument('--later', type=Path, required=True)
    parser.add_argument('--future', type=Path, required=True)
    parser.add_argument('--deviation', type=Path, required=True)
    parser.add_argument('--direct', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest = audit.load_manifest(ROOT / 'LOSAT/tests/tblastx_v010_parity_manifest.tsv', ROOT)
    fixture_paths = sorted({str(path.relative_to(ROOT)) for case in manifest for path in (case.query, case.subject)})
    fixture_diff = subprocess.run(
        ['git', 'diff', '--quiet', SOURCE_COMMIT, '--', *fixture_paths], cwd=ROOT, check=False
    )
    if fixture_diff.returncode != 0:
        raise RuntimeError('fixture inputs differ from exact candidate source commit')
    if audit.sha256_path(audit.DEFAULT_ORACLE) != audit.EXPECTED_ORACLE_SHA256:
        raise RuntimeError('current NCBI executable SHA mismatch')
    canonical = records(ROOT / 'LOSAT/tests/platform_native_v010_canonical.tsv')
    future_rows = records(args.future / 'classifications.tsv')
    default_probe = json.loads((args.future / 'default_code_equivalence.json').read_text())
    default_sha = canonical['p12_lc738874_lc738875_default']['losat_sha256']
    default_outputs = [
        args.future / 'p12_lc738874_lc738875_default' / 'ncbi.tsv',
        args.future / 'p12_lc738874_lc738875_default' / 'losat.tsv',
        args.future / 'default_code_explicit.ncbi.tsv',
        args.future / 'default_code_explicit.losat.tsv',
    ]
    default_stderr = [args.future / f'default_code_explicit.{program}.stderr' for program in ('ncbi', 'losat')]
    if (
        default_probe['status'] != 'PASS'
        or default_probe['ncbi_exit'] != 0 or default_probe['losat_exit'] != 0
        or any(audit.sha256_path(path) != default_sha for path in default_outputs)
        or any(audit.sha256_path(path) != EMPTY_SHA for path in default_stderr)
    ):
        raise RuntimeError('implicit/explicit code-1 control failed')
    oracle = json.loads((args.prefix / 'oracle_identity.json').read_text())
    if oracle['sha256'] != audit.EXPECTED_ORACLE_SHA256:
        raise RuntimeError('oracle SHA mismatch')
    if audit.sha256_path(ROOT / 'LOSAT/target/release/LOSAT') != CANDIDATE_SHA:
        raise RuntimeError('candidate binary SHA mismatch')
    original_log = Path(str(args.prefix) + '.log').read_text()
    deviation_log = Path(str(args.deviation) + '.log').read_text()
    later_log = Path(str(args.later) + '.log').read_text()
    rows = []
    for case in manifest:
        source, directory = source_dir(case.case_id, args)
        ncbi = directory / 'ncbi.tsv'
        losat = directory / 'losat.tsv'
        classification, metrics = audit.classify_outputs(ncbi, losat)
        expected = canonical[case.case_id]['losat_sha256']
        stderr_paths = [directory / 'ncbi.stderr', directory / 'losat.stderr']
        repeat_paths = [losat]
        if case.repeatability_required:
            repeat_paths += [directory / 'losat.run2.tsv', directory / 'losat.run3.tsv']
            stderr_paths += [directory / 'losat.run2.stderr', directory / 'losat.run3.stderr']
        stderr_hashes = [audit.sha256_path(path) for path in stderr_paths]
        repeat_hashes = [audit.sha256_path(path) for path in repeat_paths]
        if source == 'direct_parallel':
            reported = json.loads((directory / 'record.json').read_text())
            exit_ok = reported['ncbi_exit'] == reported['losat_exit'] == 0 and reported['pass']
            commands = {'ncbi': reported['ncbi_command'], 'losat': reported['losat_command']}
        elif source == 'future_controls':
            reported = future_rows[case.case_id]
            exit_ok = reported['ncbi_exit'] == reported['losat_exit'] == '0' and reported['contract_status'] == 'PASS'
            commands = json.loads((directory / 'commands.json').read_text())
        else:
            log = {'original_prefix': original_log, 'deviation_prefix': deviation_log, 'later_prefix': later_log}[source]
            expected_line = f'{case.case_id}: {classification} rows='
            exit_ok = any(expected_line in line and 'contract=PASS' in line for line in log.splitlines())
            commands = json.loads((directory / 'commands.json').read_text())
        expected_ncbi, expected_losat = audit.build_commands(
            case, audit.DEFAULT_ORACLE, ROOT / 'LOSAT/target/release/LOSAT', ncbi, losat
        )
        command_ok = commands == {'ncbi': expected_ncbi, 'losat': expected_losat}
        expected_class = 'EXACT_TEXT' if case.contract == audit.PARITY_CONTRACT else 'HSP_SET_DIFF'
        accepted = (
            exit_ok and command_ok and classification == expected_class
            and metrics['losat_sha256'] == expected
            and all(value == EMPTY_SHA for value in stderr_hashes)
            and len(repeat_hashes) == (3 if case.repeatability_required else 1)
            and len(set(repeat_hashes)) == 1
            and audit.output_ids_are_valid(case, losat)
        )
        row = {
            'case_id': case.case_id, 'source': source, 'contract': case.contract,
            'classification': classification, 'expected_classification': expected_class,
            'candidate_binary_sha256': CANDIDATE_SHA, 'oracle_sha256': oracle['sha256'],
            'commands': commands, 'commands_match_manifest': command_ok, 'query_sha256': audit.sha256_path(case.query),
            'subject_sha256': audit.sha256_path(case.subject),
            'ncbi_sha256': metrics['ncbi_sha256'], 'losat_sha256': metrics['losat_sha256'],
            'frozen_losat_sha256': expected, 'ncbi_rows': metrics['ncbi_rows'],
            'losat_rows': metrics['losat_rows'], 'stderr_sha256': stderr_hashes,
            'repeatability_sha256': repeat_hashes, 'output_ids_valid': audit.output_ids_are_valid(case, losat),
            'source_reports_zero_exit_and_pass': exit_ok, 'accepted': accepted,
        }
        rows.append(row)
        print(f"{case.case_id}: {'PASS' if accepted else 'FAIL'}", flush=True)
    counts = {
        'unique_cases': len({row['case_id'] for row in rows}),
        'exact_ncbi': sum(row['classification'] == 'EXACT_TEXT' and row['accepted'] for row in rows),
        'approved_subject_code': sum(row['classification'] == 'HSP_SET_DIFF' and row['accepted'] for row in rows),
    }
    status = 'PASS' if len(rows) == 20 and counts == {
        'unique_cases': 20, 'exact_ncbi': 14, 'approved_subject_code': 6,
    } and all(row['accepted'] for row in rows) else 'FAIL'
    args.output.write_text(json.dumps({
        'candidate_source_commit': SOURCE_COMMIT,
        'candidate_binary_sha256': CANDIDATE_SHA, 'oracle': oracle,
        'default_code_equivalence': default_probe, 'counts': counts,
        'status': status, 'rows': rows,
    }, indent=2) + '\n')
    if status != 'PASS':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
