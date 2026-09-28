#!/usr/bin/env python3
"""Fail-closed exact official-Windows fingerprint lookup for authority v3.

NCBI c++/src/objtools/align_format/tabular.cpp:1100-1108:
// void CBlastTabularInfo::Print() {
//     ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
//         x_PrintField(*iter);
//     }
//     m_Ostream << "\n";
// }
The output byte fingerprint is decisive; no newline conversion is accepted.
"""
import copy
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
NEW_FIELDS = ('program', 'argv', 'cwd', 'effective_environment', 'inputs',
              'exit', 'loaded_module_profile_sha256', 'loaded_modules',
              'raw_stdout', 'raw_stderr')
LEGACY_FIELDS = ('code', 'outfmt', 'threads', 'argv', 'cwd', 'environment',
                 'inputs', 'official_executable_sha256',
                 'declared_runtime_inventory_profile_sha256',
                 'raw_stdout_sha256', 'raw_stderr_sha256')


def lookup(registry, family, observed, platform='windows-x86_64', runtime=None):
    if platform != 'windows-x86_64' or family not in ('m1_n1', 'm1_n4', 'm2_m8', 'pr5_oracle'):
        raise ValueError('UNKNOWN_PLATFORM_OR_FAMILY')
    if family == 'm1_n1':
        # Retain the previously approved v2 declared-inventory contract for
        # these 78 rows. New 27-module debug closure is not inferred backwards.
        matches = [row for row in registry['m1_n1']
                   if row['code'] == observed.get('code') and row['outfmt'] == observed.get('outfmt')
                   and row['threads'] == observed.get('threads')]
        if len(matches) != 1:
            raise ValueError('UNKNOWN_LEGACY_CASE')
        approved = matches[0]
        if any(approved.get(field) != observed.get(field) for field in LEGACY_FIELDS):
            raise ValueError('UNKNOWN_LEGACY_CONTEXT_OR_FINGERPRINT')
        profile = registry['legacy_m1_n1_runtime_profiles'].get(
            approved['declared_runtime_inventory_profile_sha256'])
        if runtime is None or runtime != profile:
            raise ValueError('UNKNOWN_LEGACY_RUNTIME')
        return approved
    identity = 'index' if family == 'pr5_oracle' else 'id'
    fields = NEW_FIELDS + (('raw_report',) if family == 'pr5_oracle' else ())
    matches = [row for row in registry[family] if row[identity] == observed.get(identity)]
    if len(matches) != 1:
        raise ValueError('UNKNOWN_CASE')
    approved = matches[0]
    if any(approved.get(field) != observed.get(field) for field in fields):
        raise ValueError('UNKNOWN_CONTEXT_OR_FINGERPRINT')
    return approved


def rejects(fn):
    try:
        fn()
    except ValueError:
        return 1
    raise AssertionError('unknown selector was accepted')


def main():
    registry = json.loads((HERE / 'registry.json').read_text())
    comparison = json.loads((HERE / 'comparison.json').read_text())
    accepted = 0
    rejected = 0
    for family, expected_count in (('m1_n4', 78), ('m2_m8', 691), ('pr5_oracle', 6)):
        assert len(registry[family]) == expected_count
        for row in registry[family]:
            assert lookup(registry, family, row) == row
            accepted += 1
        sample = registry[family][0]
        for change in ('raw_stdout', 'cwd', 'loaded_modules'):
            bad = copy.deepcopy(sample)
            bad[change] = 'unknown'
            rejected += rejects(lambda bad=bad: lookup(registry, family, bad))
        rejected += rejects(lambda: lookup(registry, family, sample, platform='linux-x86_64'))
    for row in comparison['rows']:
        family = row['family']
        identity = 'index' if family == 'pr5_oracle' else 'id'
        if identity == 'id':
            matches = [case for case in registry[family] if case['id'] == row['id']]
        else:
            matches = [case for case in registry[family] if case['index'] == int(row['id'].removeprefix('step_'))]
        assert len(matches) == 1
        assert lookup(registry, family, matches[0]) == matches[0]
        assert matches[0]['process_record']['sha256'] == row['process_record_sha256']
        assert matches[0]['raw_stdout']['sha256'] == row['raw_stdout_sha256']
        assert matches[0]['raw_stderr']['sha256'] == row['raw_stderr_sha256']
    assert len(comparison['rows']) == 775
    assert len(registry['m1_n1']) == 78
    for row in registry['m1_n1']:
        profile = registry['legacy_m1_n1_runtime_profiles'][row['declared_runtime_inventory_profile_sha256']]
        assert lookup(registry, 'm1_n1', row, runtime=profile) == row
        accepted += 1
    sample = registry['m1_n1'][0]
    profile = registry['legacy_m1_n1_runtime_profiles'][sample['declared_runtime_inventory_profile_sha256']]
    for change in ('raw_stdout_sha256', 'cwd'):
        bad = copy.deepcopy(sample)
        bad[change] = 'unknown'
        rejected += rejects(lambda bad=bad: lookup(registry, 'm1_n1', bad, runtime=profile))
    rejected += rejects(lambda: lookup(registry, 'm1_n1', sample))
    rejected += rejects(lambda: lookup(registry, 'm1_n1', sample, runtime={'unknown':'runtime'}))
    assert accepted == 853 and rejected == 16
    path = HERE / 'lookup_probes.json'
    assert not path.exists()
    path.write_text(json.dumps({'schema': 'losatx-stage-g-authority-v3-lookup-probes-v2',
                                'accepted_exact_rows': accepted,
                                'new_comparison_rows_directly_resolved': 775,
                                'legacy_m1_n1_rows_v2_contract': 78,
                                'rejected_changed_or_unknown_contexts': rejected,
                                'normalization': 'NONE'}, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'accepted': accepted, 'comparison_resolved': 775, 'rejected': rejected}))

if __name__ == '__main__':
    main()
