#!/usr/bin/env python3
"""Bind the specifically approved Session G Windows and Sakai authorities.

NCBI c++/src/objtools/align_format/tabular.cpp:1100-1108:
// void CBlastTabularInfo::Print() {
//     ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
//         x_PrintField(*iter);
//     }
//     m_Ostream << "\\n";
// }
The authority compares the resulting raw bytes. It never normalizes newlines.
"""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = Path('/mnt/c/Users/genom/GitHub/LOSAT/docs/evidence')
PRIOR = ROOT / 'losatx_stage_g_hold_completion_clarified/run-20260928T104629Z'
V2 = ROOT / 'losatx_stage_g_authority_v2/run-20260928T013630Z'
VERSION = 'losatx-stage-g-authority-v3.0.0'
PINS = {
    'binding.json': 'a68d6e657231b20c3251597d9a1238df53df73d9f69c6abd4dcdca5f744f5122',
    'audit_raw.md': '5f912b618e6830f363aa155bce30c4211b27f5231e2468f93686639ccdd5716b',
    'windows_authority_proposal.json': '802965219a96074a34e94968064ae761fd00d5a1c1549ecb5f86010f06fe60dd',
    'windows_authority_replay_addendum.json': 'db3f46af526dd0869428d972051e13217a229056006d3666af2aa6840072ed5b',
    'windows_m1_n4_authority_proposal.json': 'fe706199f0de37e56f52da2ca505433b45663ccabedb8ea90f1efacc1f70683b',
    'windows_pr5_oracle_authority_proposal.json': '34b4b8ff348db1c550e3365ad5ea27feb93992b8775a62c1f63568195a33da23',
    'sakai_authority_continuation.json': '86ff6c32cae1ea355ff3cc2a8e99b283537a5771feaa9794ef287b4d7937c057',
}
OLD_REGISTRY_SHA = 'de0fc2fb0bbc9221d4adf66d5ee5a6ab9a7cab2dd8b3c2d71d93fdcdf78a1466'
SOURCE_CLOSURE_SHA = 'd5a3dfd3ac00afdf5ac1589f1c4d962ce455370a41d2ba86aeb14a14bd8fc1dd'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write(path, obj):
    path = Path(path)
    assert not path.exists(), f'already frozen: {path}'
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False) + '\n')


def checked(path, digest):
    path = Path(path)
    assert path.is_file(), f'missing: {path}'
    assert sha(path) == digest, f'changed: {path}'
    return path


def check_raw(item):
    path = checked(item['path'], item['sha256'])
    if 'size' in item:
        assert path.stat().st_size == item['size'], f'changed size: {path}'
    return path


def verify_case(case, family, package):
    receipt = read(checked(case['process_record']['path'], case['process_record']['sha256']))
    assert receipt['argv'] == case['argv']
    assert receipt['cwd'] == case['cwd']
    assert receipt['effective_environment'] == case['effective_environment']
    assert receipt['exit'] == case['exit'] == 0
    assert receipt['exe'] == case['argv'][0]
    assert receipt['debug_create_events'] == 1
    assert receipt['debug_unresolved_events'] == 0
    assert receipt['debug_errors'] == []
    assert receipt['unhashable_modules'] == []
    assert receipt['executable_before'] == receipt['executable_after']
    assert len(receipt['loaded_modules']) == 27
    # A repeated debugger load event is still one distinct loaded module.
    assert len(set(receipt['loaded_module_paths'])) == 27
    assert set(receipt['loaded_modules']) == set(receipt['loaded_module_paths'])
    assert receipt['loaded_modules'] == receipt['loaded_modules_after']
    assert receipt['executable_after'] == receipt['loaded_modules'][receipt['exe']]
    if family != 'pr5_oracle':
        assert receipt['executable_after'] == package['blastx_executable_sha256']
    assert receipt['stdout_sha256'] == case['raw_stdout']['sha256']
    assert receipt['stderr_sha256'] == case['raw_stderr']['sha256']
    for kind in ('raw_stdout', 'raw_stderr'):
        check_raw(case[kind])
    if 'raw_report' in case:
        check_raw(case['raw_report'])
    for pin in case['inputs'].values():
        source = pin.get('source_path', pin.get('mirror'))
        path = checked(source, pin['sha256'])
        assert path.stat().st_size == pin['size']
        assert pin['relative_name']
    if 'loaded_modules' in case:
        assert receipt['loaded_modules'] == case['loaded_modules']
    else:
        assert receipt['loaded_modules'] == family_profile
    return {
        'id': case.get('id', f"step_{case['index']}" if 'index' in case else None),
        'program': case['program'],
        'family': family,
        'process_record_sha256': case['process_record']['sha256'],
        'raw_stdout_sha256': case['raw_stdout']['sha256'],
        'raw_stderr_sha256': case['raw_stderr']['sha256'],
        'raw_report_sha256': case.get('raw_report', {}).get('sha256'),
        'module_count': len(receipt['loaded_modules']),
    }


def prepare():
    for name, digest in PINS.items():
        checked(PRIOR / name, digest)
    checked(V2 / 'windows_registry.json', OLD_REGISTRY_SHA)
    closure = read(checked(PRIOR / 'candidate_source_closure.json', SOURCE_CLOSURE_SHA))
    assert len(closure) == 485
    work = Path('/mnt/c/users/genom/github/losatx-blastx-v020-recovery')
    for relative, digest in closure.items():
        checked(work / relative, digest)

    old = read(V2 / 'windows_registry.json')
    assert len(old['entries']) == 78
    assert all(row['threads'] == 1 for row in old['entries'])
    m28 = read(PRIOR / 'windows_authority_proposal.json')
    n4 = read(PRIOR / 'windows_m1_n4_authority_proposal.json')
    pr5 = read(PRIOR / 'windows_pr5_oracle_authority_proposal.json')
    sakai = read(PRIOR / 'sakai_authority_continuation.json')
    addendum = read(PRIOR / 'windows_authority_replay_addendum.json')
    assert [len(x['cases']) for x in (m28, n4, pr5)] == [691, 78, 6]
    assert m28['package'] == n4['package'] == pr5['package']
    assert addendum['base_windows_authority_proposal_sha256'] == PINS['windows_authority_proposal.json']
    assert len(addendum['rows']) == 691
    assert {x['id'] for x in addendum['rows']} == {x['id'] for x in m28['cases']}
    assert sakai['affected_frozen_gate_a_steps'] == [11, 52, 53]
    assert sakai['historical_frozen_gate_a'] == 'HARD_FAIL_UNCHANGED'
    assert sakai['candidate_source_closure_sha256'] == SOURCE_CLOSURE_SHA
    assert sakai['candidate_raw_sha256'] != sakai['frozen_raw_sha256']
    assert read(PRIOR / 'final_result.json')['release'] == 'HOLD'

    global family_profile
    family_profile = m28['process_module_profile']
    assert len(family_profile) == 27
    rows = []
    for family, source in (('m2_m8', m28), ('m1_n4', n4), ('pr5_oracle', pr5)):
        for case in source['cases']:
            rows.append(verify_case(case, family, m28['package']))
    assert len(rows) == 775
    assert len({(r['family'], r['id']) for r in rows}) == 775
    summary = read(PRIOR / 'pr5_windows_current/summary.json')
    by_step = {x['index']: x for x in summary['rows']}
    assert summary['frozen_gate_a_fail_indices'] == [11, 52, 53]
    for step in (11, 52, 53):
        assert by_step[step]['raw_output_sha256'] == sakai['candidate_raw_sha256']
        assert by_step[step]['expected_frozen_gate_a_sha256'] == sakai['frozen_raw_sha256']
    for case in pr5['cases']:
        assert by_step[case['index']]['raw_output_sha256'] == case['raw_report']['sha256']
    rust = read(PRIOR / 'final_windows_candidate_result.json')
    assert rust['exact_raw_rows'] == rust['executed'] == 956
    assert rust['normalization_used_for_acceptance'] is False
    assert len([x for x in rust['rows'] if x['id'] in {c['id'] for c in m28['cases']}
                and x['raw_stdout_equal'] and x['raw_stderr_equal'] and x['exit_equal']]) == 691

    registry = {
        'schema': VERSION,
        'approval': 'USER_APPROVED_EXACT_FOUR_PROPOSALS_2026-09-28',
        'prior_approved_m1_n1_registry_sha256': OLD_REGISTRY_SHA,
        'prior_binding_sha256': PINS['binding.json'],
        'prior_independent_audit_sha256': PINS['audit_raw.md'],
        'proposal_sha256': {k: PINS[k] for k in PINS if 'proposal' in k or k == 'sakai_authority_continuation.json'},
        'package': m28['package'],
        'm1_n1': old['entries'],
        'legacy_m1_n1_runtime_profiles': old['runtime_profiles'],
        'm1_n4': n4['cases'],
        'm2_m8': [dict(case, loaded_modules=family_profile) for case in m28['cases']],
        'runtime_profiles': {m28['process_module_profile_sha256']: family_profile,
                             **n4['process_module_profiles']},
        'pr5_oracle': pr5['cases'],
        'sakai_new_gate_a': {str(i): sakai['candidate_raw_sha256'] for i in (11, 52, 53)},
        'frozen_pr5_gate_a': 'HISTORICAL_HARD_FAIL_UNCHANGED',
        'losat_expected': 'FIXED_LINUX_RAW_UNCHANGED',
        'normalization': 'NONE',
        'unknown_context_or_fingerprint': 'HARD_FAIL',
        'missing_actual_native_runner': 'HARD_FAIL',
    }
    write(HERE / 'registry.json', registry)
    write(HERE / 'comparison.json', {
        'schema': VERSION,
        'rows': rows,
        'registered_new_official_windows_fingerprints': 775,
        'prior_approved_windows_fingerprints': 78,
        'new_sakai_gate_a_steps': [11, 52, 53],
        'new_authority_comparison': 'PASS_FOR_775_EXACT_FIXED_RECEIPTS_AND_3_SAKAI_STEPS',
        'frozen_pr5_gate_a_history': 'HARD_FAIL_UNCHANGED',
        'rust_fixed_linux_expected_unchanged': True,
        'normalization_used_for_acceptance': False,
        'release': 'HOLD',
    })
    print(json.dumps({'registry_sha256': sha(HERE / 'registry.json'), 'rows': len(rows), 'release': 'HOLD'}))


if __name__ == '__main__':
    prepare()
