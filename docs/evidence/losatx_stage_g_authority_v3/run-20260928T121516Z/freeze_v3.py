#!/usr/bin/env python3
"""Freeze the post-approval exact-byte authority roster.

NCBI c++/src/objtools/align_format/tabular.cpp:1100-1108:
// void CBlastTabularInfo::Print() {
//     ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
//         x_PrintField(*iter);
//     }
//     m_Ostream << "\n";
// }
"""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRIOR_BINDING = 'a68d6e657231b20c3251597d9a1238df53df73d9f69c6abd4dcdca5f744f5122'
PRIOR_AUDIT = '5f912b618e6830f363aa155bce30c4211b27f5231e2468f93686639ccdd5716b'
EXCLUDED = {'MANIFEST.json', 'binding.json', 'audit_raw.md', 'audit.json', 'final_result.json', 'final_stdout', 'final_stderr'}

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    assert not (HERE / 'MANIFEST.json').exists()
    assert not (HERE / 'binding.json').exists()
    roster = {}
    for path in sorted(HERE.iterdir()):
        if path.is_dir() or path.name in EXCLUDED:
            continue
        roster[path.name] = {'sha256': sha(path), 'size': path.stat().st_size}
    assert len(roster) >= 20
    manifest = {'schema': 'losatx-stage-g-authority-v3-manifest-v1', 'files': roster}
    (HERE / 'MANIFEST.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    binding = {
        'schema': 'losatx-stage-g-authority-v3-binding-v1',
        'authority_version': 'losatx-stage-g-authority-v3.0.0',
        'prior_binding_sha256': PRIOR_BINDING,
        'prior_independent_audit_sha256': PRIOR_AUDIT,
        'failed_predecessor_binding_sha256': 'faa9416b04a8611fdc4ab4fa19699e46787277c00d7ef18986c1c5907cd6a052',
        'stale_prose_predecessor_binding_sha256': '38eb06edd0fdbb354b04c610d2f2e735a789c92ef44b60bea2b6986f6a4fabeb',
        'stale_runner_plan_predecessor_binding_sha256': 'a9f5f9a0584b8cf6f0d7e8c0134a216962e3e92f5654b4bd8a80e6ae7a4eb6bb',
        'stale_reproduction_predecessor_binding_sha256': 'c92908d0167f2c083771c054fdad707886102c05b467b840f492f0d5f823d49b',
        'manifest_sha256': sha(HERE / 'MANIFEST.json'),
        'registry_sha256': roster['registry.json']['sha256'],
        'comparison_sha256': roster['comparison.json']['sha256'],
        'portable_receipts_sha256': roster['approved_windows_receipts.tar.gz']['sha256'],
        'candidate_head_before_pr': '1c131e6027ff7fc0595102f5c0a9fcc5f17b7d56',
        'candidate_source_closure_sha256': 'd5a3dfd3ac00afdf5ac1589f1c4d962ce455370a41d2ba86aeb14a14bd8fc1dd',
        'new_official_windows_fingerprints_registered': 775,
        'approved_sakai_steps_in_new_authority': [11, 52, 53],
        'frozen_pr5_gate_a_history': 'HARD_FAIL_UNCHANGED',
        'normalization_used_for_acceptance': False,
        'artifact_integrity': 'PENDING_INDEPENDENT_AUDIT',
        'implementation_acceptance': 'HARD_FAIL',
        'release': 'HOLD',
    }
    (HERE / 'binding.json').write_text(json.dumps(binding, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'binding_sha256': sha(HERE / 'binding.json'), 'manifest_files': len(roster), 'manifest_sha256': sha(HERE / 'MANIFEST.json')}))

if __name__ == '__main__':
    main()
