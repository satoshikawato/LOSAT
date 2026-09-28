#!/usr/bin/env python3
"""Validate bound v3 authority, portable receipts and independent audit.

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
import sys
import tarfile

HERE = Path(__file__).resolve().parent
PRIOR = Path('/mnt/c/Users/genom/GitHub/LOSAT/docs/evidence/losatx_stage_g_hold_completion_clarified/run-20260928T104629Z')
WORK = Path('/mnt/c/users/genom/github/losatx-blastx-v020-recovery')
EXCLUDED = {'MANIFEST.json', 'binding.json', 'audit_raw.md', 'audit.json', 'final_result.json', 'final_stdout', 'final_stderr'}

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    binding = json.loads((HERE / 'binding.json').read_text())
    manifest = json.loads((HERE / 'MANIFEST.json').read_text())
    audit = json.loads((HERE / 'audit.json').read_text())
    assert sha(PRIOR / 'binding.json') == binding['prior_binding_sha256']
    assert sha(PRIOR / 'audit_raw.md') == binding['prior_independent_audit_sha256']
    failed = HERE.parent / 'run-20260928T113843Z/binding.json'
    assert sha(failed) == binding['failed_predecessor_binding_sha256']
    stale = HERE.parent / 'run-20260928T115500Z/binding.json'
    assert sha(stale) == binding['stale_prose_predecessor_binding_sha256']
    stale_runner = HERE.parent / 'run-20260928T120321Z/binding.json'
    assert sha(stale_runner) == binding['stale_runner_plan_predecessor_binding_sha256']
    stale_repro = HERE.parent / 'run-20260928T121215Z/binding.json'
    assert sha(stale_repro) == binding['stale_reproduction_predecessor_binding_sha256']
    assert sha(HERE / 'MANIFEST.json') == binding['manifest_sha256']
    names = {p.name for p in HERE.iterdir() if p.is_file() and p.name not in EXCLUDED}
    assert names == set(manifest['files'])
    for name, info in manifest['files'].items():
        path = HERE / name
        assert path.stat().st_size == info['size'] and sha(path) == info['sha256'], name
    assert sha(HERE / 'registry.json') == binding['registry_sha256']
    assert sha(HERE / 'comparison.json') == binding['comparison_sha256']
    assert sha(HERE / 'approved_windows_receipts.tar.gz') == binding['portable_receipts_sha256']
    assert audit['binding_sha256'] == sha(HERE / 'binding.json')
    assert audit['audit_sha256'] == sha(HERE / 'audit_raw.md')
    closure = json.loads((PRIOR / 'candidate_source_closure.json').read_text())
    assert len(closure) == 485
    assert sha(PRIOR / 'candidate_source_closure.json') == binding['candidate_source_closure_sha256']
    for relative, digest in closure.items():
        assert sha(WORK / relative) == digest, relative
    registry = json.loads((HERE / 'registry.json').read_text())
    comparison = json.loads((HERE / 'comparison.json').read_text())
    coverage = json.loads((HERE / 'coverage_delta.json').read_text())
    assert [len(registry[key]) for key in ('m1_n1', 'm1_n4', 'm2_m8', 'pr5_oracle')] == [78, 78, 691, 6]
    assert len(comparison['rows']) == binding['new_official_windows_fingerprints_registered'] == 775
    from lookup_v3 import lookup
    for row in comparison['rows']:
        family = row['family']
        if family == 'pr5_oracle':
            matches = [x for x in registry[family] if x['index'] == int(row['id'].removeprefix('step_'))]
        else:
            matches = [x for x in registry[family] if x['id'] == row['id']]
        assert len(matches) == 1 and lookup(registry, family, matches[0]) == matches[0]
        assert matches[0]['process_record']['sha256'] == row['process_record_sha256']
    probes = json.loads((HERE / 'lookup_probes.json').read_text())
    assert probes['accepted_exact_rows'] == 853 and probes['new_comparison_rows_directly_resolved'] == 775
    assert comparison['normalization_used_for_acceptance'] is False
    assert coverage['not_executed_before_and_after'] == 6912
    archive = HERE / 'approved_windows_receipts.tar.gz'
    index = json.loads((HERE / 'receipt_index.json').read_text())
    assert sha(archive) == index['archive_sha256']
    expected = set(index['source_path_to_sha256'].values())
    observed = set()
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar:
            assert member.name.startswith('sha256/') and member.isfile()
            digest = member.name.removeprefix('sha256/')
            assert len(digest) == 64 and digest not in observed
            stream = tar.extractfile(member)
            assert stream is not None
            h = hashlib.sha256()
            for chunk in iter(lambda: stream.read(1 << 20), b''):
                h.update(chunk)
            assert h.hexdigest() == digest
            observed.add(digest)
    assert observed == expected and len(observed) == index['unique_objects'] == 1289
    assert binding['frozen_pr5_gate_a_history'] == 'HARD_FAIL_UNCHANGED'
    assert binding['normalization_used_for_acceptance'] is False
    assert binding['implementation_acceptance'] == 'HARD_FAIL' and binding['release'] == 'HOLD'
    result = {'schema': 'losatx-stage-g-authority-v3-final-v1',
              'binding_sha256': sha(HERE / 'binding.json'),
              'audit_sha256': sha(HERE / 'audit_raw.md'),
              'artifact_integrity': 'PASS_FOR_EXACT_ROSTER_AND_LISTED_LIVE_PINS',
              'implementation_acceptance': 'HARD_FAIL', 'release': 'HOLD',
              'registered_new_official_windows_fingerprints': 775,
              'new_sakai_authority_steps': [11, 52, 53],
              'frozen_pr5_gate_a_history': 'HARD_FAIL_UNCHANGED',
              'portable_receipt_objects': len(observed),
              'unexecuted_unique_matrix_rows': 6912,
              'missing_actual_native_targets': coverage['missing_actual_native_targets'],
              'normalization_used_for_acceptance': False}
    (HERE / 'final_result.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, sort_keys=True))
    return 1

if __name__ == '__main__':
    sys.exit(main())
