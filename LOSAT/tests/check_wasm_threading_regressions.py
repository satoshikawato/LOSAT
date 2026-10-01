#!/usr/bin/env python3
"""Apply frozen PR5 raw-output expectations to current explicit artifacts.

This is a working-tree regression run, not a clean-SHA hosted certification.
Gate A expectations and retained Linux oracle fingerprints are never updated.
The existing hosted Gate B remains platform-specific and independently required.
Every search runs; known mismatches listed in frozen_mismatch_allowlist.json pass
only with their listed output, and the run fails on any other mismatch, on an
execution or thread-evidence failure, and on a stale allow-list entry.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import threading
from concurrent.futures import ThreadPoolExecutor
import certify_platform_native_v010 as authority
from frozen_allowlist import Allowlist
from wasm_performance import execute, digest, validate_thread_evidence

ROOT = Path(__file__).resolve().parents[2]
TESTS = ROOT / 'LOSAT/tests'

# NCBI reference: c++/src/objtools/align_format/tabular.cpp:1098-1108
# x_PrintField(*iter); m_Ostream << "\n";
# Preserve frozen fixture bytes/paths, formatting and order; no normalization.
def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['native', 'threaded', 'oracle-dir', 'output-dir']:
        p.add_argument('--' + name, type=Path, required=True)
    # NCBI reference: c++/include/algo/blast/blastinput/blast_args.hpp:1290-1296
    # m_NumThreads = CThreadable::kMinNumThreads; m_MTMode = eNotSupported;
    p.add_argument('--serial', type=Path, help='opt in to serial compatibility checks')
    p.add_argument('--node', default='node')
    p.add_argument('--jobs',type=int,default=1,choices=[1,2,3])
    p.add_argument('--timeout-seconds', type=int, default=3600, help='positive per-search correctness deadline; long genomes may exceed 900 seconds')
    args = p.parse_args()
    if args.timeout_seconds <= 0: p.error('--timeout-seconds must be positive')
    out = args.output_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    catalog = authority.load_catalog(ROOT, head)
    oracles = {program: (args.oracle_dir / program).resolve() for program in ['blastn', 'blastp', 'tblastx']}
    steps = authority.build_steps(ROOT, out, catalog, args.native.resolve(), oracles)
    authority.stage_required_fixtures(ROOT, head, steps, out)
    native_authority = authority.load_native_authority(ROOT, head, TESTS / 'ncbi_platform_variance_v010.json')
    environment = {k: v for k, v in os.environ.items() if not k.startswith(('LOSAT_', 'RAYON_')) and k != 'BL2SEQ_LEGACY'}
    environment.update(LC_ALL='C', NODE_NO_WARNINGS='1')
    records = []; record_lock=threading.Lock(); pending=[]
    # NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:173-180
    # (*thread)->Run(); (*thread)->Join(&result);
    metadata = dict(head=head, scope=__doc__, per_search_timeout_seconds=args.timeout_seconds, artifacts={name: {'path': str(value.resolve()), 'sha256': digest(value)} for name, value in [('native', args.native), ('serial', args.serial), ('threaded', args.threaded)] if value}, authority_sha256=native_authority.file_sha256)
    prefixes = [('native', [str(args.native.resolve())]),
                ('threaded', [args.node, str(TESTS/'run_losat_wasi_threads.js'), str(args.threaded.resolve())])]
    if args.serial:
        prefixes.append(('serial', [args.node, str(TESTS/'run_losat_wasi.js'), str(args.serial.resolve())]))
    (out / 'metadata.json').write_text(json.dumps(metadata, indent=2))
    allowlist = Allowlist()
    def record(program, case_id, label, command, expected, env, threads=None, kind=None):
        result = execute(command, ROOT, out / label, env, args.timeout_seconds)
        result.update(label=label, program=program, case_id=case_id, expected_sha256=expected)
        executed = result['status'] == 'PASS'
        result['raw_equal'] = executed and result['raw_output_sha256'] == expected
        result['allowlist'] = allowlist.classify(program, case_id, expected, result.get('raw_output_sha256'), executed)
        result['thread_evidence'] = None
        if threads is not None and executed:
            try:
                validate_thread_evidence((out / label / 'stderr.txt').read_text(), threads, kind)
                result['thread_evidence'] = 'PASS'
            except Exception as error:  # recorded and reported with every other failure
                result['thread_evidence'] = f'FAIL: {error}'
        ok = result['allowlist'] in ('match', 'allowed') and result['thread_evidence'] in (None, 'PASS')
        with record_lock:
            records.append(result); (out / 'runs.json').write_text(json.dumps(records, indent=2))
        print(label, 'PASS' if result['allowlist'] == 'match' and ok else 'ALLOWED (known mismatch)' if ok else f"FAIL ({result['allowlist']}, thread evidence {result['thread_evidence']}, status {result['status']})", flush=True)
    for step in steps:
        if step.kind == 'matrix':
            command = list(step.command)
            count_index = command.index('-num_threads') + 1
            original_n = int(command[count_index])
            # NCBI reference: c++/src/objtools/align_format/tabular.cpp:1098-1108
            # x_PrintField(*iter); m_Ostream << "\n";
            for kind, prefix in prefixes:
                if kind == 'serial' and original_n != 1: continue  # explicitly non-applicable frozen native-thread4 rows
                current = [*prefix, *command[1:]]
                n = 4 if kind == 'threaded' else original_n
                current[current.index('-num_threads') + 1] = str(n)
                pending.append((step.program, step.case_id, f'{kind}/{step.program}/{step.case_id}', current, step.expected_losat_sha256, {**environment, **dict(step.environment), 'LOSAT_WASI_THREADS_DEBUG':'1'}, n, kind))
        elif step.kind == 'oracle':
            references = [row['retained_linux_raw_sha256'] for row in native_authority.document['diagnostic_metadata'] if row['program'] == step.program and row['case_id'] == step.case_id]
            assert len(set(references)) == 1
            pending.append((step.program, step.case_id, f'retained-linux-oracle/{step.program}/{step.case_id}', list(step.command), references[0], environment))
        else:
            pending.append((step.program, step.case_id, f'repeatability/{step.program}/{step.case_id}/{step.run_index}', list(step.command), step.expected_losat_sha256, {**environment, **dict(step.environment)}))
    # NCBI reference: c++/src/objtools/align_format/tabular.cpp:1098-1108
    # x_PrintField(*iter); m_Ostream << "\n";
    # Independent oracle processes may overlap; each output is compared without reordering.
    # This is correctness execution, never a performance sample.
    metadata['concurrent_test_processes']=args.jobs
    (out/'metadata.json').write_text(json.dumps(metadata,indent=2))
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures=[executor.submit(record,*job) for job in pending]
        for future in futures: future.result()
    failures = [r['label'] for r in records if r['allowlist'] not in ('match', 'allowed') or r['thread_evidence'] not in (None, 'PASS')]
    allowed = [r['label'] for r in records if r['allowlist'] == 'allowed']
    executed = {(r['program'], r['case_id']) for r in records}
    unexecuted = allowlist.unexecuted(executed, {r['program'] for r in records})
    summary = dict(records=len(records), allowed_known_mismatches=allowed, failures=failures,
                   stale_unexecuted_allowlist_entries=[f'{p}/{c}' for p, c in unexecuted], allowlist=str(allowlist.path.relative_to(ROOT)))
    (out / 'summary.json').write_text(json.dumps(summary, indent=2))
    for label in allowed: print(f'allowed known mismatch: {label}', flush=True)
    for label in failures: print(f'FAILED: {label}', flush=True)
    for program, case_id in unexecuted: print(f'FAILED: allow-list entry {program}/{case_id} was not executed (remove it)', flush=True)
    if failures or unexecuted:
        raise SystemExit(f'{len(failures)} failed and {len(unexecuted)} unexecuted allow-list entries out of {len(records)} frozen regression records')
    print(f'{len(records)} frozen regression records PASS ({len(allowed)} allowed known mismatches)', flush=True)

if __name__ == '__main__': main()
