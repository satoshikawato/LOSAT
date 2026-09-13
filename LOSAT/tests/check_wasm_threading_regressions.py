#!/usr/bin/env python3
"""Apply frozen PR5 raw-output expectations to current explicit artifacts.

This is a working-tree regression run, not a clean-SHA hosted certification.
Gate A expectations and retained Linux oracle fingerprints are never updated.
The existing hosted Gate B remains platform-specific and independently required.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import threading
from concurrent.futures import ThreadPoolExecutor
import certify_platform_native_v010 as authority
from wasm_performance import execute, digest, validate_thread_evidence

ROOT = Path(__file__).resolve().parents[2]
TESTS = ROOT / 'LOSAT/tests'

# NCBI reference: c++/src/objtools/align_format/tabular.cpp:1098-1108
# x_PrintField(*iter); m_Ostream << "\n";
# Preserve frozen fixture bytes/paths, formatting and order; no normalization.
def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['native', 'serial', 'threaded', 'oracle-dir', 'output-dir']:
        p.add_argument('--' + name, type=Path, required=True)
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
    metadata = dict(head=head, scope=__doc__, per_search_timeout_seconds=args.timeout_seconds, artifacts={name: {'path': str(value.resolve()), 'sha256': digest(value)} for name, value in [('native', args.native), ('serial', args.serial), ('threaded', args.threaded)]}, authority_sha256=native_authority.file_sha256)
    (out / 'metadata.json').write_text(json.dumps(metadata, indent=2))
    def record(label, command, expected, env, threads=None, kind=None):
        result = execute(command, ROOT, out / label, env, args.timeout_seconds)
        result.update(label=label, expected_sha256=expected)
        result['raw_equal'] = result['status'] == 'PASS' and result['raw_output_sha256'] == expected
        with record_lock:
            records.append(result); (out / 'runs.json').write_text(json.dumps(records, indent=2))
        if not result['raw_equal']: raise RuntimeError(f'{label}: frozen raw output mismatch or execution failure: {result}')
        if threads is not None:
            validate_thread_evidence((out / label / 'stderr.txt').read_text(), threads, kind)
        print(label, 'PASS', flush=True)
    for step in steps:
        if step.kind == 'matrix':
            command = list(step.command)
            count_index = command.index('-num_threads') + 1
            original_n = int(command[count_index])
            for kind, prefix in [('native', [str(args.native.resolve())]), ('serial', [args.node, str(TESTS/'run_losat_wasi.js'), str(args.serial.resolve())]), ('threaded', [args.node, str(TESTS/'run_losat_wasi_threads.js'), str(args.threaded.resolve())])]:
                if kind == 'serial' and original_n != 1: continue  # explicitly non-applicable frozen native-thread4 rows
                current = [*prefix, *command[1:]]
                n = 4 if kind == 'threaded' else original_n
                current[current.index('-num_threads') + 1] = str(n)
                pending.append((f'{kind}/{step.program}/{step.case_id}', current, step.expected_losat_sha256, {**environment, **dict(step.environment), 'LOSAT_WASI_THREADS_DEBUG':'1'}, n, kind))
        elif step.kind == 'oracle':
            references = [row['retained_linux_raw_sha256'] for row in native_authority.document['diagnostic_metadata'] if row['program'] == step.program and row['case_id'] == step.case_id]
            assert len(set(references)) == 1
            pending.append((f'retained-linux-oracle/{step.program}/{step.case_id}', list(step.command), references[0], environment))
        else:
            pending.append((f'repeatability/{step.program}/{step.case_id}/{step.run_index}', list(step.command), step.expected_losat_sha256, {**environment, **dict(step.environment)}))
    # NCBI reference: c++/src/objtools/align_format/tabular.cpp:1098-1108
    # x_PrintField(*iter); m_Ostream << "\n";
    # Independent oracle processes may overlap; each output is compared without reordering.
    # This is correctness execution, never a performance sample.
    metadata['concurrent_test_processes']=args.jobs
    (out/'metadata.json').write_text(json.dumps(metadata,indent=2))
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures=[executor.submit(record,*job) for job in pending]
        for future in futures: future.result()
    print(f'{len(records)} frozen regression records PASS', flush=True)

if __name__ == '__main__': main()
