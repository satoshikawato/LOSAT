#!/usr/bin/env python3
"""Run ONLY threaded/tblastx/p11_avclpv_psclpv exactly as check_wasm_threading_regressions.py does."""
import sys, os, json, subprocess, argparse
from pathlib import Path
sys.path.insert(0, '/mnt/c/Users/genom/GitHub/LOSAT-web-gui/LOSAT/tests')
import certify_platform_native_v010 as authority
from wasm_performance import execute

ROOT = Path('/mnt/c/Users/genom/GitHub/LOSAT-web-gui')
TESTS = ROOT / 'LOSAT/tests'
ap = argparse.ArgumentParser()
ap.add_argument('--threaded', type=Path, required=True)
ap.add_argument('--out', type=Path, required=True)
ap.add_argument('--taskset', default='28-31')
ap.add_argument('--label', default='threaded/tblastx/p11_avclpv_psclpv')
ap.add_argument('--case', default='p11_avclpv_psclpv')
ap.add_argument('--timeout', type=int, default=7200)
a = ap.parse_args()
native = Path('/home/kawato/.cache/losat-web-gui-target/s07pg-native/release/LOSAT')
oracle_dir = Path('/home/kawato/micromamba/bin')
out = a.out.resolve(); out.mkdir(parents=True, exist_ok=True)
head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
catalog = authority.load_catalog(ROOT, head)
oracles = {p: (oracle_dir / p).resolve() for p in ['blastn', 'blastp', 'tblastx']}
steps = authority.build_steps(ROOT, out, catalog, native, oracles)
authority.stage_required_fixtures(ROOT, head, steps, out)
environment = {k: v for k, v in os.environ.items() if not k.startswith(('LOSAT_', 'RAYON_')) and k != 'BL2SEQ_LEGACY'}
environment.update(LC_ALL='C', NODE_NO_WARNINGS='1')
prefix = ['node', str(TESTS / 'run_losat_wasi_threads.js'), str(a.threaded.resolve())]
step = [s for s in steps if s.kind == 'matrix' and s.program == 'tblastx' and s.case_id == a.case]
assert len(step) == 1, step
step = step[0]
command = list(step.command)
ci = command.index('-num_threads') + 1
current = [*prefix, *command[1:]]
current[current.index('-num_threads') + 1] = '4'
env = {**environment, **dict(step.environment), 'LOSAT_WASI_THREADS_DEBUG': '1'}
if a.taskset:
    # pin via taskset on the python process itself would also pin children; but execute() wraps /usr/bin/time, so use affinity inheritance
    os.sched_setaffinity(0, set(range(int(a.taskset.split('-')[0]), int(a.taskset.split('-')[1]) + 1)))
print('affinity', sorted(os.sched_getaffinity(0)), flush=True)
print('cmd', json.dumps(current), flush=True)
result = execute(current, ROOT, out / a.label, env, a.timeout)
result.update(label=a.label, expected_sha256=step.expected_losat_sha256, raw_equal=result['status']=='PASS' and result['raw_output_sha256']==step.expected_losat_sha256)
(out / 'result_summary.json').write_text(json.dumps(result, indent=2))
print(json.dumps({k: result.get(k) for k in ['status','wall_seconds','cpu_user_seconds','cpu_system_seconds','peak_rss_bytes','raw_output_sha256','expected_sha256','raw_equal']}, indent=1), flush=True)
