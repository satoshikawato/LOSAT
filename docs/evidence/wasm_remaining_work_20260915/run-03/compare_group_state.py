# NCBI reference: c++/src/algo/blast/core/link_hsps.c:603-982
# while (number_of_hsps > 0) { ... H->hsp->evalue = evalue; ... }
# Full actual group transcripts; compare each target to its own frozen baseline.
from pathlib import Path
import argparse, subprocess, os, json, hashlib
e = Path(__file__).resolve().parent
p = argparse.ArgumentParser()
p.add_argument('--label', default='group-state')
p.add_argument('--kinds', nargs='+', default=['native', 'serial', 'threaded'])
p.add_argument('--candidate', default='I1')
a = p.parse_args()
out = e / a.label
out.mkdir()
env = {k:v for k,v in os.environ.items() if not k.startswith(('LOSAT_', 'NODE_', 'RAYON_'))}
env.update(NODE_NO_WARNINGS='1', LOSAT_WASI_THREAD_CAP='8', LOSAT_WASM_MEMORY_MAXIMUM_PAGES='16384')
records = []
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
for kind in a.kinds:
    for trace in ['off', 'on']:
        for version in ['baseline', a.candidate]:
            d = out / kind / trace / version
            d.mkdir(parents=True)
            artifact = e / 'work' / (version + 'state') / 'artifacts' / ('losat-' + kind + '-command' + ('.wasm' if kind != 'native' else ''))
            cmd = [str(artifact)] if kind == 'native' else ['node', str(e / 'work/baseline-source/LOSAT/tests' / ('run_losat_wasi.js' if kind == 'serial' else 'run_losat_wasi_threads.js')), str(artifact)]
            runenv = dict(env)
            if trace == 'on': runenv['LOSAT_STATE_TRACE'] = '1'
            record = dict(kind=kind, trace=trace, version=version, argv=cmd, artifact_sha256=digest(artifact), status='RUNNING')
            records.append(record)
            (out / 'runs.json').write_text(json.dumps(records, indent=2) + '\n')
            # NCBI link_hsps.c:603-982: diagnostic transcript only. Buffer pipe
            # reads to avoid a mounted-filesystem write per Rust println.
            # Preserve the same process timeout; these runs are never timings.
            try:
                r = subprocess.run(cmd, env=runenv, capture_output=True, timeout=120)
            except subprocess.TimeoutExpired as exc:
                (d/'stdout.txt').write_bytes(exc.stdout or b'')
                (d/'stderr.txt').write_bytes(exc.stderr or b'')
                record['status']='TIMEOUT'
                (out/'runs.json').write_text(json.dumps(records,indent=2)+'\n')
                raise
            (d/'stdout.txt').write_bytes(r.stdout)
            (d/'stderr.txt').write_bytes(r.stderr)
            record.update(returncode=r.returncode, stdout_sha256=digest(d/'stdout.txt'), stderr_sha256=digest(d/'stderr.txt'))
            text = (d / 'stdout.txt').read_text()
            record.update(cases=sum(line.startswith('STATE CASE ') for line in text.splitlines()), rounds=sum(line.startswith('STATE ENTER ') for line in text.splitlines()), visits=sum(line.startswith('STATE VISIT ') for line in text.splitlines()))
            assert r.returncode == 0 and record['cases'] == 352, record
            record['status'] = 'PASS'
            (out / 'runs.json').write_text(json.dumps(records, indent=2) + '\n')
        for stream in ['stdout', 'stderr']:
            assert (out/kind/trace/'baseline'/(stream+'.txt')).read_bytes() == (out/kind/trace/a.candidate/(stream+'.txt')).read_bytes(), (kind,trace,stream)
        print(kind, trace, 'EXACT', record['cases'], 'cases', record['rounds'], 'rounds', record['visits'], 'visits', flush=True)
(out / 'comparison.json').write_text(json.dumps(dict(status='PASS', targets=a.kinds, cases=352, trace=['off','on'], comparison='byte-exact complete stdout/stderr within each target; typed sentinel indices serialized as Option; every f64 retained as bits'), indent=2)+'\n')
