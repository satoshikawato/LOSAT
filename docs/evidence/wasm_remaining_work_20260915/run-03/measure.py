# NCBI reference: c++/src/app/blast/blastn_app.cpp:59-67
# m_StopWatch.Start(); m_StopWatch.Elapsed();
# NCBI tabular.cpp:1098-1108: x_PrintField(*iter); m_Ostream << "\n";
# Measure successful byte-identical searches; preserve warmups and all pairs.
import argparse, csv, json, os, re, shutil, statistics, sys
from pathlib import Path

e = Path(__file__).resolve().parent
b = e / 'work'
root = e.parents[3]
tests = b / 'baseline-source/LOSAT/tests'
sys.path.insert(0, str(tests))
from wasm_performance import execute, digest

p = argparse.ArgumentParser()
p.add_argument('--label', required=True)
p.add_argument('--candidate', default='I1')
p.add_argument('--kinds', nargs='+', default=['threaded'])
p.add_argument('--threads', nargs='+', type=int, default=[8, 1])
p.add_argument('--cases', nargs='+', default=['MjeNMV.MelaMJNV.tlosatx', 'AP027280.AP027280.tlosatx'])
p.add_argument('--diagnostic-background', help='Record a known background workload; results cannot authorize adoption')
a = p.parse_args()
out = e / a.label
out.mkdir()
node = str(Path(shutil.which('node')).resolve())
rows = list(csv.DictReader((tests / 'comparison_cases.tsv').open(), delimiter='\t'))
env = {k: v for k, v in os.environ.items() if not k.startswith(('LOSAT_', 'NODE_', 'RAYON_')) and k != 'BL2SEQ_LEGACY'}
env.update(LC_ALL='C', NODE_NO_WARNINGS='1', LOSAT_WASI_THREAD_CAP='8', LOSAT_WASM_MEMORY_MAXIMUM_PAGES='16384')
versions = ['baselinebody', a.candidate + 'body']
policy = dict(versions=versions, body_boundary='after CLI parse, immediately before Rust dispatch through preparation/search/output return', body_instrumentation='two clock reads and one post-boundary stderr print', warmup_per_version=1, pairs=3, order=['AB', 'BA', 'AB'], threads=a.threads, kinds=a.kinds, cases=a.cases, flags=[], affinity=sorted(os.sched_getaffinity(0)), timeout_seconds=300, concurrency='exclusive', time_guard='max(5%, 50 ms) for body and whole process', rss_guard='max(10%, 16 MiB)', status='RUNNING')
if a.diagnostic_background:
    policy.update(concurrency=a.diagnostic_background, purpose='diagnostic only', acceptance_eligible=False)
# NCBI blastn_app.cpp:59-67: record the complete timing environment separately.
if os.environ.get('LOSAT_BENCHMARK_ENVIRONMENT_POLICY'):
    policy.update(concurrency='monitored normal-desktop background', environment_policy=os.environ['LOSAT_BENCHMARK_ENVIRONMENT_POLICY'], acceptance_requires_environment_pass=True)
(out / 'policy.json').write_text(json.dumps(policy, indent=2) + '\n')
records = []

def summarize():
    summary = []
    for key in a.cases:
        for kind in a.kinds:
            for n in a.threads:
                groups = {v: [r for r in records if r['case'] == key and r['kind'] == kind and r['threads'] == n and r['timed'] and r['version'] == v] for v in versions}
                if any(len(group) != 3 for group in groups.values()):
                    continue
                z = dict(case=key, kind=kind, threads=n)
                for metric in ['wall_seconds', 'body_seconds', 'peak_rss_bytes']:
                    av, bv = ([r[metric] for r in groups[v]] for v in versions)
                    aa, bb = statistics.median(av), statistics.median(bv)
                    limit=max(aa*0.10,16*1024*1024) if metric=='peak_rss_bytes' else max(aa*0.05,0.05)
                    z[metric] = dict(baseline_samples=av, candidate_samples=bv, baseline_median=aa, candidate_median=bb, candidate_over_baseline=bb / aa, saved=aa - bb, regression_limit=limit, pass_guard=bb-aa<=limit)
                summary.append(z)
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    return summary

for key in a.cases:
    row = next(x for x in rows if x['losat_stem'] == key)
    common = ['-query', str(tests / 'fasta' / row['query']), '-subject', str(tests / 'fasta' / row['subject']), '-outfmt', '6', '-query_gencode', row['query_gencode'], '-db_gencode', row['db_gencode'], '-out', '{output}']
    expected = (b / 'oracles' / key / 'output.txt').read_bytes()
    for kind in a.kinds:
        for n in a.threads:
            for pair in [-1, 0, 1, 2]:
                ordered = versions[::-1] if pair == 1 else versions
                for v in ordered:
                    artifact = b / v / 'artifacts' / ('losat-' + kind + '-command' + ('.wasm' if kind != 'native' else ''))
                    if kind == 'native':
                        prefix = [str(artifact)]
                    else:
                        host = 'run_losat_wasi_threads.js' if kind == 'threaded' else 'run_losat_wasi.js'
                        prefix = [node, str(tests / host), str(artifact)]
                    d = out / key / f'{kind}-n{n}' / f'pair{pair}' / v
                    r = execute([*prefix, row['task'], *common, '-num_threads', str(n)], root, d, env, 300)
                    r.update(case=key, kind=kind, threads=n, pair=pair, timed=pair >= 0, version=v, artifact_sha256=digest(artifact))
                    r['raw_equal'] = r['status'] == 'PASS' and Path(r['output']).read_bytes() == expected
                    clocks = re.findall(r'\[BODY_SCOPE_SECONDS\] ([\d.]+)', (d / 'stderr.txt').read_text())
                    r['body_seconds'] = float(clocks[0]) if len(clocks) == 1 else None
                    if not r['raw_equal'] or r['body_seconds'] is None or r['body_seconds'] <= 0 or r['body_seconds'] > r['wall_seconds'] + 0.1:
                        r['status'] = 'GATE_FAIL'
                    records.append(r)
                    (out / 'runs.json').write_text(json.dumps(records, indent=2) + '\n')
                    summarize()
                    print(key, kind, n, pair, v, r['status'], round(r['wall_seconds'], 3), r['body_seconds'], flush=True)
                    if r['status'] != 'PASS':
                        policy['status'] = 'FAILED'
                        (out / 'policy.json').write_text(json.dumps(policy, indent=2) + '\n')
                        raise RuntimeError(r)
            condition=next(s for s in summarize() if s['case']==key and s['kind']==kind and s['threads']==n)
            if not all(condition[m]['pass_guard'] for m in ['wall_seconds','body_seconds','peak_rss_bytes']):
                policy['status']='REGRESSION_FAIL'
                (out/'policy.json').write_text(json.dumps(policy,indent=2)+'\n')
                raise SystemExit('Fixed control guard failed; no automatic extension')
policy['status'] = 'COMPLETE'
(out / 'policy.json').write_text(json.dumps(policy, indent=2) + '\n')
print('ALL PASS', len(records), flush=True)
