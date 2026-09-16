# NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:172-188
# (*thread)->Run(); (*thread)->Join(&result);
# Complete each real search and worker lifecycle before the next invocation.
# NCBI tabular.cpp:1098-1108: x_PrintField(*iter); m_Ostream << "\n";
# Only successful, raw-identical searches enter the fixed reuse comparison.
from pathlib import Path
import argparse, csv, json, os, shutil, statistics, sys
e = Path(__file__).resolve().parent
b = e / 'work'
root = e.parents[3]
tests = b / 'baseline-source/LOSAT/tests'
sys.path.insert(0, str(tests))
from wasm_performance import execute, digest
p = argparse.ArgumentParser()
p.add_argument('--label', required=True)
p.add_argument('--candidate', default='I1')
p.add_argument('--mode', choices=['compiled-module', 'same-instance'], default='compiled-module')
p.add_argument('--threads', type=int, nargs='+', default=[8])
p.add_argument('--cases', nargs='+', default=['MjeNMV.MelaMJNV.tlosatx', 'AP027280.AP027280.tlosatx'])
a = p.parse_args()
out = e / a.label
out.mkdir()
kind = 'threaded-command' if a.mode == 'compiled-module' else 'threaded-reactor'
versions = ['baseline', a.candidate]
policy = dict(versions=versions, mode=a.mode, kind=kind, cases=a.cases, threads=a.threads, sessions=['AB', 'BA'], warmup=1, samples=3, timeout_seconds=600, flags=[], concurrency='exclusive', affinity=sorted(os.sched_getaffinity(0)), status='RUNNING')
# NCBI blastn_app.cpp:59-67: record the complete timing environment separately.
if os.environ.get('LOSAT_BENCHMARK_ENVIRONMENT_POLICY'):
    policy.update(concurrency='monitored normal-desktop background', environment_policy=os.environ['LOSAT_BENCHMARK_ENVIRONMENT_POLICY'], acceptance_requires_environment_pass=True)
(out / 'policy.json').write_text(json.dumps(policy, indent=2) + '\n')
rows = list(csv.DictReader((tests / 'comparison_cases.tsv').open(), delimiter='\t'))
env = {k:v for k,v in os.environ.items() if not k.startswith(('LOSAT_', 'NODE_', 'RAYON_')) and k != 'BL2SEQ_LEGACY'}
env.update(LC_ALL='C', NODE_NO_WARNINGS='1', LOSAT_WASI_THREAD_CAP='8', LOSAT_WASM_MEMORY_MAXIMUM_PAGES='16384')
records, summaries = [], []
for case in a.cases:
    row = next(r for r in rows if r['losat_stem'] == case)
    query, subject = (tests / 'fasta' / row[name] for name in ['query', 'subject'])
    oracle = b / 'oracles' / case / 'output.txt'
    expected = oracle.read_bytes()
    for n in a.threads:
        for session in [0, 1]:
            grouped = {}
            for version in versions if session == 0 else versions[::-1]:
                d = out / case / f'n{n}' / f'session{session}' / version
                d.mkdir(parents=True)
                artifact = b / version / 'artifacts' / ('losat-' + kind + '.wasm')
                extra = ['-query_gencode', row['query_gencode'], '-db_gencode', row['db_gencode'], '-num_threads', str(n)]
                jobs = []
                for repeat in [-1, 0, 1, 2]:
                    output = d / f'raw-{repeat}.txt'
                    jobs.append(dict(case_id=case, repeat=repeat, timed=repeat >= 0, threads=n, program='tblastx', format='6', query_file=str(query), subject_file=str(subject), output=str(output), expected_sha256=digest(oracle), extra=extra, argv=['tblastx', '-query', str(query), '-subject', str(subject), '-outfmt', '6', *extra, '-out', str(output)]))
                jobs_file = d / 'jobs.json'
                jobs_file.write_text(json.dumps(jobs, indent=2) + '\n')
                cmd = [str(Path(shutil.which('node')).resolve()), str(tests / 'benchmark_wasi_reuse.js'), str(artifact), kind, a.mode, str(jobs_file), '-out', '{output}']
                result = execute(cmd, root, d / 'process', env, 600)
                record = dict(case=case, threads=n, session=session, version=version, artifact_sha256=digest(artifact), process=result)
                records.append(record)
                (out / 'runs.json').write_text(json.dumps(records, indent=2) + '\n')
                assert result['status'] == 'PASS', result
                report = json.loads(Path(result['output']).read_text())
                record['report'] = report
                for sample in report['samples']:
                    assert sample['status'] == 'PASS' and sample['raw_equal'] and sample['thread_contract'] == 'PASS', sample
                    assert Path(sample['output']).read_bytes() == expected
                timed = [s for s in report['samples'] if s['timed']]
                assert len(timed) == 3
                grouped[version] = timed
                (out / 'runs.json').write_text(json.dumps(records, indent=2) + '\n')
                print(case, n, session, version, 'PASS', [round(s['wall_seconds'], 6) for s in timed], flush=True)
            summary = dict(case=case, threads=n, session=session)
            for metric in ['wall_seconds', 'process_lifetime_peak_rss_bytes', 'memory_bytes']:
                av, bv = ([s[metric] for s in grouped[v]] for v in versions)
                aa, bb = statistics.median(av), statistics.median(bv)
                limit = max(aa * 0.05, 0.05) if metric == 'wall_seconds' else max(aa * 0.10, 16 * 1024 * 1024)
                summary[metric] = dict(baseline_samples=av, candidate_samples=bv, baseline_median=aa, candidate_median=bb, candidate_over_baseline=bb / aa, regression_limit=limit, pass_guard=bb-aa <= limit)
            summaries.append(summary)
            (out / 'summary.json').write_text(json.dumps(summaries, indent=2) + '\n')
            if not all(summary[m]['pass_guard'] for m in ['wall_seconds', 'process_lifetime_peak_rss_bytes', 'memory_bytes']):
                policy['status'] = 'REGRESSION_FAIL'
                (out / 'policy.json').write_text(json.dumps(policy, indent=2) + '\n')
                raise SystemExit('Fixed reuse guard failed; no automatic extension')
policy['status'] = 'COMPLETE'
(out / 'policy.json').write_text(json.dumps(policy, indent=2) + '\n')
