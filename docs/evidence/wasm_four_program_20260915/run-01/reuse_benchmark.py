#!/usr/bin/env python3
"""Correct-output cold-process and explicit module/instance reuse measurements."""
import argparse
import csv
import shutil
import sys
import time
import hashlib
import math
import json
import os
import platform
from pathlib import Path
import statistics
import subprocess
sys.path.insert(0, '/mnt/c/Users/genom/GitHub/LOSAT/LOSAT/tests')
from check_wasm_threading import fixtures
from wasm_performance import execute, digest, validate_thread_evidence, diagnostic_diff

ROOT = Path('/mnt/c/Users/genom/GitHub/LOSAT')
TESTS = ROOT / 'LOSAT/tests'

# NCBI reference: c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-180
# m_iTotalLength += (Int8) (*iter)->length;
# Fixture shape changes work distribution, never the search algorithm.
def make_cases(out):
    inputs = fixtures(out)
    def sequence(path): return ''.join(path.read_text().split('>',1)[1].splitlines()[1:])
    lc1 = sequence(TESTS/'fasta/LC738874.fasta'); lc2 = sequence(TESTS/'fasta/LC738875.fasta')
    (out/'lc-query.fasta').write_text('>q\n'+lc1[:900]+'\n')
    # NCBI reference: c++/include/algo/blast/core/blast_gapalign.h:54
    # #define MAX_DBSEQ_LEN 5000000
    # Cross the real subject chunk boundary, rather than only the old pool threshold.
    long_sequence = sequence(TESTS/'fasta/EDL933.fna')
    assert len(long_sequence) > 5_000_000
    (out/'long-query.fasta').write_text('>q\n'+long_sequence[:900]+'\n')
    (out/'long-single.fasta').write_text('>s0\n'+long_sequence+'\n')
    (out/'lc-two.fasta').write_text('>s0\n'+lc1+'\n>s1\n'+lc2+'\n')
    proteins = [''.join(record.splitlines()[1:]) for record in (TESTS/'fasta/PajaWSV.faa').read_text().split('>')[1:]]
    aa = max(proteins,key=len)[:1600]
    (out/'protein-query.fasta').write_text('>q\n'+aa[:300]+'\n')
    (out/'protein-biased.fasta').write_text(''.join(f'>s{i}\n{aa[:length]}\n' for i,length in enumerate([1600,40,1600,80,1600,120])))
    return [
        ('small','tblastx',inputs['nuc1'],inputs['nuc1'],[]),
        ('multi-subject','tblastx',inputs['nuc3'],inputs['nuc3'],[]),
        ('old-threshold','blastn',out/'lc-query.fasta',out/'lc-two.fasta',['-task','megablast']),
        ('long-single','blastn',out/'long-query.fasta',out/'long-single.fasta',['-task','megablast']),
        ('dense-biased','blastp',out/'protein-query.fasta',out/'protein-biased.fasta',[]),
    ]

# NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-48,75
# const string kArgQuery("query"); const string kArgSubject("subject");
# const string kArgOutput("out"); const string kArgNumThreads("num_threads");
# Select existing comparison inputs without changing search options or fixtures.
def select_cases(out, selection):
    if not selection:
        return make_cases(out / "fixtures")
    rows = list(csv.DictReader((TESTS / "comparison_cases.tsv").open(), delimiter="\t"))
    cases = []
    for key in selection:
        # NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-48
        # const string kArgQuery("query"); const string kArgSubject("subject");
        # Use the same predeclared P2 input paths as the production cold harness.
        if key == "single-query-32-matches":
            fixture = Path(__file__).resolve().parent / "single-query"
            cases.append((key, "blastp", fixture / "query.faa", fixture / "subjects.faa", []))
            continue
        matches = [r for r in rows if key == r["losat_stem"]]
        if len(matches) != 1:
            raise ValueError(f"expected one comparison case for {key!r}")
        row = matches[0]; task = row["task"]
        extra = ["-task", task] if task in {"blastn", "megablast"} else []
        if task == "tblastx":
            extra = ["-query_gencode", row["query_gencode"], "-db_gencode", row["db_gencode"]]
        cases.append((key, "blastn" if task == "megablast" else task,
                      TESTS / "fasta" / row["query"], TESTS / "fasta" / row["subject"], extra))
    return cases


# NCBI reference: c++/src/objtools/align_format/tabular.cpp:1100-1108
# x_PrintField(*iter); m_Ostream << "\\n";
# Every admitted sample has matching raw bytes and status. Failed and unfinished
# conditions remain separate from performance summaries; this is not certification.
def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['candidate-dir', 'artifacts', 'baseline-dir', 'baseline-runners', 'oracle-dir', 'output-dir']:
        p.add_argument('--' + name, type=Path, required=True)
    # NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:171-188
    # (*thread)->Run(); (*thread)->Join(&result);
    # Cold timing includes host/worker loading; allow equal runner filesystem paths.
    p.add_argument('--candidate-runners', type=Path, default=TESTS,
                   help='candidate JS runner directory; use the same path as baseline when host code is unchanged')
    p.add_argument('--node', default='node')
    p.add_argument('--node-arg', action='append', default=[], help='repeat as --node-arg=--flag; applies to both versions')
    p.add_argument('--candidate-node-arg', action='append', default=[], help='additional candidate-only Node flags (P1)')
    p.add_argument('--case', action='append', help='exact losat_stem from comparison_cases.tsv; omit for threading fixtures')
    p.add_argument('--threads', type=int, nargs='+', default=[1, 2, 4, 8])
    # NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:173-180
    # (*thread)->Run(); (*thread)->Join(&result);
    # Serial compatibility measurements require an explicit --kinds serial.
    p.add_argument('--kinds', nargs='+', choices=['native', 'serial', 'threaded'], default=['native', 'threaded'])
    p.add_argument('--warmups', type=int, default=1)
    p.add_argument('--repeats', type=int, default=5)
    p.add_argument('--timeout', type=float, default=300)
    p.add_argument('--reuse-timeout', type=float, default=3600)
    p.add_argument('--skip-reuse', action='store_true')
    p.add_argument('--skip-cold', action='store_true')
    p.add_argument('--reuse-sessions', type=int, default=2, help='alternate A/B then B/A; each session keeps its own warmup and repetitions')
    p.add_argument('--candidate-source', type=Path, default=ROOT/'LOSAT', help='crate directory used to build candidate artifacts')
    p.add_argument('--snapshot', type=Path, help='baseline snapshot/build manifest identity, not a source of expected bytes')
    args = p.parse_args()
    if args.warmups < 0 or args.repeats < 0 or min(args.threads) < 1 or any(not math.isfinite(v) or v <= 0 for v in (args.timeout, args.reuse_timeout)) or args.reuse_sessions < 1:
        p.error('invalid measurement counts or timeout')
    if len(set(args.threads)) != len(args.threads):
        p.error('duplicate thread counts')
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    try:
        run_benchmark(args, out)
    except BaseException as error:
        status_path = out / 'run-status.json'
        record = json.loads(status_path.read_text()) if status_path.is_file() else {}
        if record.get('status') not in {'PARTIAL', 'FAILED'}:
            record['status'] = 'INTERRUPTED' if isinstance(error, KeyboardInterrupt) else 'FAILED'
        record['reason'] = str(error)
        record['ended_utc'] = time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
        status_path.write_text(json.dumps(record, indent=2) + '\n')
        raise


# NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:178-188
# (*thread)->Join(&result); if (result) { retv = reinterpret_cast<Uint8>(result); }
# Record failed or partial orchestration before propagating the failure to CI.
def run_benchmark(args, out):
    args.node = str(Path(shutil.which(args.node) or args.node).resolve())
    cases = select_cases(out, args.case)
    env = {k:v for k,v in os.environ.items() if not k.startswith(('LOSAT_', 'RAYON_', 'NODE_')) and k != 'BL2SEQ_LEGACY'}
    env.update(LC_ALL='C', NODE_NO_WARNINGS='1')
    nodes = {'baseline':[args.node, *args.node_arg], 'candidate':[args.node, *args.node_arg, *args.candidate_node_arg]}
    prefixes = {}
    # NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:171-188
    # (*thread)->Run(); (*thread)->Join(&result);
    # Select the actual host paths explicitly instead of mixing runner filesystems.
    for version, directory, runners in [('baseline', args.baseline_dir, args.baseline_runners), ('candidate', args.candidate_dir, args.candidate_runners)]:
        for kind in args.kinds:
            if kind == 'native':
                prefix = [str((directory / 'native-command/release/LOSAT').resolve())]
            else:
                target = 'wasm32-wasip1' + ('-threads' if kind == 'threaded' else '')
                artifact = (args.artifacts / f'losat-{kind}-command.wasm') if version == 'candidate' else directory / f'{kind}-command' / target / 'release/LOSAT.wasm'
                runner = runners / ('run_losat_wasi_threads.js' if kind == 'threaded' else 'run_losat_wasi.js')
                prefix = [*nodes[version], str(runner.resolve()), str(artifact.resolve())]
            prefixes[(version, kind)] = prefix
    samples = []; summaries = []; excluded = []; expected = {}; common_args = {}; valid = {}
    def save():
        for filename, value in [('samples.json', samples), ('summary.json', summaries), ('excluded.json', excluded)]:
            (out / filename).write_text(json.dumps(value, indent=2) + '\n')
    metadata = dict(run_id=out.name, started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()), node=json.loads(subprocess.check_output([args.node, '-p', 'JSON.stringify(process.versions)'], text=True)),
        node_argv=nodes, node_sha256=digest(args.node), warmups=args.warmups, repeats=args.repeats,
        purpose='production fixed-window reuse; independent failed sessions are retained and other sessions continue',
        timeout_seconds=args.timeout, reuse_timeout_seconds=args.reuse_timeout,
        head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        source_snapshot={'path':str(args.snapshot.resolve()), 'sha256':digest(args.snapshot)} if args.snapshot else None,
        candidate_source=str(args.candidate_source.resolve()),
        production_source_hashes={str(p.relative_to(args.candidate_source.parent)):digest(p) for p in sorted((args.candidate_source/'src').rglob('*.rs'))},
        cargo_manifest_sha256=digest(args.candidate_source/'Cargo.toml'),
        cargo_lock_sha256=digest(args.candidate_source/'Cargo.lock'),
        build_script_sha256=digest(args.candidate_source/'build.rs'),
        cargo_config_sha256=digest(args.candidate_source/'.cargo/config.toml'),
        hardware={'os':platform.platform(), 'filesystem':subprocess.check_output(['stat','-f','-c','%T',str(out)],text=True).strip(), 'cpuinfo':Path('/proc/cpuinfo').read_text(), 'meminfo':Path('/proc/meminfo').read_text(), 'affinity':sorted(os.sched_getaffinity(0)),
                  'cgroup_limits':{name:(Path('/sys/fs/cgroup')/name).read_text() for name in ['cpu.max', 'memory.max', 'cpuset.cpus.effective'] if (Path('/sys/fs/cgroup')/name).is_file()}},
        oracle_hashes={program:digest(args.oracle_dir/program) for program in ['blastn', 'blastp', 'tblastx']},
        oracle_versions={program:subprocess.check_output([str(args.oracle_dir/program), '-version'], text=True).strip() for program in ['blastn', 'blastp', 'tblastx']},
        baseline_runner_hashes={p.name:digest(p) for p in args.baseline_runners.glob('*.js')},
        # NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:171-188
        # (*thread)->Run(); (*thread)->Join(&result);
        # Equal content hashes alone do not prove equal cold host-loading costs.
        candidate_runner_hashes={p.name:digest(p) for p in args.candidate_runners.glob('*.js')},
        runner_paths={'baseline':str(args.baseline_runners.resolve()), 'candidate':str(args.candidate_runners.resolve())},
        environment={k:v for k,v in env.items() if k.startswith(('NODE_', 'LOSAT_', 'RAYON_')) or k == 'LC_ALL'},
        artifacts={f'{version}-{kind}':{'path':prefix[-1], 'sha256':digest(prefix[-1])} for (version,kind),prefix in prefixes.items()},
        fixtures={str(p):digest(p) for _,_,q,s,_ in cases for p in [q,s]}, runners={p.name:digest(p) for p in TESTS.iterdir() if p.suffix in {'.js', '.py'}},
        # NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:177-188
        # (*thread)->Run(); (*thread)->Join(&result);
        # Native invocation has no Node or Wasm compilation boundary.
        boundary='cold process starts at process launch and ends after process exit; see boundary_by_kind',
        boundary_by_kind={kind: ('native process startup, search, output, and exit' if kind == 'native' else
            'Node startup, artifact validation, guard, Wasm compilation, workers when requested, search, output, and teardown') for kind in args.kinds},
        linear_memory='not sampled in cold process; reuse records memory after invocation',
        reuse_host='common current benchmark_wasi_reuse.js and host imports for both versions; hashes in runners',
        elapsed_clock='CLOCK_MONOTONIC; realtime/GNU elapsed are adjustable-clock diagnostics',
        reuse_sessions=args.reuse_sessions, reuse_order='AB then BA, alternating per session')
    (out/'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    runner_snapshot = out/'runner-source'; runner_snapshot.mkdir()
    for name in metadata['runners']:
        shutil.copy2(TESTS/name, runner_snapshot/name)
    (out/'run-status.json').write_text(json.dumps({'status':'RUNNING'}))
    (out/'node-v8-options.txt').write_text(subprocess.check_output([args.node, '--v8-options'], text=True))
    # Flag acceptance on main and worker isolates; this diagnostic is outside timing.
    probe = "const {Worker}=require('node:worker_threads');console.log(JSON.stringify({main:process.execArgv}));new Worker(\"const{parentPort}=require('node:worker_threads');parentPort.postMessage(process.execArgv)\",{eval:true}).on('message',v=>console.log(JSON.stringify({worker:v})));"
    for version, node in nodes.items():
        (out/f'node-flags-{version}.txt').write_text(subprocess.check_output([*node, '-e', probe], text=True))
    for name, program, q, s, extra in cases:
        common = ['-query', str(q), '-subject', str(s), '-outfmt', '6', *extra, '-out', '{output}']; common_args[name] = common
        oracle_args = common
        # NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1052-1054
        # if (m_Target == eDatabase && args[kArgDbGeneticCode] &&
        #     (program == eTblastn || program == eTblastx)) { opt.SetDbGeneticCode(...); }
        # Retain the existing database oracle for the approved local-subject exception.
        if program == 'tblastx' and '-db_gencode' in extra and extra[extra.index('-db_gencode') + 1] != '1':
            db = out/'oracle-db'/name; db.parent.mkdir(exist_ok=True)
            command = [str(args.oracle_dir/'makeblastdb'), '-in', str(s), '-dbtype', 'nucl', '-parse_seqids', '-out', str(db)]
            result = subprocess.run(command, capture_output=True, env=env)
            (db.parent/f'{name}.json').write_text(json.dumps({
                'argv':command, 'exit_status':result.returncode,
                'tool_sha256':digest(command[0]), 'input_sha256':digest(s),
                'tool_version':subprocess.check_output([command[0], '-version'], text=True).strip(),
                'database_files':{str(f):digest(f) for f in sorted(db.parent.glob(name + '.*')) if f.suffix not in {'.log', '.json'}}}, indent=2))
            (db.parent/f'{name}.log').write_bytes(result.stdout + result.stderr)
            result.check_returncode()
            oracle_args = ['-query', str(q), '-db', str(db), '-outfmt', '6', *extra, '-out', '{output}']
        ref = execute([str((args.oracle_dir/program).resolve()), *oracle_args, '-num_threads', '1'], ROOT, out/'oracle'/name, env, args.timeout)
        ref.update(case_id=name, version='oracle', phase='oracle', timed=False); samples.append(ref); save()
        if ref['status'] != 'PASS':
            excluded.append(dict(case_id=name, status=ref['status'], reason='oracle failed; no timings')); save(); continue
        expected[name] = ref['raw_output_sha256']
        for (version, kind), prefix in prefixes.items():
            for n in [1] if kind == 'serial' else args.threads:
                label = f'{name}/{version}-{kind}-n{n}'; directory = out/'diagnostic'/label
                result = execute([*prefix, program, *common, '-num_threads', str(n)], ROOT, directory, {**env, 'LOSAT_WASI_THREADS_DEBUG':'1', 'LOSAT_TIMING':'1'}, args.timeout)
                result.update(case_id=name, version=version, kind=kind, threads=n, timed=False, phase='diagnostic',
                              raw_equal=result['status']=='PASS' and result['raw_output_sha256']==expected[name])
                if result['status'] == 'PASS' and not result['raw_equal']:
                    result['status'] = 'PARITY_FAIL'; diagnostic_diff(Path(ref['output']), Path(result['output']), directory)
                if result['status'] == 'PASS':
                    try: result['thread_contract'] = validate_thread_evidence((directory/'stderr.txt').read_text(), n, kind)
                    except Exception as error: result.update(status='THREAD_CONTRACT_FAIL', reason=str(error))
                valid[(name, version, kind, n)] = result['status'] == 'PASS'
                samples.append(result); save(); print(label, result['status'], flush=True)
    configurations = [(name, program, version, kind, n) for name, program, *_ in cases for version, kind in prefixes
                      for n in ([1] if kind == 'serial' else args.threads)]
    runnable = []
    for config in configurations:
        name, program, version, kind, n = config
        if all(valid.get((name, v, kind, n), False) for v in ['baseline', 'candidate']): runnable.append(config)
        else: excluded.append(dict(case_id=name, version=version, kind=kind, threads=n, status='NOT_RUN', reason='paired diagnostic gate failed'))
    save()
    for repeat in range(0 if args.skip_cold else args.warmups + args.repeats):
        order = runnable if repeat % 2 == 0 else list(reversed(runnable))
        for name, program, version, kind, n in order:
            directory = out/'cold'/f'repeat-{repeat}'/f'{name}-{version}-{kind}-n{n}'
            result = execute([*prefixes[(version,kind)], program, *common_args[name], '-num_threads', str(n)], ROOT, directory, env, args.timeout)
            result.update(case_id=name, version=version, kind=kind, threads=n, timed=repeat >= args.warmups, repeat=repeat, phase='cold',
                          raw_equal=result['status']=='PASS' and result['raw_output_sha256']==expected[name])
            if result['status'] == 'PASS' and not result['raw_equal']: result['status'] = 'PARITY_FAIL'
            # NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:178-188
            # (*thread)->Join(&result);
            # A clock disagreement invalidates timing, independently of raw parity.
            if result['status'] == 'PASS' and not result['monotonic_clock_agreement']:
                result.update(status='CLOCK_INCONSISTENT', reason='boottime disagrees with monotonic duration')
            samples.append(result); save()
            if result['status'] != 'PASS':
                # Stop this run after preserving the failed row; incomplete series never get a median.
                raise RuntimeError(f'{name}/{version}-{kind}-n{n}: {result["status"]}; evidence: {directory}')
        print('cold repetition', repeat, 'complete', flush=True)
    for name, program, version, kind, n in runnable:
        rows = [x for x in samples if x.get('timed') and x['phase']=='cold' and (x['case_id'],x['version'],x['kind'],x['threads'])==(name,version,kind,n)]
        if not rows: continue
        if len(rows) != args.repeats or any(x['status'] != 'PASS' for x in rows): raise RuntimeError('incomplete measurement series')
        values = [x['wall_seconds'] for x in rows]; rss = [x['peak_rss_bytes'] for x in rows if x['peak_rss_bytes'] is not None]
        summaries.append(dict(case_id=name, version=version, kind=kind, threads=n, boundary='cold-process',
                              median_seconds=statistics.median(values), min_seconds=min(values), max_seconds=max(values), samples=values,
                              peak_rss_bytes=max(rss) if rss else None, raw_output_sha256=expected[name]))
    save()
    if not args.skip_reuse:
        for session in range(args.reuse_sessions):
            for version in (['baseline', 'candidate'] if session % 2 == 0 else ['candidate', 'baseline']):
                for kind, mode in [('serial-command','compiled-module'), ('threaded-command','compiled-module'), ('threaded-reactor','same-instance')]:
                    base_kind = kind.split('-')[0]
                    if base_kind not in args.kinds: continue
                    jobs = []; directory = out/'reuse'/f'session{session}-{version}-{kind}-{mode}'; directory.mkdir(parents=True)
                    for repeat in range(args.warmups + args.repeats):
                        for name, program, q, s, extra in (cases if repeat % 2 == 0 else list(reversed(cases))):
                            for n in ([1] if base_kind=='serial' else (args.threads if repeat % 2 == 0 else list(reversed(args.threads)))):
                                if not all(valid.get((name,v,base_kind,n), False) for v in ['baseline','candidate']): continue
                                output = directory/f'{name}-repeat{repeat}-n{n}.out'
                                jobs.append(dict(case_id=name, program=program, query_file=str(q), subject_file=str(s), extra=[*extra,'-num_threads',str(n)],
                                    argv=[program,*[str(output) if x=='{output}' else x for x in common_args[name]],'-num_threads',str(n)],
                                    output=str(output), threads=n, repeat=repeat, timed=repeat>=args.warmups, expected_sha256=expected[name]))
                    if not jobs: continue
                    jobs_path=directory/'jobs.json'; jobs_path.write_text(json.dumps(jobs,indent=2)+'\n')
                    target='wasm32-wasip1' + ('-threads' if base_kind=='threaded' else '')
                    artifact = (args.artifacts/f'losat-{kind}.wasm') if version=='candidate' else args.baseline_dir/kind/target/'release/LOSAT.wasm'
                    command=[*nodes[version],str(TESTS/'benchmark_wasi_reuse.js'),str(artifact.resolve()),kind,mode,str(jobs_path),'-out','{output}']
                    result=execute(command, ROOT, directory/'process', env, args.reuse_timeout)
                    # NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:178-188
                    # (*thread)->Join(&result);
                    # Apply the same elapsed-clock contract to every process lifetime.
                    if result['status'] == 'PASS' and not result['monotonic_clock_agreement']:
                        result.update(status='CLOCK_INCONSISTENT', reason='boottime disagrees with monotonic duration')
                    result.update(version=version,kind=kind,session=session,phase='reuse-process',timed=False);samples.append(result);save()
                    # NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:178-188
                    # (*thread)->Join(&result);
                    # Retain failed finite sessions; independent version/mode/session
                    # processes still run. No failed session is retried or extended.
                    if result['status'] != 'PASS':
                        excluded.append(dict(version=version, kind=kind, session=session,
                            status=result['status'], reason='fixed-window reuse process failed',
                            evidence=str(directory)))
                        save()
                        continue
                    try:
                        data=json.loads(Path(result['output']).read_text())
                        complete=(len(data['samples']) == len(jobs) and all(x['status']=='PASS' for x in data['samples']))
                    except (OSError, ValueError, KeyError, TypeError) as error:
                        result['reuse_schema_error']=str(error)
                        complete=False
                    if not complete:
                        result['status']='INCOMPLETE_REUSE'
                        excluded.append(dict(version=version, kind=kind, session=session,
                            status='INCOMPLETE_REUSE', reason='incomplete fixed-window samples',
                            evidence=str(directory)))
                        save()
                        continue
                    print(version,kind,mode,'PASS',flush=True)
    status = 'COMPLETE' if not excluded else ('PARTIAL' if runnable else 'FAILED')
    (out/'run-status.json').write_text(json.dumps({'status':status, 'excluded':len(excluded), 'samples':len(samples), 'ended_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}, indent=2)+'\n')
    print('completed',len(samples),'process samples; reuse samples separate; excluded',len(excluded),flush=True)
    if status != 'COMPLETE': raise RuntimeError(f'{status}: see excluded.json')

if __name__=='__main__': main()
