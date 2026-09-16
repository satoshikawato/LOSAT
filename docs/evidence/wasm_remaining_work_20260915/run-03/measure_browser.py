# NCBI reference: c++/src/app/blast/blastn_app.cpp:59-67
# m_StopWatch.Start(); m_StopWatch.Elapsed();
# Real consumer full-call browser comparison; fixed sessions, no inner clocks.
from pathlib import Path
import argparse,json,subprocess,statistics,time,os,hashlib
e=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--candidate',default='I2b');p.add_argument('--label',required=True)
p.add_argument('--cases',nargs='+',default=['MjeNMV.MelaMJNV.tlosatx','AP027280.AP027280.tlosatx'])
p.add_argument('--threads',nargs='+',type=int,default=[8,1]);a=p.parse_args();out=e/a.label;out.mkdir()
versions=['baseline',a.candidate];records=[];summaries=[]
consumer=json.loads((e/'browser-consumer-manifest.json').read_text());consumer_root=Path(consumer['snapshot_root'])
harness_hash=hashlib.sha256((e/'browser_smoke.py').read_bytes()).hexdigest()
for filename,digest in consumer['files'].items():assert hashlib.sha256((consumer_root/filename).read_bytes()).hexdigest()==digest
policy=dict(versions=versions,cases=a.cases,threads=a.threads,sessions=['AB','BA'],warmup=1,samples=3,timeout_seconds=600,concurrency='exclusive',affinity=sorted(os.sched_getaffinity(0)),flags=[],interval='actual gbdraw runLosatPairsParallel promise, including dispatch/input transfer/search/output and current host completion semantics',reuse='same page/module and normal consumer worker ownership within each session; fresh browser for each version/session',memory_metric='peak sampled sum of descendant RSS, shared pages per process, 50 ms observation; not Wasm linear-memory size',status='RUNNING')
# NCBI blastn_app.cpp:59-67: record the complete timing environment separately.
if os.environ.get('LOSAT_BENCHMARK_ENVIRONMENT_POLICY'):
    policy.update(concurrency='monitored normal-desktop background', environment_policy=os.environ['LOSAT_BENCHMARK_ENVIRONMENT_POLICY'], acceptance_requires_environment_pass=True)
(out/'policy.json').write_text(json.dumps(policy,indent=2)+'\n')
(out/'consumer-binding.json').write_text(json.dumps(dict(consumer=consumer,harness_sha256=harness_hash),indent=2)+'\n')
for case in a.cases:
    for n in a.threads:
        for session in [0,1]:
            grouped={}
            for version in versions if session==0 else versions[::-1]:
                label=a.label+'/'+case+'/n'+str(n)+'/session'+str(session)+'/'+version
                child=e/label;child.parent.mkdir(parents=True,exist_ok=True)
                argv=['python3',str(e/'browser_smoke.py'),'--version',version,'--label',label,'--cases',case,'--threads',str(n),'--serial-kind','reactor','--repeat','4','--timed','--consumer-root',str(consumer_root)]
                started=time.monotonic()
                # NCBI prelim_stage.cpp:173-188: Run(); Join(&result).
                # A timed-out test owns its entire browser process tree.
                process=subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=subprocess.PIPE,start_new_session=True)
                timed_out=False
                try:stdout,stderr=process.communicate(timeout=600)
                except subprocess.TimeoutExpired:
                    import psutil
                    timed_out=True
                    try:descendants=psutil.Process(process.pid).children(recursive=True)
                    except psutil.NoSuchProcess:descendants=[]
                    for descendant in reversed(descendants):
                        try:descendant.kill()
                        except psutil.NoSuchProcess:pass
                    process.kill();stdout,stderr=process.communicate()
                child.mkdir(exist_ok=True)
                (child/'process.stdout').write_bytes(stdout);(child/'process.stderr').write_bytes(stderr)
                record=dict(case=case,threads=n,session=session,version=version,argv=argv,returncode=process.returncode,timed_out=timed_out,monotonic_start=started,monotonic_end=time.monotonic());records.append(record)
                (out/'runs.json').write_text(json.dumps(records,indent=2)+'\n')
                assert not timed_out and process.returncode==0,(record,stderr[-2000:])
                report=json.loads((child/'manifest.json').read_text());assert report['status']=='PASS' and not report['blocked_external']
                assert report['harness_sha256']==harness_hash
                for filename,digest in report['served_sha256'].items():
                    path=Path(filename)
                    if path.is_relative_to(consumer_root):assert consumer['files'][str(path.relative_to(consumer_root))]==digest
                samples=[r for r in report['records'] if r['timed']];assert len(samples)==3 and all(r['raw_equal'] and r['status']=='PASS' for r in samples)
                assert all(r['wall_seconds']>0 and r['browser_tree_peak_rss_bytes']>0 for r in samples)
                grouped[version]=samples;record['browser']=report['browser'];record['artifact']=report['artifacts']['threaded'];record['samples']=samples
                (out/'runs.json').write_text(json.dumps(records,indent=2)+'\n')
                print(case,n,session,version,'PASS',[round(s['wall_seconds'],6) for s in samples],flush=True)
            summary=dict(case=case,threads=n,session=session)
            for metric in ['wall_seconds','browser_tree_peak_rss_bytes']:
                aa,bb=([r[metric] for r in grouped[v]] for v in versions);av,bv=statistics.median(aa),statistics.median(bb)
                limit=max(av*0.05,0.05) if metric=='wall_seconds' else max(av*0.10,16*1024*1024)
                summary[metric]=dict(baseline_samples=aa,candidate_samples=bb,baseline_median=av,candidate_median=bv,candidate_over_baseline=bv/av,regression_limit=limit,pass_guard=bv-av<=limit)
            summaries.append(summary);(out/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
            if not all(v['pass_guard'] for k,v in summary.items() if isinstance(v,dict)):
                policy['status']='REGRESSION_FAIL';(out/'policy.json').write_text(json.dumps(policy,indent=2)+'\n');raise SystemExit('Browser fixed guard failed; no extension')
policy['status']='COMPLETE';(out/'policy.json').write_text(json.dumps(policy,indent=2)+'\n')
