# NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1052-1054
# opt.SetDbGeneticCode(...);
# Approved local-subject code exception: compare to the frozen NCBI DB oracle.
# These concurrent correctness runs are excluded from all speed evidence.
from pathlib import Path
import argparse,concurrent.futures,json,os,shutil,sys
e=Path(__file__).resolve().parent;b=e/'work';root=e.parents[3];tests=b/'baseline-source/LOSAT/tests'
sys.path.insert(0,str(tests));from wasm_performance import execute,digest,diagnostic_diff
p=argparse.ArgumentParser();p.add_argument('--candidate',default='I2b');p.add_argument('--label',required=True);a=p.parse_args()
out=e/a.label;out.mkdir();ref=e/'recovered-run01/long-code4';m=json.loads((ref/'manifest.json').read_text());result=json.loads((ref/'oracle/result.json').read_text());expected=ref/'oracle/output.txt'
archive=e.parent/'run-01/verification.tar.gz';archives=json.loads((archive.parent/'archives.json').read_text())
assert digest(archive)==next(r['sha256'] for r in archives if r['file']==archive.name)
assert result['status']=='PASS' and digest(expected)==result['raw_output_sha256']
query=tests/'fasta/AP027131.fasta';subject=tests/'fasta/AP027133.fasta'
assert digest(query)==m['query_sha256'] and digest(subject)==m['subject_sha256']
assert digest(Path('/home/kawato/micromamba/bin/tblastx'))==m['oracle_sha256']
jobs=[('native',1),('serial',1),('threaded',1),('threaded',8)]
policy=dict(candidate=a.candidate,jobs=jobs,concurrency=2,timeout_seconds=1800,oracle_origin='frozen round-01 fresh NCBI DB oracle, not a fresh round-03 run',oracle_manifest=m,archive_sha256=digest(archive),oracle_output_sha256=digest(expected),role='correctness only, no speed inference')
(out/'policy.json').write_text(json.dumps(policy,indent=2)+'\n')
env={k:v for k,v in os.environ.items() if not k.startswith(('LOSAT_','NODE_','RAYON_'))};env.update(LC_ALL='C',NODE_NO_WARNINGS='1',LOSAT_WASI_THREAD_CAP='8',LOSAT_WASM_MEMORY_MAXIMUM_PAGES='16384')
def run(job):
    kind,n=job;artifact=b/a.candidate/'artifacts'/('losat-'+kind+'-command'+('.wasm' if kind!='native' else ''))
    prefix=[str(artifact)] if kind=='native' else [str(Path(shutil.which('node')).resolve()),str(tests/('run_losat_wasi.js' if kind=='serial' else 'run_losat_wasi_threads.js')),str(artifact)]
    argv=[*prefix,'tblastx','-query',str(query),'-subject',str(subject),'-query_gencode','4','-db_gencode','4','-outfmt','6','-num_threads',str(n),'-out','{output}']
    d=out/(kind+'-n'+str(n));r=execute(argv,root,d,env,1800)
    r.update(kind=kind,threads=n,artifact_sha256=digest(artifact),raw_equal=r['status']=='PASS' and Path(r['output']).read_bytes()==expected.read_bytes())
    if not r['raw_equal']:
        r['status']='PARITY_FAIL' if r['status']=='PASS' else r['status']
        if Path(r['output']).exists():diagnostic_diff(expected,Path(r['output']),d)
    (d/'comparison.json').write_text(json.dumps(r,indent=2)+'\n');print(kind,n,r['status'],flush=True);return r
records=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    for future in concurrent.futures.as_completed([pool.submit(run,job) for job in jobs]):
        records.append(future.result());(out/'runs.json').write_text(json.dumps(records,indent=2)+'\n')
assert len(records)==4 and all(r['status']=='PASS' for r in records)
print('ALL PASS',len(records),flush=True)
