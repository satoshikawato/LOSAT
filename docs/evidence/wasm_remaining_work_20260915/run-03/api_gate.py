# NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
# CATCH_ALL(status); return status;
# blast_args.cpp:3152-3187: thread count constraint.
# Validate real ABI identities and serial reactor error/recovery boundaries.
from pathlib import Path
import json,subprocess,os,hashlib,argparse
p=argparse.ArgumentParser();p.add_argument('--candidate',default='I2b');a=p.parse_args()
versions=['baseline',a.candidate]
e=Path(__file__).resolve().parent;b=e/'work';tests=b/'baseline-source/LOSAT/tests';out=e/(a.candidate+'-api-gate');out.mkdir()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
env={k:v for k,v in os.environ.items() if not k.startswith(('LOSAT_','NODE_','RAYON_'))};env['NODE_NO_WARNINGS']='1'
records=[];identities={}
for version in versions:
    for kind in ['serial-command','threaded-command','serial-reactor','threaded-reactor']:
        artifact=b/version/'artifacts'/('losat-'+kind+'.wasm');argv=['node',str(tests/'wasi_artifact.js'),str(artifact),kind]
        result=subprocess.run(argv,capture_output=True,env=env,timeout=60)
        stem=version+'-'+kind;(out/(stem+'.json')).write_bytes(result.stdout);(out/(stem+'.stderr')).write_bytes(result.stderr)
        records.append(dict(argv=argv,returncode=result.returncode,artifact_sha256=sha(artifact)));assert result.returncode==0,result.stderr
        identity=json.loads(result.stdout);identities[(version,kind)]=identity
for kind in ['serial-command','threaded-command','serial-reactor','threaded-reactor']:
    aa,bb=(dict(identities[(v,kind)]) for v in versions);aa.pop('sha256');bb.pop('sha256');assert aa==bb,kind
fixtures=out/'fixtures';fixtures.mkdir()
for name,source in [('nuc1.fasta','nuc1.fasta'),('aa1.fasta','query.faa')]:
    (fixtures/name).write_bytes((b/'boundary-fixtures'/source).read_bytes())
for version in versions:
    argv=['node',str(tests/'check_wasi_api_limits.js'),str(b/version/'artifacts/losat-serial-reactor.wasm'),'serial-reactor',str(fixtures),str(out/(version+'-limits'))]
    result=subprocess.run(argv,capture_output=True,env=env,timeout=60)
    (out/(version+'.stdout')).write_bytes(result.stdout);(out/(version+'.stderr')).write_bytes(result.stderr)
    records.append(dict(argv=argv,returncode=result.returncode));assert result.returncode==0,result.stderr
for file in (out/'baseline-limits').iterdir():assert file.read_bytes()==(out/(a.candidate+'-limits')/file.name).read_bytes(),file.name
(out/'manifest.json').write_text(json.dumps(dict(status='PASS',runs=records,scope='4 ABI kinds, serial reactor success/error/success for 3 programs; no unlimited reuse or memory-fix claim'),indent=2)+'\n')
print('All ABI and serial API boundaries PASS',flush=True)
