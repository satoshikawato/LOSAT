# NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:173-188
# (*thread)->Run(); (*thread)->Join(&result);
# Untimed independent correctness jobs; all must stop before adoption timings.
import concurrent.futures,hashlib,json,os,subprocess,time
from pathlib import Path
b=Path(__file__).resolve().parent;t=Path('/mnt/c/Users/genom/GitHub/LOSAT/LOSAT/tests')
while True:
 try:
  ready=True
  for v in ['C0','integrated']:
   rows=json.loads((b/v/'build.json').read_text())[-4:]
   expected=hashlib.sha256((b/(v+'-source/LOSAT/src/algorithm/blastn/blast_engine/run.rs')).read_bytes()).hexdigest()
   kinds={'native-command','serial-command','threaded-command','threaded-reactor'}
   ready=ready and len(rows)==4 and all(x['status']==0 and x['source_sha256']['src/algorithm/blastn/blast_engine/run.rs']==expected for x in rows) and {k for x in rows for k in kinds if '/'+k+'/' in x['path']}==kinds
  if ready:break
 except (FileNotFoundError,json.JSONDecodeError,KeyError):pass
 time.sleep(2)
jobs=[]
for v in ['C0','integrated']:
 jobs.append((v+'-edge',['python3',str(b/'edge_gate.py'),v]))
 jobs.append((v+'-runtime',['python3',str(t/'check_wasm_threading.py'),'--native',str(b/v/'native-command/release/LOSAT'),'--serial',str(b/v/'artifacts/losat-serial-command.wasm'),'--threaded',str(b/v/'artifacts/losat-threaded-command.wasm'),'--reactor',str(b/v/'artifacts/losat-threaded-reactor.wasm'),'--oracle-dir','/home/kawato/micromamba/bin','--output-dir',str(b/(v+'-runtime'))]))
 cases=['EDL933.Sakai.losatn.megablast','NZ_CP006932.NZ_CP006932.losatn.megablast','Sakai.MG1655.losatn.megablast','MelaMJNV.PemoMJNVA.losatn.blastn','MjPMNV.MlPMNV.losatn.blastn','NZ_CP006932.NZ_CP006932.losatn.blastn']
 if v=='integrated':cases+=['MjeNMV.MelaMJNV.tlosatx','AP027280.AP027280.tlosatx','MelaMJNV.PemoMJNVA.tlosatx','AP027078.AP027131.losatp','AP027132.NZ_CP006932.losatp','SicyWSV.CoBV.losatp','AP027131.NZ_CP006932.losatp']
 jobs.append((v+'-matrix',['python3',str(b/'gate.py'),v,'--label','matrix','--case',*cases]))
for kind in ['native','serial','threaded']:
 for v in ['C0','integrated']:
  jobs.append((v+'-long-'+kind,['python3',str(b/'long_gate.py'),v,'--kinds',kind,'--label','long-'+kind]))
records=[]
def work(item):
 index,(name,command)=item
 # Assign a fixed partition to each executor thread through thread-local state.
 import threading
 slot=int(threading.current_thread().name.rsplit('_',1)[1]);aff=list(range(slot*8,(slot+1)*8))
 print('START',name,aff,flush=True)
 with (b/(name+'.log')).open('w') as log:
  ret=subprocess.run(['taskset','-c',','.join(map(str,aff)),*command],stdout=log,stderr=subprocess.STDOUT)
 code=ret.returncode
 if (b/name/'runs.json').exists() and any(r['status'] not in ['PASS','EXPECTED_UNSUPPORTED'] for r in json.loads((b/name/'runs.json').read_text())):code=code or 1
 print('DONE',name,code,flush=True)
 return {'name':name,'argv':command,'affinity':aff,'returncode':code,'role':'untimed correctness; not performance evidence'}
with concurrent.futures.ThreadPoolExecutor(max_workers=4,thread_name_prefix='gate') as pool:
 for future in concurrent.futures.as_completed([pool.submit(work,item) for item in enumerate(jobs)]):
  records.append(future.result());(b/'final-validation.json').write_text(json.dumps(records,indent=2)+'\n')
print('ALL JOBS FINISHED',flush=True)
