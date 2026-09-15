# NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:173-188
# (*thread)->Run(); (*thread)->Join(&result);
# One CPU partition, sequential A/B jobs after all untimed validation finishes.
import argparse,json,os,subprocess,time
from pathlib import Path
b=Path(__file__).resolve().parent;r=Path('/mnt/c/Users/genom/GitHub/LOSAT');t=r/'LOSAT/tests'
p=argparse.ArgumentParser();p.add_argument('--phase',choices=['cold','body','reuse'],required=True);a=p.parse_args()
os.sched_setaffinity(0,set(range(8)))
primary=['EDL933.Sakai.losatn.megablast','NZ_CP006932.NZ_CP006932.losatn.megablast','MelaMJNV.PemoMJNVA.losatn.blastn','MjPMNV.MlPMNV.losatn.blastn','MjeNMV.MelaMJNV.tlosatx','AP027280.AP027280.tlosatx','AP027078.AP027131.losatp','AP027132.NZ_CP006932.losatp']
controls=['Sakai.MG1655.losatn.megablast','NZ_CP006932.NZ_CP006932.losatn.blastn','MelaMJNV.PemoMJNVA.tlosatx','SicyWSV.CoBV.losatp','AP027131.NZ_CP006932.losatp','blastn-long-query-no-extension']
jobs=[]
groups=[primary+['single-query-32-matches'],controls] if a.phase=='cold' else ([primary+['single-query-32-matches']] if a.phase=='body' else [[c] for c in primary+['single-query-32-matches']])
for index,cases in enumerate(groups):
 candidate='integratedbody' if a.phase=='body' else 'integrated';baseline='C0body' if a.phase=='body' else 'C0'
 name='integrated-'+a.phase+(('-'+str(index)) if a.phase=='reuse' else ('-controls' if a.phase=='cold' and index==1 else ''))
 cmd=['python3',str(b/('body_benchmark.py' if a.phase=='body' else ('reuse_benchmark.py' if a.phase=='reuse' else 'benchmark_driver.py'))),'--candidate-dir',str(b/candidate),'--artifacts',str(b/candidate/'artifacts'),'--baseline-dir',str(b/baseline),'--baseline-runners',str(t),'--candidate-runners',str(t),'--oracle-dir','/home/kawato/micromamba/bin','--output-dir',str(b/name),'--candidate-source',str(b/(candidate+'-source')/'LOSAT'),'--snapshot',str(r/'docs/evidence/wasm_four_program_20260915/run-01/baseline-inputs.json'),'--kinds',*(['native','serial','threaded'] if a.phase=='cold' else (['serial','threaded'] if a.phase=='reuse' else ['native','threaded'])),'--threads','1','8','--warmups','1','--repeats','5']
 cmd+=['--skip-cold','--reuse-sessions','2','--reuse-timeout','3600'] if a.phase=='reuse' else ['--skip-reuse']
 # The plan's n1/n8 ratios and explicit serial timings concern primary fixtures.
 # Controls use the declared primary adoption condition, plus its native pair;
 # all controls retain the full n1/n2/n4/n8 raw-correctness matrix.
 if a.phase=='cold' and index==1:
  pos=cmd.index('--kinds');end=cmd.index('--warmups');cmd[pos:end]=['--kinds','native','threaded','--threads','8']
 for c in cases:cmd+=['--case',c]
 jobs.append({'name':name,'argv':cmd,'cases':cases})
(b/('final-'+a.phase+'-declaration.json')).write_text(json.dumps({'declared_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'jobs':jobs,'affinity':list(range(8)),'runtime':'Node26.8.2 default flags','reuse_boundary':'separate processes per case; same instance per case in reactor; one warmup + five repeats per n/session, AB then BA','role':'diagnostic body scope only' if a.phase=='body' else 'production cold/reuse measurement'},indent=2)+'\n')
assert json.loads((b/'validation-review.json').read_text())['allow_scoped_measurements'] is True
results=[]
for job in jobs:
 print('START',job['name'],job['cases'],flush=True)
 with (b/(job['name']+'.log')).open('w') as log:ret=subprocess.run(job['argv'],stdout=log,stderr=subprocess.STDOUT)
 results.append({**job,'returncode':ret.returncode});(b/('final-'+a.phase+'-results.json')).write_text(json.dumps(results,indent=2)+'\n')
 print('DONE',job['name'],ret.returncode,flush=True)
 # Preserve a failed fixed-window reuse series and still run independent cases.
 if ret.returncode and a.phase!='reuse':raise SystemExit(ret.returncode)

# NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:178-188
# (*thread)->Join(&result);
# Independent cases are all collected, but any failed case fails the phase.
if any(row['returncode'] for row in results):
 raise SystemExit(1)
