# NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:173-188
# (*thread)->Run(); (*thread)->Join(&result);
# Identical pinned CPU set and runner path for every paired measurement.
import os,subprocess,json,time
from pathlib import Path
b=Path(__file__).resolve().parent;r=Path('/mnt/c/Users/genom/GitHub/LOSAT');t=r/'LOSAT/tests'
os.sched_setaffinity(0,set(range(8)))
n=['MelaMJNV.PemoMJNVA.losatn.blastn','MjPMNV.MlPMNV.losatn.blastn']
configs=[('N1',n),('N2',n),('N3',n+['NZ_CP006932.NZ_CP006932.losatn.megablast']),('M1',['EDL933.Sakai.losatn.megablast','NZ_CP006932.NZ_CP006932.losatn.megablast','Sakai.MG1655.losatn.megablast']),('P1',['AP027078.AP027131.losatp','AP027132.NZ_CP006932.losatp','SicyWSV.CoBV.losatp']),('P2',['single-query-32-matches','AP027078.AP027131.losatp']),('X2',['MjeNMV.MelaMJNV.tlosatx','AP027280.AP027280.tlosatx','MelaMJNV.PemoMJNVA.tlosatx'])]
records=[]
for version,cases in configs:
 cmd=['python3',str(b/'benchmark_driver.py'),'--candidate-dir',str(b/version),'--artifacts',str(b/version/'artifacts'),'--baseline-dir',str(b/('M0' if version=='M1' else 'baseline')),'--baseline-runners',str(t),'--candidate-runners',str(t),'--oracle-dir','/home/kawato/micromamba/bin','--output-dir',str(b/(version+'-performance')),'--candidate-source',str(b/(version+'-source')/'LOSAT'),'--snapshot',str(r/'docs/evidence/wasm_four_program_20260915/run-01/baseline-inputs.json'),'--kinds','native','threaded','--threads','8','--warmups','1','--repeats','5','--skip-reuse']
 for c in cases:cmd+=['--case',c]
 print('START',version,flush=True)
 with (b/(version+'-performance.log')).open('w') as log:ret=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
 records.append({'version':version,'argv':cmd,'returncode':ret.returncode,'affinity':sorted(os.sched_getaffinity(0))});(b/'performance-queue.json').write_text(json.dumps(records,indent=2)+'\n')
 print('DONE',version,ret.returncode,flush=True)
 if ret.returncode:raise SystemExit(ret.returncode)
