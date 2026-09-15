# NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:173-188
# (*thread)->Run(); (*thread)->Join(&result);
# This secondary cold condition changes only identical host flags in both versions.
import json
import os
import subprocess
import time
from pathlib import Path

b=Path(__file__).resolve().parent
r=Path('/mnt/c/Users/genom/GitHub/LOSAT')
t=r/'LOSAT/tests'
primary=['MjeNMV.MelaMJNV.tlosatx','AP027280.AP027280.tlosatx']
control=['MelaMJNV.PemoMJNVA.tlosatx']
jobs=[]
for name,baseline,candidate,cases,kinds,threads in [
    ('X2-turbofan-performance','baseline','X2',primary+control,['threaded'],[8]),
    ('integrated-turbofan-cold','C0','integrated',primary,['native','serial','threaded'],[1,8]),
    ('integrated-turbofan-controls','C0','integrated',control,['native','threaded'],[8]),
]:
    argv=['python3',str(b/'benchmark_driver.py'),'--candidate-dir',str(b/candidate),
          '--artifacts',str(b/candidate/'artifacts'),'--baseline-dir',str(b/baseline),
          '--baseline-runners',str(t),'--candidate-runners',str(t),
          '--oracle-dir','/home/kawato/micromamba/bin','--output-dir',str(b/name),
          '--candidate-source',str(b/(candidate+'-source')/'LOSAT'),
          '--snapshot',str(r/'docs/evidence/wasm_four_program_20260915/run-01/baseline-inputs.json'),
          '--node-arg=--no-liftoff','--node-arg=--no-wasm-tier-up',
          '--kinds',*kinds,'--threads',*[str(n) for n in threads],
          '--warmups','1','--repeats','5','--skip-reuse']
    for case in cases:
        argv += ['--case',case]
    jobs.append(dict(name=name,baseline=baseline,candidate=candidate,cases=cases,argv=argv))
declaration=dict(declared_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
                 role='secondary TurboFan-only cold condition; default Node remains primary adoption runtime',
                 flags=['--no-liftoff','--no-wasm-tier-up'], jobs=jobs,
                 scope='cold only, matching the previous explicitly limited host-flag adoption; no TurboFan body/reuse or browser claim',
                 affinity=list(range(8)), overlap='none; starts after default cold/body/reuse sequence finishes')
(b/'final-turbofan-declaration.json').write_text(json.dumps(declaration,indent=2)+'\n')
while True:
    try:
        if len(json.loads((b/'final-measurement-sequence.json').read_text()))==3:
            break
    except (FileNotFoundError,json.JSONDecodeError):
        pass
    time.sleep(2)
os.sched_setaffinity(0,set(range(8)))
records=[]
for job in jobs:
    print('START',job['name'],flush=True)
    with (b/(job['name']+'.log')).open('w') as log:
        result=subprocess.run(job['argv'],stdout=log,stderr=subprocess.STDOUT)
    records.append({**job,'returncode':result.returncode})
    (b/'final-turbofan-results.json').write_text(json.dumps(records,indent=2)+'\n')
    print('DONE',job['name'],result.returncode,flush=True)
raise SystemExit(int(any(row['returncode'] for row in records)))
