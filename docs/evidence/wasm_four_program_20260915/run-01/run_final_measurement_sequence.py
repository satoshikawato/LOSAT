# NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:177-188
# (*thread)->Join(&result);
# Start samples only after all validation jobs and the explicit record review.
import json,subprocess,time
from pathlib import Path
b=Path(__file__).resolve().parent
while True:
 try:
  if len(json.loads((b/'final-validation.json').read_text()))==12:break
 except (FileNotFoundError,json.JSONDecodeError):pass
 time.sleep(2)
subprocess.run(['python3',str(b/'review_validation.py')],check=True)
records=[]
for phase in ['cold','body','reuse']:
 command=['python3',str(b/'final_measurements.py'),'--phase',phase]
 print('START',phase,flush=True)
 with (b/('final-'+phase+'-queue.log')).open('w') as log:r=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
 records.append({'phase':phase,'argv':command,'returncode':r.returncode});(b/'final-measurement-sequence.json').write_text(json.dumps(records,indent=2)+'\n')
 print('DONE',phase,r.returncode,flush=True)
 if r.returncode:raise SystemExit(r.returncode)
