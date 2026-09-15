# NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1052-1054
# if (m_Target == eDatabase ... ) { opt.SetDbGeneticCode(...); }
# Fixed longer deadline for the expensive correctness oracle, never a timing sample.
from pathlib import Path
import sys,os,json
b=Path(__file__).resolve().parent;root=Path('/mnt/c/Users/genom/GitHub/LOSAT');sys.path.insert(0,str(root/'LOSAT/tests'));from wasm_performance import execute,digest
old=json.loads((b/'baseline-matrix/AP027131.AP027133.tlosatx/oracle/command.json').read_text());argv=old['ordered_argv'];argv[argv.index('-out')+1]='{output}'
os.sched_setaffinity(0,{31})
env={k:v for k,v in os.environ.items() if not k.startswith(('LOSAT_','RAYON_','NODE_'))};env.update(LC_ALL='C')
out=b/'long-code4-oracle-1800';out.mkdir(exist_ok=False)
(out/'policy.json').write_text(json.dumps({'reason':'first correctness-oracle attempt timed out at 300s; preserve it, retry once with fixed 1800s ceiling','adoption_timing':False,'affinity':[31],'oracle_sha256':digest(argv[0])},indent=2))
r=execute(argv,root,out/'oracle',env,1800);print(r['status'],r['wall_seconds'],flush=True)
