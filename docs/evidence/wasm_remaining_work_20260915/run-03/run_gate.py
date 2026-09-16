# NCBI reference: c++/src/objtools/align_format/tabular.cpp:1098-1108
# x_PrintField(*iter); m_Ostream << "\n";
# Fresh oracle comparison; never updates PR5/PR6 expected bytes.
import argparse,csv,hashlib,json,os,shutil,sys
from pathlib import Path
evidence=Path(__file__).resolve().parent;b=evidence/'work';root=evidence.parents[3];tests=b/'baseline-source/LOSAT/tests'
sys.path.insert(0,str(tests));from wasm_performance import execute,digest,diagnostic_diff,validate_thread_evidence
p=argparse.ArgumentParser();p.add_argument('version');p.add_argument('--label',required=True);p.add_argument('--cases',nargs='+');p.add_argument('--kinds',nargs='+',default=['native','serial','threaded']);p.add_argument('--threads',type=int,nargs='+',default=[1,8]);p.add_argument('--diagnostic',action='store_true');p.add_argument('--flags',default='[]');p.add_argument('--fresh',action='store_true');a=p.parse_args()
out=b/a.label;out.mkdir(exist_ok=False);oracles=b/'oracles';oracles.mkdir(exist_ok=True)
rows=list(csv.DictReader((tests/'comparison_cases.tsv').open(),delimiter='\t'))
selection=['MjeNMV.MelaMJNV.tlosatx','AP027280.AP027280.tlosatx','MelaMJNV.PemoMJNVA.tlosatx','EDL933.Sakai.losatn.megablast','NZ_CP006932.NZ_CP006932.losatn.megablast','MelaMJNV.PemoMJNVA.losatn.blastn','MjPMNV.MlPMNV.losatn.blastn','AP027078.AP027131.losatp','AP027132.NZ_CP006932.losatp']
for e in (10,100,10000):rows.append(dict(losat_stem=f'threshold-{e}',task='tblastx',query='LC738874.fasta',subject='LC738875.fasta',query_gencode='1',db_gencode='1',extra=['-evalue',str(e)]))
selection += ['threshold-10','threshold-100','threshold-10000']
edge=b/'boundary-fixtures';edge.mkdir(exist_ok=True)
rows += [dict(losat_stem='single32',task='blastp',query=str(edge/'query.faa'),subject=str(edge/'subjects.faa')),dict(losat_stem='word7',task='blastn',query=str(edge/'nuc3.fasta'),subject=str(edge/'unequal.fasta'),extra=['-word_size','7']),dict(losat_stem='valid-nohit',task='tblastx',query=str(edge/'nuc1.fasta'),subject=str(edge/'nohit.fasta'),query_gencode='1',db_gencode='1')]
selection += ['single32','word7','valid-nohit']
selected=[next(r for r in rows if r['losat_stem']==name) for name in (a.cases or selection)]
node=str(Path(shutil.which('node')).resolve());oracle=Path('/home/kawato/micromamba/bin');env={k:v for k,v in os.environ.items() if not k.startswith(('LOSAT_','RAYON_','NODE_')) and k!='BL2SEQ_LEGACY'};env.update(LC_ALL='C',NODE_NO_WARNINGS='1',LOSAT_WASI_THREAD_CAP='8',LOSAT_WASM_MEMORY_MAXIMUM_PAGES='16384')
records=[];meta={'version':a.version,'flags':json.loads(a.flags),'affinity':sorted(os.sched_getaffinity(0)),'fixtures':[],'oracles':{p:{'path':str(oracle/p),'sha256':digest(oracle/p)} for p in ['blastn','blastp','tblastx']}}
for row in selected:
 name=row['losat_stem'];task=row['task'];prog='blastn' if task=='megablast' else task;q=tests/'fasta'/row['query'];s=tests/'fasta'/row['subject'];extra=(['-task',task] if prog=='blastn' else ['-query_gencode',row['query_gencode'],'-db_gencode',row['db_gencode']] if prog=='tblastx' else [])+row.get('extra',[])
 common=['-query',str(q),'-subject',str(s),'-outfmt','6',*extra,'-out','{output}']
 meta['fixtures'].append({'case':name,'program':prog,'query':str(q),'subject':str(s),'query_sha256':digest(q),'subject_sha256':digest(s),'common':common});(out/'manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
 refdir=oracles/name
 if a.fresh or not refdir.exists():
  ref=execute([str(oracle/prog),*common,'-num_threads','1'],root,refdir,env,300)
 else:ref=json.loads((refdir/'result.json').read_text())
 if ref['status']!='PASS':raise RuntimeError((name,'oracle',ref))
 for kind in a.kinds:
  target={'serial':'wasm32-wasip1','threaded':'wasm32-wasip1-threads'}.get(kind);art=b/a.version/'artifacts'/(f'losat-{kind}-command'+('.wasm' if target else ''));runner=tests/('run_losat_wasi_threads.js' if kind=='threaded' else 'run_losat_wasi.js');prefix=[node,*json.loads(a.flags),str(runner),str(art)] if target else [str(art)]
  for n in ([1] if kind=='serial' else a.threads):
   d=out/name/f'{kind}-n{n}';runenv={**env,**({'LOSAT_TIMING':'1','LOSAT_WASI_THREADS_DEBUG':'1'} if a.diagnostic else {})}
   r=execute([*prefix,prog,*common,'-num_threads',str(n)],root,d,runenv,300);r.update(case=name,kind=kind,threads=n,artifact_sha256=digest(art),oracle_output=str(refdir/'output.txt'));r['raw_equal']=r['status']=='PASS' and Path(r['output']).read_bytes()==(refdir/'output.txt').read_bytes()
   if not r['raw_equal']:
    r['status']='PARITY_FAIL' if r['status']=='PASS' else r['status']
    if Path(r['output']).exists():diagnostic_diff(refdir/'output.txt',Path(r['output']),d)
   elif a.diagnostic:
    try:r['thread_contract']=validate_thread_evidence((d/'stderr.txt').read_text(),n,kind)
    except Exception as e:r.update(status='THREAD_CONTRACT_FAIL',reason=str(e))
   records.append(r);(out/'runs.json').write_text(json.dumps(records,indent=2)+'\n');print(name,kind,n,r['status'],round(r['wall_seconds'],3),flush=True)
   if r['status']!='PASS':raise RuntimeError((name,kind,n,r['status']))
print('ALL PASS',len(records),flush=True)
