# NCBI reference: c++/src/objtools/align_format/tabular.cpp:1098-1108
# x_PrintField(*iter); m_Ostream << "\n";
# Compare unmodified raw output; NCBI executables are validation oracles only.
import argparse,csv,hashlib,json,os,shutil,sys,subprocess,time
from pathlib import Path
base=Path(__file__).resolve().parent;root=Path('/mnt/c/Users/genom/GitHub/LOSAT');tests=root/'LOSAT/tests'
sys.path.insert(0,str(tests));from wasm_performance import execute,digest,diagnostic_diff,validate_thread_evidence
p=argparse.ArgumentParser();p.add_argument('version');p.add_argument('--kinds',nargs='+',default=['native','serial','threaded']);p.add_argument('--threads',nargs='+',type=int,default=[1,2,4,8]);p.add_argument('--case',nargs='+');p.add_argument('--label',default='gate');a=p.parse_args()
out=base/(a.version+'-'+a.label);out.mkdir(exist_ok=False)
rows=list(csv.DictReader((tests/'comparison_cases.tsv').open(),delimiter='\t'))
rows.append(dict(losat_stem='single-query-32-matches',task='blastp',query=str(base/'single-query/query.faa'),subject=str(base/'single-query/subjects.faa')))
selection=['EDL933.Sakai.losatn.megablast','NZ_CP006932.NZ_CP006932.losatn.megablast','Sakai.MG1655.losatn.megablast','MelaMJNV.PemoMJNVA.losatn.blastn','MjPMNV.MlPMNV.losatn.blastn','NZ_CP006932.NZ_CP006932.losatn.blastn','MjeNMV.MelaMJNV.tlosatx','AP027280.AP027280.tlosatx','MelaMJNV.PemoMJNVA.tlosatx','AP027131.AP027133.tlosatx','AP027078.AP027131.losatp','AP027132.NZ_CP006932.losatp','SicyWSV.CoBV.losatp','AP027131.NZ_CP006932.losatp']
selected=[next(r for r in rows if r['losat_stem']==name) for name in (a.case or selection)]
node=str(Path(shutil.which('node')).resolve());oracle=Path('/home/kawato/micromamba/bin');env={k:v for k,v in os.environ.items() if not k.startswith(('LOSAT_','RAYON_','NODE_')) and k!='BL2SEQ_LEGACY'};env.update(LC_ALL='C',NODE_NO_WARNINGS='1')
records=[]
meta={'node':node,'node_version':subprocess.check_output([node,'--version'],text=True),'node_sha256':digest(node),'oracles':{prog:{'path':str(oracle/prog),'sha256':digest(oracle/prog),'version':subprocess.check_output([str(oracle/prog),'-version'],text=True)} for prog in ['blastn','blastp','tblastx']},'role':'raw/source comparison, untimed; no release fingerprint updates','fixtures':[]}
for row in selected:
 name=row['losat_stem']; task=row['task'];program='blastn' if task=='megablast' else task;q=tests/'fasta'/row['query'];s=tests/'fasta'/row['subject'];extra=['-task',task] if program=='blastn' else (['-query_gencode',row['query_gencode'],'-db_gencode',row['db_gencode']] if program=='tblastx' else [])
 meta['fixtures'].append({'case':name,'query':str(q),'query_sha256':digest(q),'subject':str(s),'subject_sha256':digest(s),'extra':extra});(out/'metadata.json').write_text(json.dumps(meta,indent=2))
 common=['-query',str(q),'-subject',str(s),'-outfmt','6',*extra,'-out','{output}'];refargs=common
 if program=='tblastx' and row['db_gencode']!='1':
  db=out/'db'/name;db.parent.mkdir(exist_ok=True)
  dbcmd=[str(oracle/'makeblastdb'),'-in',str(s),'-dbtype','nucl','-parse_seqids','-out',str(db)]
  r=subprocess.run(dbcmd,capture_output=True);(db.parent/(name+'.log')).write_bytes(r.stdout+r.stderr);r.check_returncode();refargs=['-query',str(q),'-db',str(db),'-outfmt','6',*extra,'-out','{output}']
 ref=execute([str(oracle/program),*refargs,'-num_threads','1'],root,out/name/'oracle',env,300);ref.update(case=name,kind='oracle');records.append(ref)
 (out/'runs.json').write_text(json.dumps(records,indent=2)+'\n')
 if ref['status']!='PASS':
  print(name,'ORACLE',ref['status'],'dependent paths skipped',flush=True);continue
 for kind in a.kinds:
  target={'serial':'wasm32-wasip1','threaded':'wasm32-wasip1-threads'}.get(kind)
  artifact=base/a.version/(kind+'-command')/(target or '')/'release'/('LOSAT.wasm' if target else 'LOSAT')
  runner=tests/('run_losat_wasi_threads.js' if kind=='threaded' else 'run_losat_wasi.js')
  prefix=[node,str(runner),str(artifact)] if target else [str(artifact)]
  for n in ([1] if kind=='serial' else a.threads):
   directory=out/name/(kind+'-n'+str(n));r=execute([*prefix,program,*common,'-num_threads',str(n)],root,directory,{**env,'LOSAT_TIMING':'1','LOSAT_WASI_THREADS_DEBUG':'1'},300)
   r.update(case=name,kind=kind,threads=n,artifact_sha256=digest(artifact));r['raw_equal']=r['status']=='PASS' and ref['status']=='PASS' and Path(r['output']).read_bytes()==Path(ref['output']).read_bytes()
   if r['raw_equal']:
    try:r['thread_contract']=validate_thread_evidence((directory/'stderr.txt').read_text(),n,kind)
    except Exception as e:r.update(status='THREAD_CONTRACT_FAIL',reason=str(e))
   else:
    r['status']='PARITY_FAIL' if r['status']=='PASS' else r['status']
    if Path(r.get('output','/nonexistent')).exists() and Path(ref.get('output','/nonexistent')).exists():diagnostic_diff(Path(ref['output']),Path(r['output']),directory)
   records.append(r);(out/'runs.json').write_text(json.dumps(records,indent=2)+'\n');print(name,kind,n,r['status'],flush=True)
print('DONE',len(records),'failures',sum(r['status']!='PASS' for r in records),flush=True)
