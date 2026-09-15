# NCBI references: core/na_ungapped.c:262-349; core/greedy_align.c:71-125;
# core/blast_kappa.c:2228-2236; core/link_hsps.c:827-895.
# Uint1 q_byte = (q[0]<<6)|(q[1]<<4)|(q[2]<<2)|q[3];
# Preserve supported results and existing explicit unsupported-option errors.
import argparse,json,os,shutil,sys
from pathlib import Path
b=Path(__file__).resolve().parent;r=Path('/mnt/c/Users/genom/GitHub/LOSAT');t=r/'LOSAT/tests';sys.path.insert(0,str(t));from wasm_performance import execute,digest,diagnostic_diff,validate_thread_evidence
from check_wasm_threading import fixtures
p=argparse.ArgumentParser();p.add_argument('version');p.add_argument('--threads',type=int,nargs='+',default=[1,2,4,8]);a=p.parse_args();out=b/(a.version+'-edge-tail');out.mkdir(exist_ok=False);f=fixtures(out/'fixtures')
q=out/'fixtures/nohit-q.fasta';q.write_text('>q\n'+'C'*240+'\n');s=out/'fixtures/nohit-s.fasta';s.write_text('>s\n'+'A'*240+'\n')
pa=out/'fixtures/nohit-q.faa';pa.write_text('>q\n'+'W'*85+'\n');pb=out/'fixtures/nohit-s.faa';pb.write_text('>s\n'+'A'*85+'\n')
seq=''.join(f['nuc1'].read_text().splitlines()[1:]);rev=out/'fixtures/reverse.fasta';rev.write_text('>q\n'+seq[:600].translate(str.maketrans('ACGT','TGCA'))[::-1]+'\n')
cases=[('blastn-word7','blastn',f['nuc3'],f['unequal'],['-task','blastn','-word_size','7']),('blastn-word11','blastn',f['nuc3'],f['unequal'],['-task','blastn','-word_size','11']),('blastn-compact','blastn',t/'fasta/blastn_parity_compact.fasta',t/'fasta/blastn_parity_compact.fasta',['-task','blastn']),('megablast-affine','blastn',f['nuc1'],f['unequal'],['-task','megablast','-gapopen','2','-gapextend','1']),('megablast-short','blastn',f['short-context'],f['short-subject'],['-task','megablast']),('blastp-single-unequal','blastp',b/'single-query/query.faa',b/'single-query/subjects.faa',[])]
for task in ['blastn','megablast']:
 cases.extend([(task+'-nohit','blastn',q,s,['-task',task]),(task+'-minus','blastn',rev,f['nuc1'],['-task',task,'-strand','minus'])])
null=out/'fixtures/nohit-ambiguous.fasta';null.write_text('>subject\n'+'N'*512+'\n')
cases += [('tblastx-nohit','tblastx',f['nuc1'],null,[]),('blastp-nohit','blastp',f['aa1'],pb,[])]
for ev in ['10','100','10000']:cases.append(('tblastx-lc-e'+ev,'tblastx',t/'fasta/LC738874.fasta',t/'fasta/LC738875.fasta',['-evalue',ev]))
cases=cases[next(i for i,c in enumerate(cases) if c[0]=='tblastx-nohit'):]
oracle=Path('/home/kawato/micromamba/bin');node=str(Path(shutil.which('node')).resolve());env={k:v for k,v in os.environ.items() if not k.startswith(('LOSAT_','RAYON_','NODE_')) and k!='BL2SEQ_LEGACY'};env.update(LC_ALL='C',NODE_NO_WARNINGS='1');records=[]
def save(): (out/'runs.json').write_text(json.dumps(records,indent=2)+'\n')
def paths():
 for kind in ['native','serial','threaded']:
  target={'serial':'wasm32-wasip1','threaded':'wasm32-wasip1-threads'}.get(kind);artifact=b/a.version/(kind+'-command')/(target or '')/'release'/('LOSAT.wasm' if target else 'LOSAT');prefix=[node,str(t/('run_losat_wasi_threads.js' if kind=='threaded' else 'run_losat_wasi.js')),str(artifact)] if target else [str(artifact)]
  yield kind,artifact,prefix
(out/'metadata.json').write_text(json.dumps({'artifacts':{k:{'path':str(x),'sha256':digest(x)} for k,x,p in paths()},'fixtures':{str(p):digest(p) for _,_,q,s,_ in cases for p in [q,s]},'node':node,'node_sha256':digest(node),'oracles':{p:{'path':str(oracle/p),'sha256':digest(oracle/p)} for p in ['blastn','blastp','tblastx']}},indent=2)+'\n')
for name,program,query,subject,extra in cases:
 args=['-query',str(query),'-subject',str(subject),'-outfmt','6',*extra,'-out','{output}']
 ref=execute([str(oracle/program),*args,'-num_threads','1'],r,out/name/'oracle',env,300);ref.update(case=name,kind='oracle');records.append(ref);save();assert ref['status']=='PASS',ref
 for kind,artifact,prefix in paths():
  for n in ([1] if kind=='serial' else a.threads):
   d=out/name/(kind+'-n'+str(n));row=execute([*prefix,program,*args,'-num_threads',str(n)],r,d,{**env,'LOSAT_WASI_THREADS_DEBUG':'1'},300);row.update(case=name,kind=kind,threads=n,raw_equal=row['status']=='PASS' and Path(row['output']).read_bytes()==Path(ref['output']).read_bytes())
   if row['status']=='PASS':
    if row['raw_equal']:row['thread_contract']=validate_thread_evidence((d/'stderr.txt').read_text(),n,kind)
    else:row['status']='PARITY_FAIL';diagnostic_diff(Path(ref['output']),Path(row['output']),d)
   records.append(row);save();assert row['status']=='PASS',row
 print(name,'PASS',flush=True)
for name,extra in [('matrix',['-matrix','BLOSUM45']),('compo0',['-comp_based_stats','0']),('compo1',['-comp_based_stats','1']),('compo3',['-comp_based_stats','3'])]:
 for kind,artifact,prefix in paths():
  n=1 if kind=='serial' else 8;d=out/('unsupported-'+name)/kind
  row=execute([*prefix,'blastp','-query',str(f['aa1']),'-subject',str(f['aa3']),'-outfmt','6',*extra,'-num_threads',str(n),'-out','{output}'],r,d,{**env,'LOSAT_WASI_THREADS_DEBUG':'1'},60);log=(d/'stderr.txt').read_text();ok=row['exit_status'] not in [None,0] and not row.get('raw_output_bytes',0) and ('unsupported' in log.lower() or 'not implemented' in log.lower()) and 'spawn_attempt' not in log
  row.update(case='unsupported-'+name,kind=kind,threads=n,status='EXPECTED_UNSUPPORTED' if ok else 'FAIL');records.append(row);save();assert ok,row
print('DONE',len(records),flush=True)
