#!/usr/bin/env python3
"""Comparison fixtures only. Run with NEW output directory and clean NCBI archive.
NCBI source: c++/src/algo/blast/unit_tests/api/bl2seq_unit_test.cpp:1638-1654:
 CBl2Seq blaster(*query, *subj, eBlastx); TSeqAlignVector sav(blaster.Run());
NCBI source: blastfilter_unit_test.cpp:949-1035: masks.push_back(TSeqRange(0,75));
NCBI source: c++/src/objects/seqfeat/gc.prt: ncbieaa / Base1 / Base2 / Base3.
Synthetic inputs are deterministic probes, never hand-authored search expectations.
"""
import argparse,csv,hashlib,itertools,json,re,subprocess
from pathlib import Path
from Bio import SeqIO

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
CODES=[1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,33]
AA='MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQANLADKPEADKVKAKGDKVIGIDLGTTNSCVAIMDGTTPRVLENLEKKYP'
CODONS=dict(zip('ACDEFGHIKLMNPQRSTVWY','GCT TGT GAT GAA TTT GGT CAT ATT AAA CTG ATG AAT CCT CAA CGT TCT ACT GTT TGG TAT'.split()))
DNA=''.join(CODONS[a] for a in AA)
def rc(s):return s.translate(str.maketrans('ACGTacgt','TGCAtgca'))[::-1]
def fasta(records):return ''.join('>'+i+'\n'+s+'\n' for i,s in records).encode()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('output',type=Path);ap.add_argument('--ncbi-source',type=Path,required=True);a=ap.parse_args()
 a.output.mkdir();rows=[]; deriv=[]
 tools=json.loads((HERE/'authority/acquisition_tools.json').read_text())

 for line in (HERE/'authority/runtime_libraries.sha256').read_text().splitlines():
  checksum,path=line.split('  ',1);assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==checksum,path
 for name in ['blastx','blastdbcmd']:
  assert sha(Path(tools[name]['path']))==tools[name]['sha256'],name
 for line in (HERE/'ncbi_source.sha256').read_text().splitlines():
  checksum,path=line.split('  ',1);assert sha(a.ncbi_source/path)==checksum,path
 def put(name,raw,family,origin,tier='INITIAL_REQUIRED'):
  p=a.output/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(raw)
  rows.append(dict(fixture_id=name,family=family,tier=tier,path=name,bytes=len(raw),sha256=sha(p),origin=origin,oracle='PINNED_NCBI_CLI_REPORT; NCBI_SOURCE_INTERNAL',branch_status='PROBE_ONLY; later internal trace required'))
 def q(name,seq,fam):put(name+'.fna',fasta([(name,seq)]),fam,'generator v1 fixed peptide/codon map; no RNG')
 put('protein.faa',fasta([('protein '+AA[:12],AA)]),'S01','generator fixed peptide')
 for sign in [1,-1]:
  for frame in [1,2,3]:
   s='C'*(frame-1)+DNA
   q('S01_frame'+str(sign*frame),s if sign>0 else rc(s),'S01')
 for n in [0,1,2,3,8,9,10,14,15,16,237,238,239]:q('S02_len'+str(n),(DNA+'A'*300)[:n],'S02')
 gc=(a.ncbi_source/'c++/src/objects/seqfeat/gc.prt').read_text()
 tables={}
 for block in re.findall(r'\{[^{}]*\}',gc,re.S):
  m=re.search(r'\bid\s+(\d+)',block); aa=re.search(r'ncbieaa\s+"([^"]+)"',block)
  if not m or not aa:continue
  bases=[re.search(r'Base'+str(i)+r'\s+([TCAG]+)',block).group(1) for i in [1,2,3]]
  tables[int(m.group(1))]=dict(zip(map(''.join,zip(*bases)),aa.group(1)))
 codons=list(tables[1]);assert len(codons)==64 and set(CODES)<set(tables)
 # GC table is source-derived unit expected, not a search report.
 put('codon_tables.json',json.dumps(tables,indent=2,sort_keys=True).encode()+b'\n','S03','pinned NCBI gc.prt')
 for code in CODES:
  seq=DNA+''.join(codons)+DNA
  put('S03_code'+str(code)+'.fna',fasta([('S03_code'+str(code)+'_plus',seq),('S03_code'+str(code)+'_minus',rc(seq))]),'S03','all64codons selected-code probe on both strands')
  put('S03_code'+str(code)+'.faa',fasta([('code'+str(code),AA+''.join(tables[code][c] for c in codons)+AA)]),'S03','protein letters from pinned gc.prt; retains stops')
 q('S03_iupac',DNA+'NNNRAYYTRTGYTARNNNTGA'+DNA,'S03')
 for name,s in [('internal_stop',DNA[:90]+'TAATAGTGA'+DNA[90:]),('terminal_stop',DNA+'TAA'),('all_n','N'*300),('ambiguous',DNA+'NNNRAY'+DNA)]:q('S04_'+name,s,'S04')
 put('S04_ambiguous.faa',fasta([('ambiguous_protein',AA+'BZJXUO*'+AA)]),'S04','literal ambiguous NCBISTDAA alphabet')
 for name,s in [('lcase',DNA[:72]+DNA[72:144].lower()+DNA[144:]),('all_lcase',DNA.lower()),('seg',DNA+'GCT'*40+DNA),('minus_mask',rc(DNA[:72]+DNA[72:144].lower()+DNA[144:]))]:q('S05_'+name,s,'S05')
 q('S06_seed',DNA,'S06');q('S06_window',DNA[:90]+'GCT'*40+DNA[90:],'S06')
 for n in [1,2,3,9,30]:q('S07_insert'+str(n),DNA[:90]+'GGT'*n+DNA[90:],'S07')
 q('S07_delete',DNA[:90]+DNA[99:],'S07')
 for n in [0,1,2,119,120,121,2999,3000,3001]:q('S08_gap'+str(n),DNA[:120]+'N'*n+DNA[120:],'S08')
 q('S09_biased',DNA+'GCT'*80+DNA,'S09')
 put('S09_subjects.faa',fasta([('balanced',AA),('biased',AA+'A'*80+AA),('short',AA[:30])]),'S09','fixed composition probes')
 for order in ['forward','reverse']:
  ids=list(range(8));ids=ids if order=='forward' else ids[::-1]
  put('S10_'+order+'.faa',fasta([('tie'+str(i),AA) for i in ids]),'S10','8 identical subjects, OID permutation')
 for n in [10001,10002,10003]:
  put('S11_batch'+str(n)+'.fna',fasta([('hit',DNA+'N'*(n-len(DNA))),('nohit','N'*90),('invalid','A'),('second',rc(DNA))]),'S11','first record total nt '+str(n)+'; threshold append semantics')
 # 19410 = 2*(10002-297): first two-chunk quotient boundary, not 10002.
 for n in [19409,19410,19411]:
  s='N'*9600+DNA+'N'*(n-9600-len(DNA))
  q('S12_split'+str(n),s,'S12');q('S12_minus'+str(n),rc(s),'S12')
 put('S13_ids.fna',fasta([('local description',DNA),('ref|NM_000518.5| title',DNA),('duplicate first',DNA),('duplicate second',rc(DNA))]),'S13','literal local/accession/pipe/repeated identities')
 put('S13_ids.faa',fasta([('ref|NP_000509.1| subject title',AA),('local subject description',AA)]),'S13','literal ID styles; parse_deflines omitted')
 q('S14_filter',DNA+DNA,'S14')
 put('S10_251.faa',fasta([(f'tie{i:03d}',AA) for i in range(251)]),'S10','251 identical subjects exposes omitted250 vs explicit500 pairwise alignment count')
 for fam,fn in [('R01','AvCLPV.faa'),('R02','PsCLPV.faa'),('R03','AP027131.faa'),('R04','AP027133.faa'),('R05','PemoMJNVA.faa')]:
  records=list(itertools.islice(SeqIO.parse(ROOT/'LOSAT/tests/fasta'/fn,'fasta'),12))
  put(fam+'_compact.faa',fasta([(r.description,str(r.seq)) for r in records]),fam,'first12 official FAA records; preserve ID/title/sequence,order')
 # Input reduction chosen by pinned NCBI, never LOSAT. Source app/blast/blastx_app.cpp:276-281.
 argv=['/home/kawato/micromamba/bin/blastx','-query',str(ROOT/'LOSAT/tests/fasta/AP027131.fasta'),'-subject',str(a.output/'R03_compact.faa'),'-query_gencode','4','-seg','no','-comp_based_stats','0','-outfmt','6 qstart qend qframe score']
 import os
 env={k:v for k,v in os.environ.items() if k not in ['BATCH_SIZE','CHUNK_SIZE','OVERLAP_CHUNK_SIZE','CTOOLKIT_COMPATIBLE','BL2SEQ_LEGACY','OLD_FSC']};env.update(LC_ALL='C',BLAST_USAGE_REPORT='false')
 result=subprocess.run(argv,env=env,capture_output=True);assert result.returncode==0 and result.stdout
 hits=[list(map(int,line.split())) for line in result.stdout.decode().splitlines()];hit=sorted(hits,key=lambda x:(-x[3],min(x[:2]),max(x[:2]),x[2]))[0]
 lo,hi=sorted(hit[:2]);rec=next(SeqIO.parse(ROOT/'LOSAT/tests/fasta/AP027131.fasta','fasta'));seq=str(rec.seq)[lo-1:hi];seq=seq if hit[2]>0 else rc(seq)
 put('R03_subset.fna',fasta([(f'{rec.id}:{lo}-{hi}:frame{hit[2]}',seq)]),'R03','highest score pinned NCBI compact subject selection;original nt interval and strand preserved')
 # Original bytes report is invariant;argv paths recorded relative to source identities in provenance.
 put('R03_subset_selection.tsv',result.stdout,'R03','pinned BLASTX code4/SEGno/composition0 full-query selection report')
 deriv.append(dict(source='AP027131.fasta',selection_report_sha256=hashlib.sha256(result.stdout).hexdigest(),selection_options=['query_gencode4','segno','composition0','outfmt6 qstart qend qframe score'],original_start_1based=lo,original_end_1based=hi,frame=hit[2],score=hit[3]))
 put('S15_empty.fna',b'','S15','zero byte input');put('S15_invalid.fna',b'not a FASTA nucleotide record\n','S15','literal malformed input')
 put('S02_crlf.fna',fasta([('crlf',DNA)]).replace(b'\n',b'\r\n'),'S02','CRLF equivalent of DNA')
 put('S02_no_final_lf.fna',fasta([('no_final_lf',DNA)]).rstrip(b'\n'),'S02','no final newline')
 # Real bytes remain at pinned Git paths. Only compact derivatives are duplicated.
 for fam,fn in [('R01','AvCLPV.fasta'),('R02','PsCLPV.fasta'),('R03','AP027131.fasta'),('R04','AP027133.fasta'),('R05','MelaMJNV.fasta')]:
  rec=next(SeqIO.parse(ROOT/'LOSAT/tests/fasta'/fn,'fasta'));seq=str(rec.seq)
  for n in [150,300,1000]:
   for sign in [1,-1]:
    subseq=seq[:n] if sign==1 else rc(seq[:n]); name=f'R06_{fam}_{n}_{sign}.fna'
    put(name,fasta([(f'{rec.id}:1-{n}:{sign}',subseq)]),'R06',f'{fn},1-based inclusive1..{n},strand{sign}')
 for fam,fn in [('R01','AvCLPV.gb'),('R02','PsCLPV.gb'),('R07','EDL933.gb'),('R07','Sakai.gb')]:
  rec=next(SeqIO.parse(ROOT/'LOSAT/tests/fasta'/fn,'genbank'));proteins=[]
  for i,f in enumerate(rec.features):
   if f.type!='CDS':continue
   d=dict(index=i,location=str(f.location),qualifiers={k:v for k,v in f.qualifiers.items() if k in ['protein_id','codon_start','transl_table','pseudo','pseudogene','partial','gene']},translation_sha256=hashlib.sha256(f.qualifiers.get('translation',[''])[0].encode()).hexdigest(),selected='translation' in f.qualifiers and 'pseudo' not in f.qualifiers and 'pseudogene' not in f.qualifiers)
   deriv.append(dict(source=fn,**d))
   if d['selected']:proteins.append((f.qualifiers.get('protein_id',[f'CDS{i}'])[0],f.qualifiers['translation'][0]))
  if fam=='R07':put(Path(fn).stem+'.faa',fasta(proteins),fam,'GenBank /translation; no retranslation; omit pseudo/missing translation; preserve CDS order','LATER_LARGE')
  else:
   f=next(f for f in rec.features if f.type=='CDS' and 'translation' in f.qualifiers)
   put(fam+'_cds.fna',fasta([(rec.id+':'+str(f.location),str(f.extract(rec.seq)))]),fam,'first annotated translated CDS extracted in biological strand; annotation manifest')
   put(fam+'_protein.faa',fasta([(f.qualifiers['protein_id'][0],f.qualifiers['translation'][0])]),fam,'same CDS GenBank translation, not LOSAT')
 for acc in ['NM_000518.5','NP_000509.1','NG_000007.3']:
  rec=next(SeqIO.parse(HERE/'inputs/r08'/f'{acc}.gb','genbank'));assert rec.id==acc
  if acc.startswith('NG_'):
   feats=[f for f in rec.features if f.type=='gene' and f.qualifiers.get('gene')==['HBB']];assert len(feats)==1
   f=feats[0];assert len(f.location.parts)==1
   put('R08_HBB_genomic.fna',fasta([(acc+':'+str(f.location),str(f.extract(rec.seq)))]),'R08','annotated gene HBB continuous interval; biological orientation; retains introns')
   deriv.append(dict(source=acc+'.gb',location=str(f.location),start_1based=int(f.location.start)+1,end_1based_inclusive=int(f.location.end),strand=f.location.strand,qualifiers=f.qualifiers))
  else:put('R08_'+acc+('.faa' if acc.startswith('NP_') else '.fna'),fasta([(acc, str(rec.seq))]),'R08','archived accession.version GenBank ORIGIN')
 nroot=a.ncbi_source/'c++/src/algo/blast/unit_tests'
 raw=subprocess.check_output(['/home/kawato/micromamba/bin/blastdbcmd','-db',str(nroot/'seqdb_reader/data/f555'),'-entry','all'])
 put('N01_gi555.fna',raw,'N01','NCBI same-commit f555 database via pinned blastdbcmd -entry all')
 rec=next(SeqIO.parse(HERE/'inputs/ncbi/gi129295.gb','genbank'))
 put('N01_gi129295.faa',fasta([(rec.id,str(rec.seq))]),'N01','archived GI129295 resolved accession.version '+rec.id)
 for gi,length in [(63122693,122347),(112817621,5567),(112585373,5987),(112585216,5531),(112585119,5046)]:
  rec=SeqIO.read(HERE/'inputs/ncbi'/('gi'+str(gi)+'.gb'),'genbank')
  assert len(rec.seq)==length,(gi,rec.id,len(rec.seq),length)
  put('N04_gi'+str(gi)+'.fna',fasta([(rec.id,str(rec.seq))]),'N04','original split_query_unit_test GI'+str(gi)+' -> archived '+rec.id)
  deriv.append(dict(family='N04',gi=gi,accession_version=rec.id,sequence_length=len(rec.seq),archive='inputs/ncbi/gi'+str(gi)+'.gb'))
 for fn in ['bl2seq_unit_test.cpp','blastfilter_unit_test.cpp','linkhsp_unit_test.cpp','split_query_unit_test.cpp']:
  put('ncbi/'+fn,(nroot/'api'/fn).read_bytes(),'N_SOURCE','same-commit original unit source; public domain; preserves internal expected and API profile')
 put('ncbi/split_query.ini',(nroot/'api/data/split_query.ini').read_bytes(),'N04','same-commit original query/context/bounds expected')
 put('derivations.json',json.dumps(deriv,indent=2).encode()+b'\n','PROVENANCE','all CDS qualifiers including partial/codon_start/transl_table and HBB extraction')
 with (a.output/'fixtures.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 print(len(rows),'materialized fixture artifacts')
if __name__=='__main__':main()
