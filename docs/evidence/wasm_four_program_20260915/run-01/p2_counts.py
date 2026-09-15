# NCBI reference: c++/src/algo/blast/core/blast_kappa.c:3329-3334,3493-3504
# NRrecord_tld[i] = Blast_CompositionWorkspaceNew(); query_info = query_info_tld[tid];
# Diagnostic-only work-unit initialization and retained scratch measurements.
from pathlib import Path
import shutil,subprocess,json
b=Path(__file__).resolve().parent;dest=b/'P2counts-source';shutil.copytree(b/'P2-source',dest,dirs_exist_ok=True)
p=dest/'LOSAT/src/algorithm/blastp/blast_engine.rs';subprocess.run(['rustfmt','--edition','2021',str(p)],check=True);s=p.read_text();s=s.replace('let query_workspace = query_workspace.get_or_init(|| {','let query_workspace = query_workspace.get_or_init(|| { eprintln!("P2_QUERY_INIT");',1)
needle='.map_init(\n                        || {';assert s.count(needle)==1;s=s.replace(needle,needle+'\n                            eprintln!("P2_SLOT_INIT");')
start=s.index('let query_workspace = std::sync::OnceLock::new();');i=s.index('.map(Some)',start);s=s[:i]+s[i:].replace('.map(Some)','.map(|r| {eprintln!("P2_RETAIN {} {}",kappa_gap_scratch.diagnostic_bytes(),kappa_preliminary_hits.capacity());Some(r)})',1);p.write_text(s)
p=dest/'LOSAT/src/algorithm/blastp/gapalign.rs';s=p.read_text();needle='impl GapAlignScratch {';s=s.replace(needle,needle+'''
    pub(crate) fn diagnostic_bytes(&self)->usize {
        self.dp_mem.capacity()*std::mem::size_of::<BlastGapDp>()
        +self.trace_rows.capacity()*std::mem::size_of::<Vec<u8>>()
        +self.trace_rows.iter().map(Vec::capacity).sum::<usize>()
        +self.trace_offsets.capacity()*std::mem::size_of::<usize>()
        +self.trace_ops_reversed.capacity()*std::mem::size_of::<(u8,u32)>()
    }
''',1);p.write_text(s)
subprocess.run(['python3',str(b/'build.py'),'P2counts','--kinds','native'],check=True)
rows=[]
for n in [1,2,4,8]:
 cmd=[str(b/'P2counts/native-command/release/LOSAT'),'blastp','-query',str(b/'single-query/query.faa'),'-subject',str(b/'single-query/subjects.faa'),'-outfmt','6','-num_threads',str(n),'-out',str(b/f'P2counts-n{n}.out')]
 r=subprocess.run(cmd,capture_output=True,check=True);(b/f'P2counts-n{n}.stderr').write_bytes(r.stderr)
 text=r.stderr.decode();retain=[list(map(int,l.split()[1:])) for l in text.splitlines() if l.startswith('P2_RETAIN ')]
 rows.append(dict(n=n,argv=cmd,query_init=text.count('P2_QUERY_INIT'),slot_init=text.count('P2_SLOT_INIT'),retained_observations=retain,raw_equal=(b/f'P2counts-n{n}.out').read_bytes()==(b/'single-query/oracle.out').read_bytes()))
(b/'P2-counts.json').write_text(json.dumps(rows,indent=2)+'\n');print(rows)
