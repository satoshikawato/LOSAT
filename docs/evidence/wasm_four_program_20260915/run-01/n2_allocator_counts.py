# NCBI reference: c++/src/algo/blast/core/blast_extend.c:219-224,337-341
# malloc(MIN_INIT_HITLIST_SIZE * sizeof(BlastInitHSP)); realloc(...);
# Count actual native System allocator calls, outside performance samples.
import hashlib,json,re,shutil,subprocess
from pathlib import Path
b=Path(__file__).resolve().parent;t=Path('/mnt/c/Users/genom/GitHub/LOSAT/LOSAT/tests/fasta');records=[]
for v in ['baseline','N2']:
 name=v+'alloc';dest=b/(name+'-source');shutil.copytree(b/(v+'-source'),dest)
 shutil.copyfile(b/'count_allocations.rs',dest/'LOSAT/src/count_allocations.rs')
 p=dest/'LOSAT/src/main.rs';s=p.read_text();s=s.replace('use anyhow::Result;', '// NCBI reference: c++/src/algo/blast/core/blast_extend.c:219-224,337-341\n// malloc(MIN_INIT_HITLIST_SIZE * sizeof(BlastInitHSP)); realloc(...);\nmod count_allocations;\nuse anyhow::Result;')
 assert 'mod count_allocations;' in s
 assert s.count('    Ok(())')==1;s=s.replace('    Ok(())','    count_allocations::report();\n    Ok(())');p.write_text(s)
 subprocess.run(['python3',str(b/'build.py'),name,'--kinds','native'],check=True)
 for case,q,s in [('short','LC738874.fasta','LC738870.fasta'),('large','AP027202.fasta','LC738875.fasta')]:
  output=b/f'{name}-{case}.out';cmd=[str(b/name/'native-command/release/LOSAT'),'blastn','-task','blastn','-query',str(t/q),'-subject',str(t/s),'-num_threads','8','-outfmt','6','-out',str(output)]
  ret=subprocess.run(cmd,capture_output=True,check=True);(b/f'{name}-{case}.stderr').write_bytes(ret.stderr)
  line=next(x for x in ret.stderr.decode().splitlines() if x.startswith('[ALLOCATOR]'))
  records.append({'version':v,'case':case,'argv':cmd,'sha256':hashlib.sha256(output.read_bytes()).hexdigest(),'counts':{key:int(value) for key,value in re.findall(r'(\w+)=(\d+)',line)},'purpose':'untimed allocator diagnostic; requested bytes are not RSS or actual allocator usable bytes'})
  (b/'N2-allocator-counts.json').write_text(json.dumps(records,indent=2)+'\n')
