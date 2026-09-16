# NCBI reference: c++/src/algo/blast/core/link_hsps.c:827-895
# for (H2_index=H_index-1; H2_index>1;) { ... }
# Read-only emitted-code inspection; no timing or artifact rewriting.
from pathlib import Path
import argparse,subprocess,json,hashlib
p=argparse.ArgumentParser();p.add_argument('version');p.add_argument('--kind',choices=['threaded','serial'],default='threaded');a=p.parse_args()
e=Path(__file__).resolve().parent;d=e/'work'/a.version
artifact=d/'artifacts'/('losat-'+a.kind+'-command.wasm');data=artifact.read_bytes()
if a.kind=='serial':d=d/'serial-emission';d.mkdir(exist_ok=True)
pos=8;names={}
def uint():
    global pos
    v=shift=0
    while True:
        c=data[pos];pos+=1;v|=(c&127)<<shift
        if not c&128:return v
        shift+=7
def string():
    global pos
    n=uint();v=data[pos:pos+n].decode();pos+=n;return v
while pos<len(data):
    sid=data[pos];pos+=1;n=uint();end=pos+n
    if sid==0 and string()=='name':
        while pos<end:
            sub=data[pos];pos+=1;size=uint();nextpos=pos+size
            if sub==1:
                for _ in range(uint()):
                    index=uint();name=string()
                    if any(x in name for x in ['link_hsp_group_ncbi','scan_large_gap_predecessors']):names[index]=name
            pos=nextpos
    pos=end
subprocess.run(['wasm-dis',str(artifact),'-o',str(d/'emitted.wat')],check=True)
lines=(d/'emitted.wat').read_text().splitlines();records=[]
for index,name in names.items():
    start=next(i for i,l in enumerate(lines) if l.startswith(' (func $'+name+' '))
    end=next(i for i in range(start+1,len(lines)) if lines[i].startswith(' (func '))
    body=lines[start:end];(d/f'function-{index}.wat').write_text('\n'.join(body)+'\n')
    records.append(dict(index=index,name=name,wat_lines=len(body),wat_start_line=start+1,locals=sum('(local $' in l for l in body),callee_lines=[l.strip() for l in body if '(call $' in l and 'scan_large_gap_predecessors' in l],wasm_sha256=hashlib.sha256(data).hexdigest()))
(d/'function-indices.json').write_text(json.dumps(records,indent=2)+'\n')
print(json.dumps(records,indent=2))
