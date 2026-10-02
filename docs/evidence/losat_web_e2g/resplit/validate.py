import random, subprocess, os, sys, json
from pathlib import Path
W=Path("."); NCBI="/home/kawato/micromamba/bin/blastn"; L="/home/kawato/.cache/losat-web-gui-target/e2g-bins/AC1/LOSAT"
FA="/mnt/c/Users/genom/GitHub/LOSAT-web-gui/LOSAT/tests/fasta"
def genome(n): return "".join(l.strip() for l in open(f"{FA}/{n}").read().splitlines()[1:])
rng=random.Random(20261002)
def mutate(s,r):
    return "".join(c if rng.random()>r else rng.choice([b for b in "ACGT" if b!=c]) for c in s)
inputs={}
for tag,gn,off,qlen in [("edl30k","EDL933.fna",2000000,30000),("mjnv60k","MejoMJNV.fasta",50000,60000),("sakai45k","Sakai.fna",3500000,45000)]:
    g=genome(gn); q=g[off:off+qlen]
    parts=[]
    for k in range(14):
        a=rng.randrange(0,qlen-800); ln=rng.randint(200,800)
        parts.append("".join(rng.choice("ACGT") for _ in range(rng.randint(50,300))))
        w=mutate(q[a:a+ln], rng.choice((0.0,0.02,0.05)))
        parts.append(w)
    (W/f"{tag}_q.fa").write_text(f">{tag}\n{q}\n"); (W/f"{tag}_s.fa").write_text(f">{tag}_s\n{''.join(parts)}\n")
    inputs[tag]=qlen
def calc(C,O,Ln):
    n=Ln//(C-O) if C>O else 0
    if n<=1: return 1,Ln
    c=(Ln+(n-1)*O)//n
    if n<c-O: c+=1
    return n,c
def ranges(C,O,Ln):
    n,c=calc(C,O,Ln); out=[]; s=0
    for k in range(n):
        e=s+c
        if e>=Ln or (e<Ln and k+1==n): e=Ln
        out.append((s,e)); s+=c-O
        if s>Ln or e==Ln: break
    return out
def run(cmd,env):
    p=subprocess.run(cmd,capture_output=True,env={**{k:v for k,v in os.environ.items() if k not in ("BATCH_SIZE","CHUNK_SIZE","OVERLAP_CHUNK_SIZE")},**env})
    return p.returncode,p.stdout.decode(),p.stderr.decode()
def resplit(C,O,Ln):
    return any(calc(C,O,e-s)[0]>1 for s,e in ranges(C,O,Ln))
rows=[]
for tag,Ln in inputs.items():
    q=str(W/f"{tag}_q.fa"); s=str(W/f"{tag}_s.fa")
    for task in ("blastn","megablast"):
        base=["-query",q,"-subject",s,"-task",task,"-outfmt","6"]
        n_unsplit=run([NCBI,*base],{"CHUNK_SIZE":"100000000"})
        for C,O in [(3000,1526),(3000,2200),(3000,2900),(8000,4100),(8000,7000),(20000,11250),(20000,15000)]:
            assert resplit(C,O,Ln), (C,O,Ln)
            l=run([L,"blastn",*base],{"CHUNK_SIZE":str(C),"OVERLAP_CHUNK_SIZE":str(O)})
            n=run([NCBI,*base],{"CHUNK_SIZE":str(C),"OVERLAP_CHUNK_SIZE":str(O)})
            # NCBI with a split that does not re-split (largest overlap below the boundary), for comparison
            Ob=O
            while Ob>0 and resplit(C,Ob,Ln): Ob-=1
            nb=run([NCBI,*base],{"CHUNK_SIZE":str(C),"OVERLAP_CHUNK_SIZE":str(Ob)})
            zones=[]
            rg=ranges(C,O,Ln)
            for k in range(len(rg)-1): zones.append((rg[k+1][0],rg[k][1]))
            def hsps(t): return set(t.splitlines())
            U=hsps(n_unsplit[1]); LS=hsps(l[1]); NB=hsps(nb[1])
            def touches(line):
                f=line.split("\t"); a,b=sorted((int(f[6])-1,int(f[7])))
                return any(a<z1 and b>z0-0 for z0,z1 in zones) or any(a<e<b for _,e in rg) or any(a<st<b for st,_ in rg)
            lonly=LS-U; uonly=U-LS
            bad=[x for x in lonly|uonly if not touches(x)]
            rows.append(dict(tag=tag,task=task,C=C,O=O,ncbi_exit=n[0],losat_exit=l[0],losat_hsps=len(LS),unsplit_hsps=len(U),losat_only=len(lonly),unsplit_only=len(uonly),diff_not_at_chunk_zone=len(bad),ncbi_nonresplit_O=Ob,ncbi_nonresplit_vs_unsplit=len(NB^U),losat_stderr=l[2][:80]))
            print(json.dumps(rows[-1]),flush=True)
json.dump(rows,open("results.json","w"),indent=1)
