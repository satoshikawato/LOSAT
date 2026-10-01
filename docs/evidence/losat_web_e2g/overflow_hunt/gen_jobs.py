#!/usr/bin/env python3
"""Deterministic job generator for the BLASTN overflow hunt. Writes jobs.jsonl."""
import json, random, os
from pathlib import Path
R = Path("/mnt/c/Users/genom/GitHub/LOSAT-web-gui/LOSAT")
W = Path("/home/kawato/.cache/losat-web-gui-target/e2g-overflow-hunt")
F0 = R/"tests/fasta/outfmt0"; FX = R/"tests/fixtures/blastn_regression/inputs"; FA = R/"tests/fasta"
jobs = []
def add(group, q, s, opts):
    jobs.append({"group": group, "query": str(q), "subject": str(s), "opts": opts})
# 1
PAIRS = [("strand_query","strand_subject"),("mask_query","mask_subject"),("longdef_query","longdef_subject"),
         ("many_query","many_subject"),("multi_query","multi_subject"),("width_query","width_subject"),
         ("lcase_minus_query","lcase_minus_subject"),("iupac_query","multi_subject")]
for q,s in PAIRS:
    q=F0/f"{q}.fasta"; s=F0/f"{s}.fasta"
    for ws in (4,5,6,7,8,11):
        for ev in ("10","1e5"):
            add("sweep", q, s, ["-task","blastn","-word_size",str(ws),"-evalue",ev])
        add("sweep_sc11", q, s, ["-task","blastn","-word_size",str(ws),"-reward","1","-penalty","-1","-gapopen","3","-gapextend","2"])
        add("sweep_sc34", q, s, ["-task","blastn","-word_size",str(ws),"-reward","3","-penalty","-4","-gapopen","10","-gapextend","3"])
    for ws in (11,16,28):
        for ev in ("10","1e5"):
            add("sweep_mb", q, s, ["-task","megablast","-word_size",str(ws),"-evalue",ev])
# 3
FP = [("q3k","s600"),("q3k","s600ties"),("mq","s60k"),("ambiguity_query_c","ambiguity_subject"),
      ("lcase_island_blastn_query","lcase_island_blastn_subject"),("lcase_island_megablast_query","lcase_island_megablast_subject"),
      ("mq","sakai_60k"),("s60k","sakai_60k"),("q3k","sakai_60k")]
for q,s in FP:
    for ws in (4,7):
        add("fixture", FX/f"{q}.fa", FX/f"{s}.fa", ["-task","blastn","-word_size",str(ws),"-evalue","1e5"])
        if q.startswith("lcase"):
            add("fixture_lc", FX/f"{q}.fa", FX/f"{s}.fa", ["-task","blastn","-word_size",str(ws),"-evalue","1e5","-lcase_masking"])
# 2
def load(p):
    seq=[]
    with open(p) as f:
        for l in f:
            if l.startswith(">"):
                if seq: break
                continue
            seq.append(l.strip())
    return "".join(seq).upper()
genomes = {n: load(FA/n) for n in ["EDL933.fna","Sakai.fna","MG1655.fna","LC738868.fasta","AvCLPV.fasta","AP027078.fasta","NZ_CP006932.fasta","MeenMJNV.fasta","LvMJNV.fasta"]}
gnames = sorted(genomes)
COMP = str.maketrans("ACGTNRYKMSWBDHV","TGCANYRMKSWVHDB")
def rc(s): return s.translate(COMP)[::-1]
def mutate(rng, s, sub, indel):
    out=[]; 
    for c in s:
        r=rng.random()
        if r<sub: out.append(rng.choice("ACGT"))
        elif r<sub+indel/2: continue
        elif r<sub+indel:
            out.append(c); out.append("".join(rng.choice("ACGT") for _ in range(rng.choice((1,1,1,2,3,5)))))
        else: out.append(c)
    return "".join(out)
def slice_(rng, g, n):
    st=rng.randrange(0,len(g)-n); return st, g[st:st+n]
def fasta(recs):
    o=[]
    for name,s in recs:
        o.append(">"+name)
        for i in range(0,len(s),70): o.append(s[i:i+70])
    return "\n".join(o)+"\n"
rng = random.Random(20261002)
(W/"rand").mkdir(exist_ok=True)
def rand_case(i, kind):
    gs = genomes[rng.choice(gnames)]
    n_s = int(10**rng.uniform(2.3,4.3)); n_s=max(200,min(20000,n_s))
    st, sub = slice_(rng, gs, n_s)
    if kind=="plain":
        n_q = max(200,min(20000,int(10**rng.uniform(2.3,4.3))))
        if rng.random()<0.75:
            ln=min(n_q,n_s); off=rng.randrange(0,n_s-ln+1) if n_s>ln else 0
            q = sub[off:off+ln]
            q = mutate(rng,q,rng.uniform(.01,.10),rng.uniform(0,.02))
            if rng.random()<.5: q=rc(q)
            if rng.random()<.3:
                _,fl=slice_(rng,genomes[rng.choice(gnames)],rng.randrange(50,1000)); q=fl+q
        else:
            g2=genomes[rng.choice(gnames)]; _,q=slice_(rng,g2,n_q)
            q=mutate(rng,q,rng.uniform(.01,.10),rng.uniform(0,.02))
        qrecs=[("q%d"%i,q)]
    elif kind=="rcpal":
        n_x=max(200,min(10000,int(10**rng.uniform(2.3,4.0))))
        _,x=slice_(rng,gs,n_x)
        x=mutate(rng,x,rng.uniform(.0,.10),rng.uniform(0,.02))
        y=rc(x)
        if rng.random()<.5: y=mutate(rng,y,rng.uniform(.01,.10),rng.uniform(0,.02))
        qrecs=[("q%d"%i,x+y)]
        if rng.random()<.5: sub=mutate(rng,x[:n_s],rng.uniform(.01,.1),rng.uniform(0,.02))
    else:
        nq=rng.randint(2,20); qrecs=[]
        for k in range(nq):
            ln=rng.randint(50,400)
            if rng.random()<.7 and len(sub)>=ln:
                o=rng.randrange(0,len(sub)-ln+1); q=sub[o:o+ln]
            else: _,q=slice_(rng,gs,ln)
            q=mutate(rng,q,rng.uniform(.01,.10),rng.uniform(0,.02))
            if rng.random()<.5: q=rc(q)
            qrecs.append(("q%d_%d"%(i,k),q))
    if rng.random()<.1:
        qrecs=[(n,"".join(rng.choice("NRY") if rng.random()<.003 else c for c in s)) for n,s in qrecs]
    if rng.random()<.1:
        sub=sub.lower() if rng.random()<.5 else sub[:len(sub)//3]+sub[len(sub)//3:2*len(sub)//3].lower()+sub[2*len(sub)//3:]
    qp=W/"rand"/f"{kind}{i:03d}_q.fa"; sp=W/"rand"/f"{kind}{i:03d}_s.fa"
    qp.write_text(fasta(qrecs)); sp.write_text(fasta([("s%d"%i,sub)]))
    for ws in (4,5,6,7):
        add("rand_"+kind, qp, sp, ["-task","blastn","-word_size",str(ws),"-evalue","1e5"])
    add("rand_"+kind, qp, sp, [])
for i in range(400): rand_case(i,"plain")
for i in range(50): rand_case(i,"rcpal")
for i in range(50): rand_case(i,"multi")
with open(W/"jobs.jsonl","w") as f:
    for j in jobs: f.write(json.dumps(j)+"\n")
print(len(jobs))
