import random,sys
seed=int(sys.argv[1]); nseq=int(sys.argv[2]); fn=sys.argv[3]
r=random.Random(seed)
sizes=[0,1,2,3,7,100,1000,4095,4096,4097,8180,8185,8186,8187,8188,8189,8190,8191,8192,8193,8194,8195,12000,16380,16381,16382,16383,16384,16385,20000]
ops=[("W",30),("I",5),("P",15),("F",10),("E",6),("S",4),("C",3),("X",5),("T",4),("O",2),("OB",1),("OF",1),("OM",0.5),("K",2)]
tot=sum(w for _,w in ops)
out=[]
for _ in range(nseq):
    lim=-1 if r.random()<0.35 else r.choice([0,1,5,50,100,1000,4096,8190,8191,8192,8193,9000,16381,16382,16383,20000,30000,r.randrange(0,40000)])
    mode=r.randrange(2)
    out.append(f"R {lim} {mode}")
    out.append(r.choice(["O"]*20+["OB","OF","OM"]))
    for _ in range(r.randrange(1,40)):
        x=r.random()*tot
        for o,w in ops:
            x-=w
            if x<=0: break
        if o in("W","I"): out.append(f"{o} {r.choice(sizes) if r.random()<0.85 else r.randrange(0,25000)}")
        elif o=="S": out.append(f"S {r.choice([0,5,8,100,4096,8192,65536])}")
        elif o=="T": out.append(f"T {r.choice([1,2,4,5,0])}")
        elif o=="K": out.append(f"K {r.choice([5,28,122,0,4])}")
        else: out.append(o)
    out.append("Z")
open(fn,"w").write("\n".join(out)+"\n")
