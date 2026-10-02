import json,subprocess,os
NCBI="/home/kawato/micromamba/bin/blastn"; L="/home/kawato/.cache/losat-web-gui-target/e2g-bins/AC1/LOSAT"
rows=json.load(open("results.json"))
def run(cmd,env):
    e={k:v for k,v in os.environ.items() if k not in ("CHUNK_SIZE","OVERLAP_CHUNK_SIZE","BATCH_SIZE")}; e.update(env)
    return subprocess.run(cmd,capture_output=True,env=e).stdout
out=[]
for r in rows:
    for fmt in ("6","0"):
        base=["-query",f"{r['tag']}_q.fa","-subject",f"{r['tag']}_s.fa","-task",r["task"],"-outfmt",fmt]
        l=run([L,"blastn",*base],{"CHUNK_SIZE":str(r["C"]),"OVERLAP_CHUNK_SIZE":str(r["O"])})
        n=run([NCBI,*base],{"CHUNK_SIZE":str(r["C"]),"OVERLAP_CHUNK_SIZE":str(r["ncbi_nonresplit_O"])})
        out.append((r["tag"],r["task"],r["C"],r["O"],r["ncbi_nonresplit_O"],fmt,l==n))
        print(*out[-1],flush=True)
json.dump(out,open("exact.json","w"))
print("equal",sum(x[-1] for x in out),"of",len(out))
