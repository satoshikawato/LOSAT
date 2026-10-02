import subprocess,os
NCBI="/home/kawato/micromamba/bin/blastn"; L="/home/kawato/.cache/losat-web-gui-target/e2g-bins/AC1/LOSAT"
def run(cmd,env):
    e={k:v for k,v in os.environ.items() if k not in ("CHUNK_SIZE","OVERLAP_CHUNK_SIZE")}; e.update(env)
    return subprocess.run(cmd,capture_output=True,env=e).stdout.decode()
base=["-query","mjnv60k_q.fa","-subject","mjnv60k_s.fa","-task","blastn","-outfmt","6"]
U=set(run([NCBI,*base],{"CHUNK_SIZE":"100000000"}).splitlines())
for C,O,Ob in [(8000,4100,4090),(20000,11250,10587)]:
    LS=set(run([L,"blastn",*base],{"CHUNK_SIZE":str(C),"OVERLAP_CHUNK_SIZE":str(O)}).splitlines())
    NB=set(run([NCBI,*base],{"CHUNK_SIZE":str(C),"OVERLAP_CHUNK_SIZE":str(Ob)}).splitlines())
    LB=set(run([L,"blastn",*base],{"CHUNK_SIZE":str(C),"OVERLAP_CHUNK_SIZE":str(Ob)}).splitlines())
    print("==",C,O,"unsplit-only:"); [print("  ",x) for x in sorted(U-LS)]
    print("  losat-only:"); [print("  ",x) for x in sorted(LS-U)]
    print("  ncbi nonresplit O=%d unsplit-only:"%Ob); [print("  ",x) for x in sorted(U-NB)]
    print("  losat==ncbi at nonresplit O:", LB==NB)
