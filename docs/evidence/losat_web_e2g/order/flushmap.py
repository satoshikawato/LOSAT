import subprocess, sys, re, shlex
G="/tmp/claude-1000/-mnt-c-Users-genom-GitHub-LOSAT/19b36304-5bfc-441a-88bc-4d1f633795a3/scratchpad"
args=sys.argv[1:]
out=G+"/fm.txt"
script=open(G+"/wtrace.gdb").read().replace("\nrun\n", "\nrun " + " ".join(shlex.quote(a) for a in args) + f" > {out}\n")
open(G+"/wtrace_run.gdb","w").write(script)
p=subprocess.run(["gdb","-q","-batch","-x",G+"/wtrace_run.gdb","/home/kawato/micromamba/bin/blastn"],capture_output=True,text=True)
ev=[]
lines=[l for l in p.stdout.splitlines() if l.startswith("W fd=")]
for l in lines[::2]:
    fd,n=map(int,re.findall(r"-?\d+",l))
    ev.append((fd,n))
data=open(out,"rb").read()
pos=0
for fd,n in ev:
    if fd==2: print(f"   <stderr {n} bytes>"); continue
    chunk=data[pos:pos+n]; pos+=n
    print(f"[{n:5d}] start={chunk[:40]!r} ... end={chunk[-50:]!r}")
print("total",pos,"file",len(data))
