#!/usr/bin/env python3
import json, re, subprocess, os, shlex, signal, threading, time
from concurrent.futures import ThreadPoolExecutor
W="/home/kawato/.cache/losat-web-gui-target/e2g-overflow-hunt"
BIN="/home/kawato/.cache/losat-web-gui-target/e2g-bins/pre-T2-overflow/LOSAT"
jobs=[json.loads(l) for l in open(f"{W}/jobs.jsonl")]
lock=threading.Lock()
PAN=re.compile(r"panicked at ([^\n]*?):?\n(.*)")
PAN2=re.compile(r"panicked at (\S+?)(?::\n|\n|$)")
def run(j):
    argv=[BIN,"blastn","-query",j["query"],"-subject",j["subject"],*j["opts"],"-outfmt","6"]
    cmd=" ".join(shlex.quote(a) for a in argv)
    size=os.path.getsize(j["query"])+os.path.getsize(j["subject"])
    p=subprocess.Popen(argv,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,start_new_session=True,env={**os.environ,"RUST_BACKTRACE":"0"})
    try:
        _,err=p.communicate(timeout=300)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid,signal.SIGKILL); p.communicate()
        with lock:
            open(f"{W}/TIMEOUTS.tsv","a").write(f"{size}\t{j['group']}\t{cmd}\n")
        return
    err=err.decode(errors="replace")
    m=re.search(r"panicked at ([^\n]*)\n([^\n]*)",err)
    if m or p.returncode==101:
        if m:
            loc=m.group(1).rstrip(":"); msg=m.group(2)
        else:
            loc="UNKNOWN"; msg=err.strip().replace("\n"," | ")[:300]
        with lock:
            open(f"{W}/PANICS.tsv","a").write(f"{loc}\t{msg}\t{cmd}\t{size}\n")
if not os.path.exists(f"{W}/PANICS.tsv"):
    open(f"{W}/PANICS.tsv","w").write("")
B=100
with ThreadPoolExecutor(16) as ex:
    for b in range(0,len(jobs),B):
        list(ex.map(run,jobs[b:b+B]))
        with lock:
            open(f"{W}/PROGRESS.txt","a").write(f"{time.strftime('%H:%M:%S')} batch {b//B+1}/{(len(jobs)+B-1)//B} done ({min(b+B,len(jobs))}/{len(jobs)} commands)\n")
open(f"{W}/PROGRESS.txt","a").write("ALL DONE\n")
