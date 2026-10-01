import collections, os, re
W="/home/kawato/.cache/losat-web-gui-target/e2g-overflow-hunt"
rows=[l.rstrip("\n").split("\t") for l in open(f"{W}/PANICS.tsv") if l.strip()]
g=collections.defaultdict(list)
for loc,msg,cmd,size in rows: g[(loc,msg)].append((int(size),cmd))
nj=sum(1 for _ in open(f"{W}/jobs.jsonl"))
to=[l for l in open(f"{W}/TIMEOUTS.tsv")] if os.path.exists(f"{W}/TIMEOUTS.tsv") else []
out=[f"# Overflow hunt summary\n\nCommands run: {nj}; panics: {len(rows)}; timeouts: {len(to)}\nBinary: pre-T2-overflow LOSAT. Inputs of random cases kept in {W}/rand/\n"]
for (loc,msg),v in sorted(g.items(),key=lambda x:-len(x[1])):
    v.sort()
    out.append(f"\n## {loc}\n\nMessage: {msg}\n\nCount: {len(v)}\n\n3 shortest (by total query+subject bytes):\n")
    seen=set()
    for s,c in v:
        if c in seen: continue
        seen.add(c)
        out.append(f"- [{s} bytes] `{c.replace('/home/kawato/.cache/losat-web-gui-target/e2g-bins/pre-T2-overflow/','')}`")
        if len(seen)==3: break
out.append("\n## Timeouts\n"+("".join(f"- {l}" for l in to) if to else "none\n"))
open(f"{W}/SUMMARY.md","w").write("\n".join(out)+"\n")
print("\n".join(out)[:4000])
