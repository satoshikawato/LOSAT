#!/usr/bin/env python3
# parse.py NAME  -> one line per operation: "seq# op-text | events ; state" read from NAME.log and NAME.state
import sys, re
def fix(e, name):
    pre='out.'+name+' '
    e=e[len(pre):] if e.startswith(pre) else e
    return 'OPEN' if e=='open' else e

name = sys.argv[1]
log = open(f"{name}.log").read().split("\n")
states = [l.split()[1] for l in open(f"{name}.state").read().split("\n") if l.startswith("STATE")]
files = [l for l in open(f"{name}.state").read().split("\n") if l.startswith("FILE")]
out = []
si = 0; fi = 0
cur = None
ev = []
seq = None
def flush():
    global cur, ev
    if cur is not None:
        out.append((cur, ev))
    cur = None; ev = []
for l in log:
    if l.startswith("SEQ"):
        flush(); seq = l; cur = ("SEQ", l)
    elif l.startswith("OP "):
        flush(); cur = ("OP", l)
    elif l == "DESTROY":
        flush(); cur = ("DESTROY", l)
    elif l == "END":
        flush()
    elif l:
        ev.append(l)
flush()
# attach states: each OP has a state, in order
res = []
for kind, l in [(c[0], c[1]) for c, _ in out]:
    pass
si = 0; fi = 0
for (c, ev) in out:
    kind, l = c
    if kind == "OP":
        st = states[si]; si += 1
        res.append(f"{l[3:]} | {' ; '.join(fix(e,name) for e in ev)} | st={st}")
    elif kind == "SEQ":
        res.append(l)
    else:
        res.append(f"DESTROY | {' ; '.join(fix(e,name) for e in ev)} | {files[fi]}"); fi += 1
print("\n".join(res))
