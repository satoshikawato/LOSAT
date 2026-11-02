#!/usr/bin/env python3
"""Line WO: table of an ab run. usage: ab_table.py runs.tsv [first round used, default 1]"""
import csv, statistics as st, sys
from collections import defaultdict
first = int(sys.argv[2]) if len(sys.argv) > 2 else 1
rows = [r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")]
by = defaultdict(dict)
bad = [r for r in rows if r["rc"] != "0"]
for r in rows:
    if r["rc"] == "0" and int(r["round"]) >= first:
        by[(r["case"], int(r["round"]))][r["arm"]] = r
cases = []
for r in rows:
    if r["case"] not in cases: cases.append(r["case"])
print("case | pairs | wall dev s | fix/dev wall median (min-max) | null/dev wall median (min-max) | RSS fix/dev | RSS null/dev | sha same")
gf, gn = [], []
for c in cases:
    pairs = [v for (cc, _), v in by.items() if cc == c and {"dev", "null", "fix"} <= set(v)]
    if not pairs: continue
    ratio = lambda arm, key: [float(p[arm][key]) / float(p["dev"][key]) for p in pairs]
    wf, wn, mf, mn = ratio("fix", "wall"), ratio("null", "wall"), ratio("fix", "maxrss_kb"), ratio("null", "maxrss_kb")
    shas = {p[a]["sha16"] for p in pairs for a in ("dev", "null", "fix")}
    gf.append(st.median(wf)); gn.append(st.median(wn))
    print(f"{c} | {len(pairs)} | {st.median(float(p['dev']['wall']) for p in pairs):.2f} | x{st.median(wf):.3f} ({min(wf):.3f}-{max(wf):.3f}) | x{st.median(wn):.3f} ({min(wn):.3f}-{max(wn):.3f}) | x{st.median(mf):.3f} | x{st.median(mn):.3f} | {'yes' if len(shas) == 1 else 'NO ' + str(shas)}")
if gf:
    print(f"geometric mean of medians: fix x{st.geometric_mean(gf):.3f}, null x{st.geometric_mean(gn):.3f}")
print(f"failed runs: {len(bad)}", [(r['case'], r['arm'], r['round'], r['rc']) for r in bad][:8])
