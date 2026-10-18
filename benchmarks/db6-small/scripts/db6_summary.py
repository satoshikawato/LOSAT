#!/usr/bin/env python3
"""dbsearch-20261012 session 14, DB-6: summary of ws_db6_bench.sh's runs.tsv.

usage: db6_summary.py RUNS_TSV
Per case and thread count: median wall time of the three timed runs (1-3) with min-max, LOSAT / NCBI of the
medians, median user CPU and max RSS. Exactness: every LOSAT output of a case (warmup, timed, its own database)
against NCBI -num_threads 1 (all NCBI n1 runs must agree); NCBI n8 / n16 against NCBI n1. Builds: median and
min-max per tool and database.
"""
import csv
import statistics
import sys
from collections import defaultdict


def fnum(x):
    try:
        return float(x)
    except ValueError:
        return float("nan")


rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
cases = []
by = defaultdict(list)
for r in rows:
    if r["case"] not in cases:
        cases.append(r["case"])
    by[r["case"], r["tool"], r["threads"]].append(r)


def timed(c, tool, t):
    return [r for r in by[c, tool, t] if r["run"] in ("1", "2", "3")]


def med(rs, key):
    v = [fnum(r[key]) for r in rs]
    return statistics.median(v) if v else float("nan"), min(v) if v else float("nan"), max(v) if v else float("nan")


bad = 0
print("case\tthreads\tncbi_wall_med[min-max]\tlosat_wall_med[min-max]\tlosat/ncbi\tncbi_user\tlosat_user\tncbi_rss_mib\tlosat_rss_mib\tbusy_max\trepeated")
for c in cases:
    if c.startswith("build-") or c.endswith("-pilot"):
        continue
    for t in ("1", "8", "16"):
        n, l = timed(c, "ncbi", t), timed(c, "losat", t)
        if not n and not l:
            continue
        nw, lw = med(n, "wall_s"), med(l, "wall_s")
        print(f"{c}\t{t}\t{nw[0]:.2f} [{nw[1]:.2f}-{nw[2]:.2f}]\t{lw[0]:.2f} [{lw[1]:.2f}-{lw[2]:.2f}]\t"
              f"{lw[0] / nw[0]:.3f}\t{med(n, 'user_s')[0]:.1f}\t{med(l, 'user_s')[0]:.1f}\t"
              f"{med(n, 'maxrss_kib')[0] / 1024:.0f}\t{med(l, 'maxrss_kib')[0] / 1024:.0f}\t"
              f"{max(int(r[k]) for r in n + l for k in ('busy0', 'busy1'))}\t"
              f"{sum(1 for tool in ('ncbi', 'losat') for r in by[c, tool, t] if 'd' in r['run'])}")
print()
for c in cases:
    if c.startswith("build-") or c.endswith("-pilot"):
        continue
    ref = {r["sha16"] for r in by[c, "ncbi", "1"]}
    rcs = {r["rc"] for k, v in by.items() if k[0] == c for r in v}
    msg = []
    if len(ref) != 1:
        msg.append(f"NCBI n1 runs disagree {sorted(ref)}")
        bad += 1
    los = {(r["tool"], r["threads"], r["sha16"]) for k, v in by.items() if k[0] == c and k[1] != "ncbi" for r in v}
    off = sorted(x for x in los if x[2] not in ref)
    if off:
        msg.append(f"LOSAT != NCBI n1: {off}")
        bad += 1
    else:
        msg.append(f"LOSAT n1/n8/n16 + own db = NCBI n1 ({next(iter(ref))}, {by[c, 'ncbi', '1'][0]['rows']} rows)")
    for t in ("8", "16"):
        s = {r["sha16"] for r in by[c, "ncbi", t]}
        msg.append(f"NCBI n{t} {'= n1' if s <= ref else 'differs from n1 ' + str(sorted(s))}")
    if rcs != {"0"}:
        msg.append(f"exit codes {sorted(rcs)}")
        bad += 1
    print(f"{c}: " + "; ".join(msg))
print()
print("build\ttool\twall_med[min-max]\tmaxrss_mib\tsq_sha16")
for c in cases:
    if not c.startswith("build-"):
        continue
    for tool in ("ncbi5", "ncbi4", "losat"):
        rs = timed(c, tool, "1")
        if rs:
            w = med(rs, "wall_s")
            print(f"{c[6:]}\t{tool}\t{w[0]:.2f} [{w[1]:.2f}-{w[2]:.2f}]\t{med(rs, 'maxrss_kib')[0] / 1024:.0f}\t"
                  f"{sorted({r['sha16'] for r in rs})}")
print(f"\nexactness problems: {bad}")
