#!/usr/bin/env python3
"""dbsearch-20261012 session 14: per-CPU use on ws from /proc/stat (Owner 2026-10-17 08:30: decide by per-CPU use,
not by markers alone).

usage:
  cpu_busy.py count SECS NODE MINE   number of CPUs of NODE (0: 0-31,64-95; 1: 32-63,96-127) outside MINE
                                     (a list like 16-31,80-95) busier than 50 % over SECS seconds
  cpu_busy.py pick SECS N STATE_DIR  N physical cores of one node (node 0 first: 24-31, 16-23, 0-15; then node 1:
                                     48-63, 32-47; node 1's untimed gates float over the whole node) whose cores and HT siblings were under 10 % busy over SECS seconds
                                     and are not named by another marker, SHARED ones included (a "WHOLE node" one
                                     takes its node); prints "NODE CORES SIBLINGS" (comma lists) or nothing
  cpu_busy.py pick SECS N STATE_DIR polite MY_QUEUE_MARKER   the same, but (1) skip a node whose TIMED marker asks to
                                     keep heavy work off it, and (2) while an older NODE_QUEUED_ marker waits (those
                                     queues take the first node with no NODE<n>_ marker), take only a node that already
                                     has a marker, so the older queue still gets the next free node
"""
import glob
import os
import re
import sys
import time

NODES = {0: list(range(0, 32)) + list(range(64, 96)), 1: list(range(32, 64)) + list(range(96, 128))}


def cpus(spec):
    out = set()
    for part in filter(None, spec.split(",")):
        a, _, b = part.partition("-")
        out.update(range(int(a), int(b or a) + 1))
    return out


def stat():
    s = {}
    for line in open("/proc/stat"):
        m = re.match(r"cpu(\d+) (.*)", line)
        if m:
            f = [int(x) for x in m.group(2).split()[:8]]
            s[int(m.group(1))] = (sum(f), f[3] + f[4])
    return s


def busy(secs):
    a = stat()
    time.sleep(secs)
    b = stat()
    return {c: 1 - (b[c][1] - a[c][1]) / max(1, b[c][0] - a[c][0]) for c in b}


def booked(state):
    taken = set()
    keep_off.clear()
    marked.clear()
    for p in glob.glob(os.path.join(state, "NODE[01]_*")):
        if "dbsearch-s14-db6" in os.path.basename(p):
            continue
        marked.add(int(os.path.basename(p)[4]))
        text = open(p, errors="replace").read()
        node = int(re.match(r"NODE(\d)", os.path.basename(p)).group(1))
        if re.search(r"keep [^.;]*off", text, re.I) and re.search(r"\bTIMED\b", text):
            keep_off.add(node)
        if re.search(r"whole", text, re.I):
            taken.update(NODES[node])
        for m in re.finditer(r"CPUs? ([0-9][0-9,\-]*)", text):
            for c in cpus(m.group(1).rstrip(",-")):
                taken.update({c, c + 64 if c < 64 else c - 64})
    return taken


keep_off = set()
marked = set()
mode = sys.argv[1]
if mode == "overlap":  # overlap STATE_DIR CPUS MY_MARKER: other markers no newer than mine that name one of CPUS
    state, mine, me = sys.argv[2], cpus(sys.argv[3]), sys.argv[4]
    for p in glob.glob(os.path.join(state, "NODE[01]_*")):
        if os.path.basename(p) == os.path.basename(me) or os.path.getmtime(p) > os.path.getmtime(me):
            continue
        text = open(p, errors="replace").read()
        named = set(NODES[int(os.path.basename(p)[4])]) if re.search(r"whole", text, re.I) else set()
        for m in re.finditer(r"CPUs? ([0-9][0-9,\-]*)", text):
            named |= cpus(m.group(1).rstrip(",-"))
        if named & mine:
            print(os.path.basename(p))
if mode == "count":
    secs, node, mine = float(sys.argv[2]), int(sys.argv[3]), cpus(sys.argv[4])
    u = busy(secs)
    print(sum(1 for c in NODES[node] if c not in mine and u[c] > 0.5))
elif mode == "pick":
    secs, n, state = float(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    taken = booked(state)
    u = busy(secs)
    order = {1: list(range(48, 64)) + list(range(32, 48)), 0: list(range(24, 32)) + list(range(16, 24)) + list(range(0, 16))}
    polite = len(sys.argv) > 6 and sys.argv[5] == "polite"
    older_queue = False
    if polite:
        mine = os.path.getmtime(sys.argv[6])
        older_queue = any(os.path.getmtime(q) < mine for q in glob.glob(os.path.join(state, "NODE_QUEUED_*"))
                          if os.path.basename(q) != os.path.basename(sys.argv[6]))
    for node in (0, 1):
        if polite and (node in keep_off or (older_queue and node not in marked)):
            continue
        free = [c for c in order[node] if c not in taken and c + 64 not in taken and u[c] < 0.10 and u[c + 64] < 0.10]
        if len(free) >= n:
            block = free[:n]
            print(node, ",".join(map(str, block)), ",".join(str(c + 64) for c in block))
            break
