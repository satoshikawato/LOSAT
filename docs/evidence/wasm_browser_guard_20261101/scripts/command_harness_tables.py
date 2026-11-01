#!/usr/bin/env python3
# Tables of the gbdraw-type harness from its run logs: passes (with the load when the pass started),
# one output hash per job of the dev arms, traps with the bytes at the reported offset.
import collections, csv, json, re, sys
from pathlib import Path

H = Path(__file__).resolve().parent
L = H / "logs"
OUT = Path(sys.argv[1])
SOURCES = [("a", L / "b/run-gshipped-aborted.jsonl", L / "rounds-a-aborted.log"),
           ("a2", L / "b/run-gshipped-cliv1.jsonl", L / "rounds-a-cliv1.log"),
           ("b", L / "b/run-gshipped.jsonl", L / "rounds-b.log"),
           ("b", L / "b/run-devraw.jsonl", L / "rounds-b.log"),
           ("b", L / "b/run-devguarded.jsonl", L / "rounds-b.log"),
           ("c", L / "run-devraw.jsonl", L / "rounds-c.log"),
           ("c", L / "run-devguarded.jsonl", L / "rounds-c.log")]
def loads(log, arm):
    return [float(m.group(1)) for m in re.finditer(rf"^== round \d+ arm {arm} \S+ load (\S+)", log.read_text(), re.M)]
passes, jobs = [], collections.defaultdict(lambda: collections.defaultdict(set))
for run, path, log in SOURCES:
    rows = [json.loads(line) for line in path.open() if line.strip()]
    rows = [r for r in rows if r["out"]]
    ls = loads(log, rows[0]["arm"])
    for i, r in enumerate(rows):
        passes.append([run, r["arm"], r["chromium"], len(r["out"]), r["traps"], r["s"], round(max(r["sizes"]) / 1048576), ls[i] if i < len(ls) else ""])
        if r["arm"] != "gshipped":
            for o in r["out"]:
                jobs[(f"jobs{len(r['out'])}", o["index"])][(o["sha"], o["len"])].add(r["arm"])
with (OUT / "command_harness_passes.tsv").open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow("run arm chromium jobs traps seconds max_main_memory_MiB load_at_start".split()); w.writerows(passes)
with (OUT / "command_harness_outputs.tsv").open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow("job_list job output_sha256 output_bytes arms".split())
    for (name, index), seen in sorted(jobs.items()):
        for (sha, n), arms in sorted(seen.items()):
            w.writerow([name, index, sha, n, ",".join(sorted(arms))])
print("passes", len(passes), "jobs", len(jobs), "jobs with >1 hash", sum(len(v) > 1 for v in jobs.values()))
l = [p[-1] for p in passes if p[-1] != ""]
print("load", min(l), max(l), "passes without load", sum(p[-1] == "" for p in passes))
