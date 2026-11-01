#!/usr/bin/env python3
# usage: summarize_ab.py FILE.jsonl  - per case: median wall per arm, guarded/raw (median of per-round ratios), hash agreement
import json, statistics as st, sys
rows = [json.loads(l) for l in open(sys.argv[1]) if l.strip()]
hdr = next((r["header"] for r in rows if "header" in r), {})
rows = [r for r in rows if "case" in r and r["round"] >= 0]   # round -1 = discarded warm-up
print("header:", {k: hdr.get(k) for k in ("threads", "rounds", "skipped", "guarded_has_record", "raw_has_record", "guard_of_raw_equals_guarded_file")})
order = ["bp", "bpff", "mb", "bn", "dc", "bns", "tbn", "bx", "tbx"]
cases = [c for c in order if any(r["case"] == c for r in rows)]
print(f"{'case':6} {'n(g/r)':7} {'med guarded':>11} {'med raw':>9} {'g/r (med of rounds)':>20} {'guard_s(g)':>10} {'out hashes':>12} {'fail g/r':>9} {'compiled_ok':>11}")
for c in cases:
    d = {a: {r["round"]: r for r in rows if r["case"] == c and r["arm"] == a} for a in ("guarded", "raw")}
    ok = {a: [r for r in d[a].values() if r["exit"] == 0 and not r["error"]] for a in d}
    med = {a: (st.median(r["wall_s"] for r in ok[a]) if ok[a] else float("nan")) for a in d}
    ratios = [d["guarded"][k]["wall_s"] / d["raw"][k]["wall_s"] for k in d["guarded"] if k in d["raw"]
              and d["guarded"][k] in ok["guarded"] and d["raw"][k] in ok["raw"]]
    shas = {r["out_sha256"] for a in ok for r in ok[a]}
    allsha = {r["out_sha256"] for a in d for r in d[a].values()}
    agree = "AGREE" if len(shas) == 1 and None not in shas else ("DIFFER(%d)" % len(shas) if shas else "none")
    fails = (len(d["guarded"]) - len(ok["guarded"]), len(d["raw"]) - len(ok["raw"]))
    gs = st.median(r["guard_s"] for r in ok["guarded"]) if ok["guarded"] else float("nan")
    cok = all(r["compiled_ok"] for a in d for r in d[a].values())
    print(f"{c:6} {len(ok['guarded'])}/{len(ok['raw']):<5} {med['guarded']:11.3f} {med['raw']:9.3f} "
          f"{(st.median(ratios) if ratios else float('nan')):20.4f} {gs:10.4f} {agree:>12} {fails[0]}/{fails[1]:<7} {str(cok):>11}")
bad = [r for r in rows if r["exit"] != 0 or r["error"]]
for r in bad: print("FAILED", r["case"], r["arm"], "round", r["round"], "exit", r["exit"], r["error"][:120])
