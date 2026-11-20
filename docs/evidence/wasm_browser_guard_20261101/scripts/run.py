"""WG step 1: run the captured gbdraw jobs through gbdraw's threaded host in Chromium, one pass per
call, and count traps (a "running the pair again" warning). Adapted from gbdraw's analysis/ov2369/run.py.
usage: run.py <arm> <passes> <jobs file under data/> [port]"""
import sys, json, time
from pathlib import Path
from playwright.sync_api import sync_playwright
H = Path(__file__).parent
arm, passes, jobs = sys.argv[1], int(sys.argv[2]), sys.argv[3]
port = int(sys.argv[4]) if len(sys.argv) > 4 else 9893
out = open(H / f"logs/run-{arm}.jsonl", "a")
with sync_playwright() as p:
    b = p.chromium.launch(headless=True)
    page = b.new_page()
    page.on("crash", lambda *_: print("PAGE CRASH", flush=True))
    page.goto(f"http://127.0.0.1:{port}/arm/{arm}/harness.html")
    page.wait_for_function("() => window.__ready", timeout=60000)
    print("chromium", b.version, "arm", arm, flush=True)
    for i in range(passes):
        t = time.time()
        r = page.evaluate("(o) => window.runPass(o)", {"dataUrl": f"data/{jobs}", "totalThreadBudget": "safe", "execution": "threaded"})
        traps = [w for w in r["warns"] if "running the pair again" in w]
        sizes = sorted(p["bytes"] for p in r["probes"] if p.get("kind") == "size")
        probes = [p for p in r["probes"] if p.get("kind") != "size"]
        rec = {"arm": arm, "pass": i, "chromium": b.version, "s": round(time.time() - t, 1), "err": r["err"], "traps": len(traps), "warns": r["warns"],
               "probes": probes, "sizes": sizes, "out": r["out"]}
        out.write(json.dumps(rec) + "\n"); out.flush()
        print(f"pass {i} s={rec['s']} traps={len(traps)} probes={len(probes)} maxMiB={max(sizes or [0])/2**20:.0f} err={r['err']!r:.200}", flush=True)
    b.close()
