#!/usr/bin/env python3
"""Bounded correctness-only audit retries; no concurrent latency benchmarks.

NCBI 2.17.0 c++/src/objtools/align_format/tabular.cpp:1100-1108:
  x_PrintField(*iter); ... m_Ostream << "\n";
Independent oracle/LOSAT processes retain the same complete ordered outputs.
Concurrency here is test scheduling and never a LOSAT search implementation.
"""

from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import Counter
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3] / "LOSAT/tests"))
import wasm_performance as perf
from retry_audit import read_result


def main():
    mp = HERE / "manifest-concurrent.json"
    manifest = json.loads(mp.read_text())
    perf.preflight(manifest)
    dest = HERE / "concurrent-audit"
    dest.mkdir(exist_ok=False)
    old = HERE.parent / "initial-audit"
    rows = [json.loads(s) for s in (old / "samples.jsonl").read_text().splitlines()]
    failed = [r for r in rows if r["status"] != "PASS"]
    cases = {c["case_id"]: c for c in manifest["cases"]}
    selected = sorted({r["case_id"] for r in failed}, key=lambda c: (0 if c.startswith("p10") else 1 if c.startswith("p11") else 2, c))
    paths, tasks = {}, []
    for cid in selected:
        case = cases[cid]
        path = old / case["program"] / cid / "oracle/result.json"
        if cid == "p11_avclpv_psclpv":
            path = HERE / "retry-audit" / cid / "oracle/result.json"
        oracle = read_result(path)
        if oracle["status"] == "PASS" or cid == "p11_avclpv_psclpv":
            paths[(cid, "oracle")] = path
        else:
            tasks.append((cid, "oracle", case["oracle_argv"]))
        for row in [r for r in failed if r["case_id"] == cid]:
            path = Path(row["output"]).parent / "result.json"
            result = read_result(path)
            if result["status"] == "PASS":
                paths[(cid, row["kind"])] = path
            else:
                tasks.append((cid, row["kind"], perf.thread_command(case, row["kind"], row["threads_requested"], manifest, exact=True)))
    assert len(tasks) == manifest["followup"]["maximum_new_searches"]
    perf.dump(dest / "invocation.json", {"argv": sys.argv, "manifest_sha256": perf.digest(mp),
              "driver_sha256": perf.digest(__file__), "harness_sha256": perf.digest(perf.__file__),
              "max_concurrent_searches": 3, "timing_eligible": False,
              "tasks": [{"case_id": c, "kind": k, "ordered_argv": a} for c, k, a in tasks]})

    def execute(task):
        cid, kind, argv = task
        print(f"start {cid} {kind}", flush=True)
        directory = dest / cid / kind
        result = perf.execute(argv, manifest["cwd"], directory, perf.environment(manifest), 600)
        print(f"done {cid} {kind}: {result['status']} {result['wall_seconds']:.3f}s", flush=True)
        return (cid, kind), directory / "result.json"

    with ThreadPoolExecutor(max_workers=3) as pool:
        for future in as_completed([pool.submit(execute, task) for task in tasks]):
            key, path = future.result()
            paths[key] = path
    completed = []
    for row in failed:
        cid, kind = row["case_id"], row["kind"]
        rp, op = paths[(cid, kind)], paths[(cid, "oracle")]
        result, oracle = read_result(rp), read_result(op)
        canonical = perf.raw_gate(result, cases[cid]["losat_sha256"])
        behavioral = perf.oracle_contract(cases[cid], oracle, result)
        status = canonical if canonical != "PASS" else "PASS" if behavioral == "PASS" else "ORACLE_FAIL"
        if result["peak_rss_bytes"] is not None and result["peak_rss_bytes"] > manifest["measurement"]["memory_budget_bytes"]:
            status = "MEMORY_BUDGET_FAIL"
        completed.append({"program": row["program"], "case_id": cid, "kind": kind,
                          "status": status, "canonical_gate": canonical, "oracle_contract": behavioral,
                          "result_source": str(rp), "result_sha256": perf.digest(rp),
                          "oracle_source": str(op), "oracle_result_sha256": perf.digest(op),
                          "reused_losat": old in rp.parents, "timing_eligible": False})
    perf.dump(dest / "summary.json", {"counts": dict(Counter(r["status"] for r in completed)),
              "rows": completed, "complete_catalog": False, "timing_eligible": False})


if __name__ == "__main__":
    main()
