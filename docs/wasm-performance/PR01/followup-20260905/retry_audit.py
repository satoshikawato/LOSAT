#!/usr/bin/env python3
"""Retry incomplete PR01 audit executions; retain successful evidence by reference.

NCBI 2.17.0 c++/src/algo/blast/format/blast_format.cpp:794-808:
  tabinfo.PrintHeader(strProgVersion, ..., dbname, ...);
c++/src/objtools/align_format/tabular.cpp:1100-1108:
  x_PrintField(*iter); ... m_Ostream << "\n";
The existing harness preserves these ordered bytes and lexical paths. This
test-only driver changes execution budget, never search or output authority.
"""

import json
from collections import Counter
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "LOSAT/tests"))
import wasm_performance as perf


def read_result(path):
    result = json.loads(path.read_text())
    if result["raw_output_sha256"] is not None:
        perf.checked_file(result["output"], result["raw_output_sha256"])
        assert Path(result["output"]).stat().st_size == result["raw_output_bytes"]
    return result


def main():
    manifest = json.loads((HERE / "manifest.json").read_text())
    perf.preflight(manifest)
    dest = HERE / "retry-audit"
    dest.mkdir(exist_ok=False)
    old = HERE.parent / "initial-audit"
    rows = [json.loads(line) for line in (old / "samples.jsonl").read_text().splitlines()]
    cases = {c["case_id"]: c for c in manifest["cases"]}
    incomplete = {r["case_id"] for r in rows if r["status"] != "PASS"}
    order = sorted(incomplete, key=lambda c: (0 if c.startswith("p11") else 1 if c.startswith("p10") else 2, c))
    perf.dump(dest / "invocation.json", {
        "argv": sys.argv, "manifest_sha256": perf.digest(HERE / "manifest.json"),
        "driver_sha256": perf.digest(__file__), "harness_sha256": perf.digest(perf.__file__),
        "selected_cases": order, "timing_eligible": False,
    })
    completed = []
    for cid in order:
        case = cases[cid]
        oracle_path = old / case["program"] / cid / "oracle/result.json"
        oracle = read_result(oracle_path)
        if oracle["status"] != "PASS":
            print(f"retry oracle {cid}", flush=True)
            directory = dest / cid / "oracle"
            oracle = perf.execute(case["oracle_argv"], manifest["cwd"], directory,
                                  perf.environment(manifest), manifest["measurement"]["timeout_seconds"])
            oracle_path = directory / "result.json"
            print(f"oracle {cid}: {oracle['status']} {oracle['wall_seconds']:.3f}s", flush=True)
        for row in [r for r in rows if r["case_id"] == cid]:
            source = Path(row["output"]).parent / "result.json"
            result = read_result(source)
            reused = result["status"] == "PASS"
            if not reused:
                print(f"retry {cid} {row['kind']}", flush=True)
                command = perf.thread_command(case, row["kind"], row["threads_requested"], manifest, exact=True)
                directory = dest / cid / row["kind"]
                result = perf.execute(command, manifest["cwd"], directory,
                                      perf.environment(manifest), manifest["measurement"]["timeout_seconds"])
                source = directory / "result.json"
            canonical = perf.raw_gate(result, case["losat_sha256"])
            behavioral = perf.oracle_contract(case, oracle, result)
            status = canonical if canonical != "PASS" else "PASS" if behavioral == "PASS" else "ORACLE_FAIL"
            if result["peak_rss_bytes"] is not None and result["peak_rss_bytes"] > manifest["measurement"]["memory_budget_bytes"]:
                status = "MEMORY_BUDGET_FAIL"
            record = {"case_id": cid, "program": case["program"], "kind": row["kind"],
                      "status": status, "canonical_gate": canonical, "oracle_contract": behavioral,
                      "result_source": str(source), "result_sha256": perf.digest(source),
                      "oracle_source": str(oracle_path), "oracle_result_sha256": perf.digest(oracle_path),
                      "reused_losat": reused, "timing_eligible": False}
            completed.append(record)
            with (dest / "coverage.jsonl").open("a") as handle:
                handle.write(json.dumps(record, sort_keys=True) + "\n")
            print(f"{cid} {row['kind']}: {status}; raw={canonical} oracle={behavioral}", flush=True)
    perf.dump(dest / "summary.json", {"counts": dict(Counter(r["status"] for r in completed)),
              "rows": completed, "complete_catalog": False, "timing_eligible": False})


if __name__ == "__main__":
    main()
