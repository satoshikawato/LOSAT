#!/usr/bin/env python3
"""Read-only evidence reconciliation; never replace canonical output.

NCBI 2.17.0 c++/src/objtools/align_format/tabular.cpp:1100-1108:
  x_PrintField(*iter); ... m_Ostream << "\n";
Use the existing complete ordered-byte and exception validators. Each joined
row names its actual process evidence; an earlier timeout is retained on disk.
"""

import json
from collections import Counter
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3] / "LOSAT/tests"))
import wasm_performance as perf
from retry_audit import read_result


def main():
    manifest = json.loads((HERE / "manifest-concurrent.json").read_text())
    perf.preflight(manifest)
    old = HERE.parent / "initial-audit"
    original = [json.loads(s) for s in (old / "samples.jsonl").read_text().splitlines()]
    updates = json.loads((HERE / "concurrent-audit/summary.json").read_text())["rows"]
    update_count = len(updates)
    updates = {(r["program"], r["case_id"], r["kind"]): r for r in updates}
    assert len(updates) == update_count, "duplicate retry selector"
    cases = {(c["program"], c["case_id"]): c for c in manifest["cases"]}
    expected_keys = {(c["program"], c["case_id"], kind)
                     for c in manifest["cases"] if not c.get("extra")
                     for kind in ("native", "serial")
                     if kind == "native" or c["serial_applicable"]}
    original_keys = [(r["program"], r["case_id"], r["kind"]) for r in original]
    assert len(original_keys) == len(set(original_keys)), "duplicate original selector"
    assert set(original_keys) == expected_keys, "missing/extra original selector"
    assert set(updates) == {(r["program"], r["case_id"], r["kind"])
                            for r in original if r["status"] != "PASS"}, "missing/extra retry selector"
    rows = []
    for row in original:
        key = (row["program"], row["case_id"], row["kind"])
        case = cases[key[:2]]
        update = updates.get(key)
        result_path = Path(update["result_source"]) if update else Path(row["output"]).parent / "result.json"
        oracle_path = Path(update["oracle_source"]) if update else old / row["program"] / row["case_id"] / "oracle/result.json"
        result, oracle = read_result(result_path), read_result(oracle_path)
        if update:
            perf.checked_file(result_path, update["result_sha256"])
            perf.checked_file(oracle_path, update["oracle_result_sha256"])
        for path, command in ((result_path, perf.thread_command(case, row["kind"], row["threads_requested"], manifest, exact=True)), (oracle_path, case["oracle_argv"])):
            saved = json.loads((path.parent / "command.json").read_text())
            assert saved["cwd"] == manifest["cwd"] and saved["shell"] is False
            assert perf.authority.replace_output(saved["ordered_argv"], Path("SINK")) == perf.authority.replace_output(command, Path("SINK"))
        canonical = perf.raw_gate(result, case["losat_sha256"])
        behavioral = perf.oracle_contract(case, oracle, result)
        status = canonical if canonical != "PASS" else "PASS" if behavioral == "PASS" else "ORACLE_FAIL"
        if result["peak_rss_bytes"] is not None and result["peak_rss_bytes"] > manifest["measurement"]["memory_budget_bytes"]:
            status = "MEMORY_BUDGET_FAIL"
        rows.append({"program": key[0], "case_id": key[1], "kind": key[2], "status": status,
                     "canonical_gate": canonical, "oracle_contract": behavioral,
                     "original_status": row["status"], "followup_selected": bool(update),
                     "result_source": str(result_path), "oracle_source": str(oracle_path),
                     "raw_output_sha256": result["raw_output_sha256"], "raw_output_bytes": result["raw_output_bytes"]})
    assert len(rows) == 84 and sum(r["kind"] == "native" for r in rows) == 43
    perf.dump(HERE / "combined-audit.json", {
        "method": "Current preflight, saved command binding, raw SHA/length checks and existing behavioral/exception validators; completed evidence reused, incomplete evidence replaced only in this joined view.",
        "counts": dict(Counter(r["status"] for r in rows)),
        "canonical_counts": dict(Counter(r["canonical_gate"] for r in rows)),
        "oracle_contract_counts": dict(Counter(r["oracle_contract"] for r in rows)),
        "rows": rows,
    })
    print(Counter(r["status"] for r in rows))


if __name__ == "__main__":
    main()
