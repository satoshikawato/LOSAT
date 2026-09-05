#!/usr/bin/env python3
"""Recheck retained raw evidence and export PR01 tables without changing gates.

NCBI c++/src/objtools/align_format/tabular.cpp:1098-1108:
  x_PrintField(*iter); ... m_Ostream << "\n";
Hashes here cover those complete ordered bytes; TSV tables summarize evidence
and never replace the saved search output or its acceptance decision.
"""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import statistics

import wasm_performance as perf


def series_result(rows, measurement):
    """Require every ordered raw-output gate, including the diagnostic output.

    NCBI tabular.cpp:1098-1108: x_PrintField(*iter); ... m_Ostream << "\n";
    A successful earlier output cannot stand in for a later failed command.
    """
    failed = next((s for s in rows if s["status"] != "PASS"), None)
    if failed:
        return failed["status"], failed, None
    warmups = measurement["warmup_count"]
    count = measurement["timed_repetitions"]
    if (
        len(rows) != warmups + count + 1
        or [s.get("repeat") for s in rows] != list(range(warmups + count + 1))
        or any(
            s.get("timed") != (warmups <= i < warmups + count)
            or s.get("diagnostics_enabled") != (i == warmups + count)
            for i, s in enumerate(rows)
        )
    ):
        return "NOT_RUN", None, None
    timed = [s for s in rows if s["timed"]]
    vals = [s["wall_seconds"] for s in timed]
    return (
        "PASS",
        rows[0],
        {
            "median_seconds": statistics.median(vals),
            "min_seconds": min(vals),
            "max_seconds": max(vals),
            "samples": count,
            "peak_rss_bytes": max(s["peak_rss_bytes"] for s in timed)
            if all(s["peak_rss_bytes"] is not None for s in timed)
            else None,
        },
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence", type=Path)
    args = parser.parse_args()
    root = args.evidence.resolve()
    manifest = json.loads((root / "contract/manifest.json").read_text())
    perf.preflight(manifest)
    cases = {(c["program"], c["case_id"]): c for c in manifest["cases"]}
    samples, summaries = [], {}
    for name in ("initial-audit", "cold", "warm", "protected-profile"):
        path = root / name / "samples.jsonl"
        if path.exists():
            for line in path.read_text().splitlines():
                row = json.loads(line)
                if row.get("raw_output_sha256") is not None:
                    perf.checked_file(row["output"], row["raw_output_sha256"])
                    if Path(row["output"]).stat().st_size != row["raw_output_bytes"]:
                        raise perf.GateFailure("retained output byte length changed")
                if row["status"] == "PASS":
                    expected = cases[(row["program"], row["case_id"])]["losat_sha256"]
                    if row.get("raw_output_sha256") != expected:
                        raise perf.GateFailure("PASS sample differs from frozen output")
                    if row.get("exit_status", 0) != 0:
                        raise perf.GateFailure("PASS sample has nonzero process exit")
                    if name != "warm":
                        result = json.loads(
                            (Path(row["output"]).parent / "result.json").read_text()
                        )
                        for field in (
                            "status",
                            "exit_status",
                            "wall_seconds",
                            "cpu_user_seconds",
                            "cpu_system_seconds",
                            "peak_rss_bytes",
                        ):
                            if row.get(field) != result.get(field):
                                raise perf.GateFailure(
                                    f"sample differs from process result: {field}"
                                )
                samples.append({"run": name, **row})
        summary = root / name / "summary.json"
        summaries[name] = (
            json.loads(summary.read_text())
            if summary.exists()
            else {"status": "NOT_RUN"}
        )
    fields = [
        "run",
        "program",
        "case_id",
        "kind",
        "threads_requested",
        "effective_compute_threads",
        "spawned_helpers",
        "active_workers",
        "repeat",
        "timed",
        "diagnostics_enabled",
        "boundary",
        "status",
        "canonical_gate",
        "oracle_status",
        "oracle_contract_status",
        "exit_status",
        "started_utc",
        "wall_seconds",
        "cpu_user_seconds",
        "cpu_system_seconds",
        "peak_rss_bytes",
        "rss_after_bytes",
        "process_lifetime_peak_rss_bytes",
        "linear_memory_after_bytes",
        "raw_output_sha256",
        "raw_output_bytes",
        "output",
        "reason",
    ]
    with (root / "samples.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(samples)
    revalidated = []
    for row in samples:
        if row["run"] != "initial-audit":
            continue
        case = cases[(row["program"], row["case_id"])]
        result = json.loads((Path(row["output"]).parent / "result.json").read_text())
        oracle = json.loads(
            (Path(row["output"]).parents[2] / "oracle/result.json").read_text()
        )
        revalidated.append(
            {
                "program": row["program"],
                "case_id": row["case_id"],
                "kind": row["kind"],
                "canonical_gate": perf.raw_gate(result, case["losat_sha256"]),
                "oracle_contract": perf.oracle_contract(case, oracle, result),
            }
        )
    perf.dump(
        root / "audit-revalidation.json",
        {
            "method": "retained raw SHA/length verification plus current existing exception validators; no new searches",
            "rows": revalidated,
        },
    )
    profiles = [
        {
            k: s.get(k)
            for k in (
                "run",
                "program",
                "case_id",
                "kind",
                "threads_requested",
                "effective_compute_threads",
                "spawned_helpers",
                "stage_seconds",
                "thread_stages",
                "output",
                "status",
            )
        }
        for s in samples
        if s.get("diagnostics_enabled")
    ]
    perf.dump(root / "stage-profiles.json", profiles)
    cold_groups = {}
    for sample in samples:
        if sample["run"] == "cold":
            key = "/".join(
                str(sample[k])
                for k in ("program", "case_id", "kind", "threads_requested")
            )
            cold_groups.setdefault(key, []).append(sample)
    cold = {}
    for key, rows in cold_groups.items():
        _, _, median = series_result(rows, manifest["measurement"])
        if median is not None:
            cold[key] = median
    if cold != summaries["cold"].get("medians", {}):
        raise perf.GateFailure("saved cold summary differs from retained sample series")
    comparisons = []
    for case in manifest["cases"]:
        if not case["focus"] or case["focus"].get("benchmark") != "true":
            continue
        prefix = f"{case['program']}/{case['case_id']}"
        native = cold.get(prefix + "/native/1")
        serial = cold.get(prefix + "/serial/1")
        threaded = cold.get(prefix + "/threaded/1")
        row = {
            "program": case["program"],
            "case_id": case["case_id"],
            "native1_seconds": native["median_seconds"] if native else None,
            "serial_seconds": serial["median_seconds"] if serial else None,
            "threaded1_seconds": threaded["median_seconds"] if threaded else None,
            "serial_native_ratio": serial["median_seconds"] / native["median_seconds"]
            if native and serial
            else None,
        }
        for n in (2, 4, 8):
            nw = cold.get(prefix + f"/threaded/{n}")
            nn = cold.get(prefix + f"/native/{n}")
            row[f"requested_{n}_speedup"] = (
                threaded["median_seconds"] / nw["median_seconds"]
                if threaded and nw
                else None
            )
            row[f"requested_{n}_wasm_native_ratio"] = (
                nw["median_seconds"] / nn["median_seconds"] if nw and nn else None
            )
        comparisons.append(row)
    perf.dump(
        root / "comparisons.json",
        {
            "boundary": "cold process",
            "note": "N ratios compare requested N; native effective N and active utilization are unknown. No efficiency or target attainment claim.",
            "rows": comparisons,
        },
    )
    # Fill the supplied per-run template as well as the normalized contract.
    template_text = (root / "run-manifest-template.json").read_text()
    platform_authority = json.loads(
        (perf.TESTS / "ncbi_platform_variance_v010.json").read_text()
    )
    cpu_model = next(
        (
            line.split(":", 1)[1].strip()
            for line in manifest["hardware"]["cpuinfo"].splitlines()
            if line.startswith("model name")
        ),
        None,
    )
    run_records, seen = [], set()
    for sample in samples:
        if sample["run"] != "cold":
            continue
        key = (
            sample["program"],
            sample["case_id"],
            sample["kind"],
            sample["threads_requested"],
        )
        if key in seen:
            continue
        seen.add(key)
        series = [
            s
            for s in samples
            if s["run"] == "cold"
            and (s["program"], s["case_id"], s["kind"], s["threads_requested"]) == key
        ]
        diagnostic = next((s for s in series if s.get("diagnostics_enabled")), {})
        series_status, representative, _ = series_result(
            series, manifest["measurement"]
        )
        representative = representative or {}
        record = json.loads(template_text)
        record.pop("required_before_execution")
        case = cases[key[:2]]
        artifact = manifest["artifacts"][sample["kind"]]
        command = (
            json.loads((Path(sample["output"]).parent / "command.json").read_text())
            if sample.get("output")
            else None
        )
        case_dir = root / "cold" / sample["program"] / sample["case_id"]
        oracle_command = json.loads((case_dir / "oracle/command.json").read_text())
        inputs = {
            role: path
            for role, _, path in perf.authority._input_arguments(
                case["losat_argv"], sample["program"]
            )
        }
        policy_path = perf.TESTS / (
            {
                "blastn": "blastn_parity_manifest.tsv",
                "blastp": "blastp_v010_parity_manifest.tsv",
                "tblastx": "tblastx_v010_parity_manifest.tsv",
            }[sample["program"]]
            if not case.get("extra")
            else "wasm_performance_extra_cases.json"
        )
        record.update(
            {
                "status": "RECORDED_LOCAL_SCOPE_NOT_CERTIFIED",
                "pr_id": "PR01",
                "run_id": "/".join(map(str, key)),
                "baseline_or_candidate": "baseline",
            }
        )
        record["source"].update(
            {
                "head": manifest["source"]["head"],
                "dirty_identity": manifest["source"]["dirty_identity"],
                "source_snapshot_sha256": perf.digest(root / "source-identity.json"),
            }
        )
        record["build"].update(
            {
                "rustc_verbose": manifest["toolchain"]["rustc_vV"],
                "cargo_version": manifest["toolchain"]["cargo_V"],
                "cargo_lock_sha256": artifact["build"]["cargo_lock_sha256"],
                "target": artifact["build"]["target"],
                "features": artifact["build"]["features"],
                "rustflags": artifact["build"]["rustflags"],
                "artifact_path": artifact["path"],
                "artifact_sha256": artifact["sha256"],
                "expected_exports": artifact.get("exports", []),
                "expected_imports": artifact.get("imports", []),
            }
        )
        record["hardware"].update(
            {
                "cpu_model": cpu_model,
                "logical_cpus": manifest["hardware"]["logical_cpus"],
                "memory_bytes": next(
                    (
                        int(line.split()[1]) * 1024
                        for line in manifest["hardware"]["meminfo"].splitlines()
                        if line.startswith("MemTotal:")
                    ),
                    None,
                ),
                "os": manifest["hardware"]["os"],
            }
        )
        runner = manifest["runtime"]["runners"].get(sample["kind"])
        record["runtime"].update(
            {
                "name": "native" if sample["kind"] == "native" else "Node.js",
                "version": None
                if sample["kind"] == "native"
                else manifest["runtime"]["version"],
                "argv_prefix": []
                if runner is None
                else [
                    manifest["runtime"]["node"],
                    *manifest["runtime"]["flags"],
                    runner,
                ],
                "flags": [] if runner is None else manifest["runtime"]["flags"],
                "host_source_sha256": manifest["runtime"]["runner_hashes"].get(runner),
            }
        )
        record["fixture"].update(
            {
                "id": sample["case_id"],
                "program": sample["program"],
                "task": case["losat_argv"][case["losat_argv"].index("--task") + 1]
                if "--task" in case["losat_argv"]
                else sample["program"],
                "query_path": inputs["query"],
                "query_sha256": manifest["inputs"][inputs["query"]]["sha256"],
                "subject_path": inputs["subject"],
                "subject_sha256": manifest["inputs"][inputs["subject"]]["sha256"],
                "manifest_path": str(policy_path),
                "manifest_sha256": manifest["authorities"][str(policy_path)],
                "contract_class": case["classification"],
                "expected_raw_sha256": case["losat_sha256"],
            }
        )
        record["command"].update(
            {
                "ordered_argv": command["ordered_argv"] if command else None,
                "cwd": manifest["cwd"],
                "controlled_lexical_paths": list(inputs.values()),
                "environment_overrides": command["environment_overrides"]
                if command
                else {},
                "output_sink": sample.get("output"),
            }
        )
        record["threads"].update(
            {
                "requested": sample["threads_requested"],
                "host_cap": manifest["measurement"]["host_cap"],
                "effective_compute": diagnostic.get(
                    "effective_compute_threads", sample.get("effective_compute_threads")
                ),
                "spawned_helpers": diagnostic.get("spawned_helpers"),
                "observation_source": diagnostic.get("output"),
                "observation_scope": "Separate diagnostic; configured threads, not active utilization",
            }
        )
        record["oracle"].update(
            {
                "official_version": manifest["oracle"][sample["program"]]["version"],
                "source_sha256": manifest["oracle"]["official_source_archive"][
                    "sha256"
                ],
                "executable_sha256": manifest["oracle"][sample["program"]]["sha256"],
                "ordered_argv": oracle_command["ordered_argv"],
                "platform_authority_version": platform_authority["authority_version"],
                "platform_authority_sha256": manifest["authorities"][
                    str(perf.TESTS / "ncbi_platform_variance_v010.json")
                ],
            }
        )
        record["measurement"].update(
            {
                "boundary": sample["boundary"],
                "sample_file": str(root / "samples.tsv"),
                "eligible_case_set_sha256": manifest["authorities"][
                    str(perf.TESTS / "wasm_performance_cases.tsv")
                ],
                "memory_budget_bytes": manifest["measurement"]["memory_budget_bytes"],
            }
        )
        record["result"].update(
            {
                "exit_status": representative.get("exit_status"),
                "raw_output_sha256": representative.get("raw_output_sha256"),
                "raw_output_bytes": representative.get("raw_output_bytes"),
                "parity_status": series_status,
                "result_source": representative.get("output"),
                "scope": "All required warmup, timed and diagnostic rows; hash refers to result_source",
                "performance_status": "INCONCLUSIVE",
                "missing_reason": representative.get("reason")
                or ("incomplete series" if series_status == "NOT_RUN" else None),
            }
        )
        record["unmeasured_field_reasons"] = (
            "Native effective N/active workers, power mode, physical core topology and temperature not independently measured. Frozen canonical authority supplies SHA only; expected byte count is not learned from candidate. Final full-scope audit not accepted. See EVIDENCE.md."
        )
        run_records.append(record)
    with (root / "RUN_MANIFEST.jsonl").open("w") as handle:
        for record in run_records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    lines = [
        "# PR01 measured results",
        "",
        "Status: NOT_PROVEN. These are baseline measurements, not optimization or release certification.",
        "",
        "| Run | Status | Recorded rows |",
        "| --- | --- | ---: |",
    ]
    for name, summary in summaries.items():
        lines.append(
            f"| {name} | {summary['status']} | {sum(s['run'] == name for s in samples)} |"
        )
    lines += [
        "",
        "Initial audit raw Gate A revalidation: "
        + str(dict(Counter(r["canonical_gate"] for r in revalidated)))
        + ".",
        "Initial audit oracle-contract revalidation: "
        + str(dict(Counter(r["oracle_contract"] for r in revalidated)))
        + ".",
        "",
        "Cold medians below each require one successful warmup, five successful timed outputs, and the separate diagnostic gate. Missing cells are NOT_RUN. Times are seconds.",
        "",
        "| Case | Native 1 | Plain serial Wasm | Threaded Wasm 1 | Serial/native | Requested 4 speedup |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]

    def number(value):
        return "NOT_RUN" if value is None else f"{value:.3f}"

    for row in comparisons:
        lines.append(
            "| "
            + row["case_id"]
            + " | "
            + " | ".join(
                number(row[k])
                for k in (
                    "native1_seconds",
                    "serial_seconds",
                    "threaded1_seconds",
                    "serial_native_ratio",
                    "requested_4_speedup",
                )
            )
            + " |"
        )
    lines += [
        "",
        "Complete median ranges and RSS: [cold/summary.json](cold/summary.json). All 2/4/8 requested-N comparisons: [comparisons.json](comparisons.json). Warm serial results: [warm/summary.json](warm/summary.json).",
        "",
        "[samples.tsv](samples.tsv) retains failures, warmups and diagnostics as separate rows. [stage-profiles.json](stage-profiles.json) links diagnostic output and logs. Smoke runs are excluded from these tables. Stage timers are nested and may sum work across threads; they are not an additive wall-time partition.",
        "",
        "No whole-eligible-set geometric mean or performance attainment is accepted while an eligible case is missing. Requested-N speedup is not evidence that N computation threads were active. Platform Gate B, real browsers, resident API, threaded warm lifecycle, and native active-worker measurements remain NOT_RUN.",
    ]
    (root / "RESULTS.md").write_text("\n".join(lines) + "\n")
    print(root / "RESULTS.md")


if __name__ == "__main__":
    main()
