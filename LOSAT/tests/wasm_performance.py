#!/usr/bin/env python3
"""PR01 explicit-artifact baseline; frozen hashes are never learned from LOSAT.

NCBI 2.17.0 c++/src/algo/blast/format/blast_format.cpp:794-808,828-832:
  tabinfo.PrintHeader(strProgVersion, ..., dbname, ...);
  ITERATE(CSeq_align_set::Tdata, itr, copy_aln_set.Get()) { ... tabinfo.Print(); }
c++/src/objtools/align_format/tabular.cpp:1098-1108:
  x_PrintField(*iter); ... m_Ostream << "\n";
This test-only orchestration preserves ordered bytes and lexical input paths.
It supplies process/runtime measurement, never search behavior or a fallback.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import statistics
import subprocess
import sys
import threading
import time

import certify_platform_native_v010 as authority
import certify_integrated_runtime_v010 as integrated
import certify_blastn_v010 as blastn_cert
import audit_tblastx_v010 as tblastx_audit

ROOT = Path(__file__).resolve().parents[2]
TESTS = ROOT / "LOSAT/tests"
SCHEMA = "losat-wasm-perf-run-v1"
ORACLE_HASHES = {
    "blastn": "33b64bc67d3149cee2459b2f7766b363323df632cf12c099546de00aea9698b5",
    "blastp": "5ce267c04e4988c265357bfbedc64e809545b6fcfae7ff6775266fabbee8ba0e",
    "tblastx": "583e5d60bbd444ac455d20e0956c5aa0aeef675da8daee8204d8f9376ddb8804",
}


class GateFailure(RuntimeError):
    pass


def digest(path):
    return integrated.sha256_path(Path(path))


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def checked_file(path, expected):
    path = Path(path)
    if not path.is_file():
        raise GateFailure(f"missing artifact/authority/input: {path}")
    if digest(path) != expected:
        raise GateFailure(f"SHA256 mismatch: {path}")


def capture(argv):
    return subprocess.check_output(
        argv, cwd=ROOT, text=True, stderr=subprocess.STDOUT
    ).strip()


def wasm_identity(path, node):
    script = "const f=require('fs');const m=new WebAssembly.Module(f.readFileSync(process.argv[1]));console.log(JSON.stringify({imports:WebAssembly.Module.imports(m),exports:WebAssembly.Module.exports(m)}));"
    return json.loads(capture([node, "-e", script, str(path)]))


def validate_kind(kind, inspection):
    exports = {x["name"] for x in inspection["exports"]}
    imports = {(x["module"], x["name"]) for x in inspection["imports"]}
    if "_start" not in exports:
        raise GateFailure("command artifact missing _start")
    threaded = ("wasi", "thread-spawn") in imports and "wasi_thread_start" in exports
    if kind == "threaded" and not threaded:
        raise GateFailure(
            "threaded label on artifact without thread-spawn/wasi_thread_start"
        )
    if kind == "serial" and (threaded or any(m == "wasi" for m, _ in imports)):
        raise GateFailure("serial label on threaded artifact")


def read_cases(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def freeze(args):
    dest = args.output.resolve()
    dest.mkdir(parents=True, exist_ok=False)
    head = capture(["git", "rev-parse", "HEAD"])
    catalog = authority.load_catalog(ROOT, head)
    canonical = ROOT / authority.CANONICAL_MANIFEST_REPO_RELATIVE
    # Binding the worktree file to the committed authority is independent of B/C.
    if canonical.read_bytes() != subprocess.check_output(
        ["git", "show", f"{head}:{canonical.relative_to(ROOT)}"], cwd=ROOT
    ):
        raise GateFailure("canonical worktree authority mismatch")
    oracles = {p: args.oracle_dir.resolve() / p for p in ORACLE_HASHES}
    specs = [
        integrated.OracleSpec(
            p, path, ORACLE_HASHES[p], (f"{p}: 2.17.0+", "Package: blast 2.17.0")
        )
        for p, path in oracles.items()
    ]
    oracle_info = integrated.validate_oracles(
        ROOT,
        dest,
        specs,
        args.source_archive.resolve(),
        "502057a88e9990e34e62758be21ea474cc0ad68d6a63a2e37b2372af1e5ea147",
    )
    artifacts = {}
    for kind in ("native", "serial", "threaded"):
        path = getattr(args, kind).resolve()
        if not path.is_file():
            raise GateFailure(f"missing explicit {kind} artifact: {path}")
        item = {"path": str(path), "sha256": digest(path), "kind": kind}
        if kind != "native":
            item.update(wasm_identity(path, args.node))
            validate_kind(kind, item)
        elif path.read_bytes()[:4] != b"\x7fELF":
            raise GateFailure(
                "this measured host requires native ELF; not command Wasm"
            )
        artifacts[kind] = item
    toolchain = integrated.record_toolchain(ROOT, dest)
    for kind, artifact in artifacts.items():
        target = {
            "native": toolchain["host_triple"],
            "serial": "wasm32-wasip1",
            "threaded": "wasm32-wasip1-threads",
        }[kind]
        target_dir = Path(artifact["path"]).parents[1 if kind == "native" else 2]
        build_argv = ["cargo", "build", "--locked", "--release", "--bin", "LOSAT"]
        if kind != "native":
            build_argv += ["--target", target]
            build_argv += (
                ["--no-default-features"]
                if kind == "serial"
                else ["--features", "wasm-threads"]
            )
        build_argv += ["--target-dir", str(target_dir)]
        artifact["build"] = {
            "target": target,
            "features": []
            if kind == "serial"
            else ["parallel"] + (["wasm-threads"] if kind == "threaded" else []),
            "profile": "release",
            "ordered_argv": build_argv,
            "cwd": str(ROOT / "LOSAT"),
            "cargo_lock_sha256": digest(ROOT / "LOSAT/Cargo.lock"),
            "rustflags": [] if kind == "native" else ["-C", "target-feature=+simd128"],
        }
        key = {
            "native": "native",
            "serial": "serial_wasm",
            "threaded": "threaded_wasm",
        }[kind]
        toolchain[key] = artifact["build"]
    dump(dest / "cert_toolchain.json", toolchain)
    focuses = read_cases(args.cases)
    focus = {(x["program"], x["case_id"]): x for x in focuses}
    if not set(focus) <= set(catalog.canonical):
        raise GateFailure("focus case missing from frozen canonical")
    cases, inputs = [], {}
    for row in catalog.canonical_rows:
        program, case_id = row["program"], row["case_id"]
        ncbi, losat = authority._base_commands(
            ROOT, catalog, program, case_id, args.native.resolve(), oracles, dest
        )
        for _, _, lexical in authority._input_arguments(losat, program):
            relative = authority._fixture_relative_from_lexical(lexical)
            source = ROOT / "LOSAT" / str(relative)
            data = source.read_bytes()
            target = Path(lexical)
            # Never overwrite another run's controlled fixture or normalize FASTA.
            if target.exists() and target.read_bytes() != data:
                raise GateFailure(f"controlled input conflict: {target}")
            if not target.exists():
                target.parent.mkdir(parents=True, exist_ok=True)
                with target.open("xb") as handle:
                    handle.write(data)
            inputs[lexical] = {
                "source": str(source),
                "sha256": digest(target),
                "bytes": len(data),
                "sequences": sum(line.startswith(b">") for line in data.splitlines()),
            }
        cases.append(
            {
                **row,
                "losat_argv": losat,
                "oracle_argv": ncbi,
                "serial_applicable": program != "blastp"
                or catalog.blastp_cases[case_id].num_threads == 1,
                "focus": focus.get((program, case_id)),
            }
        )
    if len(cases) != 43 or sum(c["serial_applicable"] for c in cases) != 41:
        raise GateFailure("43/41 contract mismatch")
    extra_path = TESTS / "wasm_performance_extra_cases.json"
    for extra in json.loads(extra_path.read_text())["cases"]:
        base = next(
            c
            for c in cases
            if c["program"] == extra["program"]
            and c["case_id"] == extra["base_case_id"]
        )
        item = {
            **base,
            "case_id": extra["case_id"],
            "losat_sha256": extra["expected_raw_sha256"],
            "extra": True,
            "focus": {
                "role": extra["role"],
                "gap_eligible": "false",
                "protected": "true",
                "benchmark": "true",
            },
        }
        item["losat_argv"] = [*base["losat_argv"], "--evalue", extra["evalue"]]
        item["oracle_argv"] = [*base["oracle_argv"], "-evalue", extra["evalue"]]
        if extra.get("input_relative"):
            lexical = authority.historical_lexical_path(extra["input_relative"])
            if lexical not in inputs:
                raise GateFailure("extra input must already be staged and hashed")
            for field in ("losat_argv", "oracle_argv"):
                for _, option, _ in authority._input_arguments(
                    item[field], extra["program"]
                ):
                    item[field][item[field].index(option) + 1] = lexical
        cases.append(item)
    authorities = {
        str(p): digest(p)
        for p in [
            canonical,
            TESTS / "ncbi_platform_variance_v010.json",
            TESTS / "blastn_parity_manifest.tsv",
            TESTS / "blastp_v010_parity_manifest.tsv",
            TESTS / "tblastx_v010_parity_manifest.tsv",
            args.cases.resolve(),
        ]
    }
    authorities.update(
        {
            str(p): digest(p)
            for p in [extra_path, TESTS / "blastn_v010_source_exceptions.tsv"]
        }
    )
    sources = {str(p): digest(p) for p in sorted((ROOT / "LOSAT/src").rglob("*.rs"))}
    sources.update(
        {
            str(p): digest(p)
            for p in [
                ROOT / "LOSAT/Cargo.toml",
                ROOT / "LOSAT/Cargo.lock",
                ROOT / "LOSAT/.cargo/config.toml",
            ]
        }
    )
    runners = {
        k: str(TESTS / name)
        for k, name in [
            ("serial", "run_losat_wasi.js"),
            ("threaded", "run_losat_wasi_threads.js"),
        ]
    }
    runners["warm_serial"] = str(TESTS / "run_wasm_performance_warm.js")
    cpus = os.cpu_count()
    manifest = {
        "schema_version": SCHEMA,
        "pr_id": "PR01",
        "baseline_or_candidate": "baseline",
        "status": "FROZEN_NOT_CERTIFIED",
        "source": {
            "head": head,
            "files": sources,
            "dirty_identity": hashlib.sha256(
                json.dumps(sources, sort_keys=True).encode()
            ).hexdigest(),
        },
        "toolchain": toolchain,
        "artifacts": artifacts,
        "oracle": oracle_info,
        "authorities": authorities,
        "runtime": {
            "node": shutil.which(args.node),
            "version": capture([args.node, "--version"]),
            "flags": [],
            "runners": runners,
            "runner_hashes": {p: digest(p) for p in runners.values()},
        },
        "hardware": {
            "os": platform.platform(),
            "cpuinfo": Path("/proc/cpuinfo").read_text(),
            "meminfo": Path("/proc/meminfo").read_text(),
            "logical_cpus": cpus,
            "affinity": sorted(os.sched_getaffinity(0)),
            "power_mode": None,
            "temperature": None,
        },
        "inputs": inputs,
        "cases": cases,
        "cwd": str(ROOT),
        "measurement": {
            "boundary": "cold-process-through-teardown",
            "warmup_count": 1,
            "timed_repetitions": 5,
            "thread_matrix": [1, 2, 4, 8],
            "host_cap": min(8, len(os.sched_getaffinity(0))),
            "memory_budget_bytes": 4 * 1024**3,
            "wasm_linear_memory_limit_bytes": 1024**3,
            "timeout_seconds": args.timeout,
            "diagnostics_enabled": False,
            "targets": {
                "geomean_wasm_native": 1.5,
                "per_case_wasm_native": 2.0,
                "speedup4": 2.0,
                "speedup8": 3.0,
                "protected_candidate_baseline": 1.05,
            },
            "parallel_eligibility": "baseline engine >=1s and proven independent jobs >=8 (N4), >=16 (N8); otherwise NOT_PROVEN",
            "warm_module": "NOT_RUN: requires separate host boundary",
            "resident": "NOT_RUN: PR03 lifecycle proof required",
            "browser": "NOT_RUN: Chromium/Firefox real host deferred to PR12",
        },
    }
    manifest["measurement"]["wasm_linear_memory_limit_bytes"] = {
        "threaded": 1024**3,
        "serial": None,
    }
    manifest["measurement"]["memory_limit_note"] = (
        "Threaded host caps linear memory; serial maximum is module-defined and not enforced by this harness. RSS budget is a post-run rejection threshold."
    )
    manifest["measurement"]["warm_module"] = (
        "Separate serial cached-module/new-instance series, warmup 1 + repeats 5; returnOnExit=true; threaded NOT_RUN pending host lifecycle proof"
    )
    manifest["measurement"]["protected_case_policy"] = (
        "All 43/41 initial contracts; focus benchmark=false rows receive profile/matrix gates without repeated timing. Timing targets include gap_eligible=true only, selected before candidate optimization."
    )
    manifest["runtime"]["node_sha256"] = digest(manifest["runtime"]["node"])
    dump(dest / "manifest.json", manifest)
    print(dest / "manifest.json", flush=True)


def environment(manifest, diagnostic=False):
    env = {
        k: v
        for k, v in os.environ.items()
        if not k.startswith(("LOSAT_", "RAYON_", "NODE_"))
    }
    env.update(
        {
            "NODE_NO_WARNINGS": "1",
            "LOSAT_WASI_THREAD_CAP": str(manifest["measurement"]["host_cap"]),
            "LOSAT_WASM_MEMORY_MAXIMUM_PAGES": "16384",
        }
    )
    if diagnostic:
        env.update({"LOSAT_TIMING": "1", "LOSAT_WASI_THREADS_DEBUG": "1"})
    return env


def execute(argv, cwd, directory, env, timeout):
    directory.mkdir(parents=True, exist_ok=False)
    output = directory / "output.txt"
    argv = authority.replace_output(argv, output)
    dump(
        directory / "command.json",
        {
            "ordered_argv": argv,
            "cwd": str(cwd),
            "shell": False,
            "sink": str(output),
            "environment_overrides": {
                k: v
                for k, v in env.items()
                if k.startswith(("LOSAT_", "RAYON_", "NODE_"))
            },
        },
    )
    usage = directory / "usage.txt"
    measured = ["/usr/bin/time", "-f", "%U\t%S\t%M", "-o", str(usage), *argv]
    start = time.monotonic()
    started_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    status, exit_code, reason = "PASS", None, None
    with (
        (directory / "stdout.txt").open("wb") as stdout,
        (directory / "stderr.txt").open("wb") as stderr,
    ):
        try:
            process = subprocess.Popen(
                measured,
                cwd=cwd,
                env=env,
                stdout=stdout,
                stderr=stderr,
                start_new_session=True,
            )
            expired = threading.Event()

            def expire():
                if process.poll() is None:
                    expired.set()
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass

            watchdog = threading.Timer(timeout, expire)
            watchdog.start()
            try:
                exit_code = process.wait()  # Blocking waitpid: no timeout polling bias.
            finally:
                watchdog.cancel()
                watchdog.join()
            if expired.is_set():
                status, reason = "TIMEOUT", f"exceeded {timeout}s"
            if status == "PASS" and exit_code != 0:
                status, reason = "BAD_EXIT", f"subprocess exit {exit_code}"
        except OSError as error:
            status, reason = "BAD_EXIT", str(error)
    elapsed = time.monotonic() - start
    if status == "PASS" and not output.is_file():
        status, reason = "MISSING_OUTPUT", "successful subprocess did not write output"
    user = system = rss = None
    if usage.exists():
        fields = (
            usage.read_text().strip().splitlines()[-1].split("\t")
            if usage.read_text().strip()
            else []
        )
        if len(fields) == 3:
            user, system, rss = (
                float(fields[0]),
                float(fields[1]),
                int(fields[2]) * 1024,
            )
    result = {
        "status": status,
        "reason": reason,
        "exit_status": exit_code,
        "wall_seconds": elapsed,
        "cpu_user_seconds": user,
        "cpu_system_seconds": system,
        "peak_rss_bytes": rss,
        "output": str(output),
        "raw_output_sha256": digest(output) if output.is_file() else None,
        "raw_output_bytes": output.stat().st_size if output.is_file() else None,
    }
    result["started_utc"] = started_utc
    dump(directory / "result.json", result)
    return result


# NCBI tabular.cpp:1098-1108: m_Ostream << m_FieldDelimiter; ... << "\n";
# No whitespace/header/order normalization is permissible in this acceptance gate.
def raw_gate(result, expected):
    if result["status"] != "PASS":
        return result["status"]
    return "PASS" if result["raw_output_sha256"] == expected else "BASELINE_PARITY_FAIL"


def oracle_contract(case, oracle, result):
    if oracle["status"] != "PASS":
        return oracle["status"]
    if result["status"] != "PASS":
        return "NOT_RUN"
    ncbi, losat = Path(oracle["output"]), Path(result["output"])
    if case["classification"] == "EXACT_TEXT":
        return "PASS" if ncbi.read_bytes() == losat.read_bytes() else "ORACLE_FAIL"
    if case["classification"] == "SOURCE_UNDETERMINED_ACCEPTED":
        try:
            exceptions = blastn_cert.read_source_exceptions(
                TESTS / "blastn_v010_source_exceptions.tsv"
            )
            blastn_cert.certify_source_exception(
                exceptions[case["case_id"]], ncbi, losat
            )
            return "PASS"
        except blastn_cert.CertificationError:
            return "ORACLE_FAIL"
    cases = {
        c.case_id: c
        for c in tblastx_audit.load_manifest(
            TESTS / "tblastx_v010_parity_manifest.tsv", ROOT
        )
    }
    classification, _ = tblastx_audit.classify_outputs(ncbi, losat)
    return (
        "PASS"
        if classification == "HSP_SET_DIFF"
        and tblastx_audit.contract_accepts(cases[case["case_id"]], classification)
        else "ORACLE_FAIL"
    )


def diagnostic_diff(expected, observed, directory):
    import difflib

    left, right = Path(expected).read_bytes(), Path(observed).read_bytes()
    a, b = (
        left.decode(errors="replace").splitlines(True),
        right.decode(errors="replace").splitlines(True),
    )
    Path(directory, "ordered.diff").write_text(
        "".join(difflib.unified_diff(a, b, fromfile="oracle", tofile="LOSAT"))
    )
    Path(directory, "sorted-diagnostic.diff").write_text(
        "".join(difflib.unified_diff(sorted(a), sorted(b)))
    )


def thread_command(case, kind, threads, manifest, exact=False):
    argv = list(case["losat_argv"])
    argv[0] = manifest["artifacts"][kind]["path"]
    options = ("--num-threads", "--num_threads", "-num_threads", "-n")
    option = next((x for x in options if x in argv), None)
    if not exact:
        if option is None:
            raise GateFailure(f"thread option not found: {argv}")
        argv[argv.index(option) + 1] = str(threads)
    if kind != "native":
        argv = [
            manifest["runtime"]["node"],
            *manifest["runtime"]["flags"],
            manifest["runtime"]["runners"][kind],
            *argv,
        ]
    return argv


def preflight(manifest):
    if manifest["schema_version"] != SCHEMA or manifest["cwd"] != str(ROOT):
        raise GateFailure("schema/cwd contract mismatch")
    catalog = authority.load_catalog(ROOT, manifest["source"]["head"])
    primary = [c for c in manifest["cases"] if not c.get("extra")]
    if len(primary) != 43 or sum(c["serial_applicable"] for c in primary) != 41:
        raise GateFailure("manifest lost 43/41 contracts")
    if {(c["program"], c["case_id"]) for c in primary} != set(catalog.canonical):
        raise GateFailure("manifest canonical selectors changed")
    for case in primary:
        expected = catalog.canonical[(case["program"], case["case_id"])]
        if any(
            case[k] != expected[k]
            for k in ("losat_sha256", "classification", "contract")
        ):
            raise GateFailure(
                "manifest authority mismatch; expected cannot come from candidate"
            )
        applicable = (
            case["program"] != "blastp"
            or catalog.blastp_cases[case["case_id"]].num_threads == 1
        )
        if case["serial_applicable"] != applicable:
            raise GateFailure("serial applicability mismatch")
        oracle, losat = authority._base_commands(
            ROOT,
            catalog,
            case["program"],
            case["case_id"],
            Path(manifest["artifacts"]["native"]["path"]),
            {p: Path(manifest["oracle"][p]["path"]) for p in ORACLE_HASHES},
            Path("unused"),
        )
        for field, command in (("oracle_argv", oracle), ("losat_argv", losat)):
            if authority.replace_output(
                case[field], Path("SINK")
            ) != authority.replace_output(command, Path("SINK")):
                raise GateFailure(f"ordered argv authority mismatch: {case['case_id']}")
    focus = {
        (c["program"], c["case_id"]): c
        for c in read_cases(TESTS / "wasm_performance_cases.tsv")
    }
    for case in primary:
        if case["focus"] != focus.get((case["program"], case["case_id"])):
            raise GateFailure("eligible/protected case policy mismatch")
    extras = {
        c["case_id"]: c
        for c in json.loads((TESTS / "wasm_performance_extra_cases.json").read_text())[
            "cases"
        ]
    }
    extra_cases = [c for c in manifest["cases"] if c.get("extra")]
    if len(extra_cases) != len(extras) or {c["case_id"] for c in extra_cases} != set(
        extras
    ):
        raise GateFailure("extra case catalog mismatch")
    for case in manifest["cases"]:
        if (
            case.get("extra")
            and case["losat_sha256"] != extras[case["case_id"]]["expected_raw_sha256"]
        ):
            raise GateFailure("extra no-hit authority mismatch")
        if case.get("extra"):
            extra = extras[case["case_id"]]
            base = next(
                c
                for c in primary
                if c["program"] == extra["program"]
                and c["case_id"] == extra["base_case_id"]
            )
            if any(
                case[k] != base[k]
                for k in ("program", "classification", "contract", "serial_applicable")
            ):
                raise GateFailure("extra case contract mismatch")
            if case["focus"] != {
                "role": extra["role"],
                "gap_eligible": "false",
                "protected": "true",
                "benchmark": "true",
            }:
                raise GateFailure("extra case measurement policy mismatch")
            for field, option in (
                ("losat_argv", "--evalue"),
                ("oracle_argv", "-evalue"),
            ):
                command = [*base[field], option, extra["evalue"]]
                if extra.get("input_relative"):
                    for _, flag, _ in authority._input_arguments(
                        command, extra["program"]
                    ):
                        command[command.index(flag) + 1] = (
                            authority.historical_lexical_path(extra["input_relative"])
                        )
                if command != case[field]:
                    raise GateFailure("extra ordered argv authority mismatch")
    for group in (
        manifest["authorities"],
        manifest["source"]["files"],
        manifest["runtime"]["runner_hashes"],
    ):
        for path, expected in group.items():
            checked_file(path, expected)
    checked_file(manifest["runtime"]["node"], manifest["runtime"]["node_sha256"])
    for path, item in manifest["inputs"].items():
        checked_file(path, item["sha256"])
        checked_file(item["source"], item["sha256"])
    for kind, item in manifest["artifacts"].items():
        checked_file(item["path"], item["sha256"])
        if kind != "native":
            validate_kind(
                kind, wasm_identity(item["path"], manifest["runtime"]["node"])
            )
    for p in ORACLE_HASHES:
        checked_file(manifest["oracle"][p]["path"], ORACLE_HASHES[p])


def run(args):
    if args.phase == "audit" and args.threads:
        raise GateFailure(
            "audit preserves canonical thread arguments; no thread override"
        )
    manifest = json.loads(args.manifest.read_text())
    dest = args.output.resolve()
    dest.mkdir(parents=True, exist_ok=False)
    preflight(manifest)
    dump(dest / "manifest.json", manifest)
    dump(
        dest / "invocation.json",
        {
            "run_id": dest.name,
            "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "argv": sys.argv,
            "manifest_sha256": digest(args.manifest),
            "phase": args.phase,
            "harness_sha256": digest(__file__),
        },
    )
    if args.phase == "warm":
        return run_warm(args, manifest, dest)
    samples = []

    def save_sample(sample):
        samples.append(sample)
        with (dest / "samples.jsonl").open("a") as handle:
            handle.write(json.dumps(sample, sort_keys=True) + "\n")

    timeout = manifest["measurement"]["timeout_seconds"]
    selected = [
        c
        for c in manifest["cases"]
        if (
            not c.get("extra")
            if args.phase == "audit"
            else c["focus"]
            and (args.phase == "profile" or c["focus"]["benchmark"] == "true")
        )
    ]
    if args.case:
        selected = [c for c in selected if c["case_id"] in args.case]
        if not selected or set(args.case) != {c["case_id"] for c in selected}:
            raise GateFailure("unknown selected case")
    for case in selected:
        p, cid = case["program"], case["case_id"]
        case_dir = dest / p / cid
        print(f"{args.phase}: {p}/{cid}", flush=True)
        oracle = execute(
            case["oracle_argv"],
            manifest["cwd"],
            case_dir / "oracle",
            environment(manifest),
            timeout,
        )
        # Orthogonal platform Gate B is deliberately never inferred from this run.
        oracle_status = oracle["status"]
        exact_oracle = case["classification"] == "EXACT_TEXT"
        if exact_oracle:
            oracle_status = raw_gate(oracle, case["losat_sha256"])
        kinds = args.targets or (
            ["native", "serial"]
            if args.phase == "audit"
            else ["native", "serial", "threaded"]
        )
        for kind in kinds:
            if kind == "serial" and not case["serial_applicable"]:
                continue
            counts = (
                [1]
                if args.phase == "audit" or kind == "serial"
                else args.threads or manifest["measurement"]["thread_matrix"]
            )
            for threads in counts:
                entry = {
                    "program": p,
                    "case_id": cid,
                    "kind": kind,
                    "threads_requested": threads,
                    "effective_compute_threads": 1 if kind == "serial" else None,
                    "host_cap": manifest["measurement"]["host_cap"],
                    "spawned_helpers": None,
                    "active_workers": None,
                    "oracle_status": oracle_status,
                    "platform_gate_b": "NOT_RUN",
                    "platform_reason": "Linux local baseline is outside the registered Windows/macOS Gate B",
                    "contract_class": case["classification"],
                    "boundary": manifest["measurement"]["boundary"],
                }
                if args.phase != "audit" and oracle_status != "PASS":
                    save_sample(
                        {
                            **entry,
                            "status": "NOT_RUN",
                            "reason": f"oracle gate {oracle_status}; timing forbidden",
                        }
                    )
                    continue
                if threads > manifest["measurement"]["host_cap"]:
                    save_sample(
                        {**entry, "status": "NOT_RUN", "reason": "hardware cap"}
                    )
                    continue
                command = thread_command(
                    case, kind, threads, manifest, exact=args.phase == "audit"
                )
                if args.phase == "audit":
                    option = next(
                        x
                        for x in (
                            "--num-threads",
                            "--num_threads",
                            "-num_threads",
                            "-n",
                        )
                        if x in command
                    )
                    entry["threads_requested"] = int(command[command.index(option) + 1])
                warmups = manifest["measurement"]["warmup_count"]
                timed = manifest["measurement"]["timed_repetitions"]
                if warmups < 1 or timed < 5:
                    raise GateFailure(
                        "benchmark needs at least one warmup and five timed repeats"
                    )
                repeats = (
                    1
                    if args.phase == "audit"
                    else (2 if args.phase == "profile" else warmups + timed + 1)
                )
                for i in range(repeats):
                    diagnostic = args.phase != "audit" and i == repeats - 1
                    label = "diagnostic" if diagnostic else f"repeat-{i}"
                    result = execute(
                        command,
                        manifest["cwd"],
                        case_dir / f"{kind}-{threads}" / label,
                        environment(manifest, diagnostic),
                        timeout,
                    )
                    status = raw_gate(result, case["losat_sha256"])
                    behavioral_status = oracle_contract(case, oracle, result)
                    if status == "PASS" and (
                        oracle_status != "PASS" or behavioral_status != "PASS"
                    ):
                        status = "ORACLE_FAIL"
                    if (
                        result["peak_rss_bytes"] is not None
                        and result["peak_rss_bytes"]
                        > manifest["measurement"]["memory_budget_bytes"]
                    ):
                        status = "MEMORY_BUDGET_FAIL"
                    sample = {
                        **entry,
                        **result,
                        "status": status,
                        "oracle_contract_status": behavioral_status,
                        "canonical_gate": raw_gate(result, case["losat_sha256"]),
                        "repeat": i,
                        "timed": args.phase == "benchmark"
                        and warmups <= i < warmups + timed,
                        "diagnostics_enabled": diagnostic,
                    }
                    if diagnostic:
                        log = Path(result["output"]).parent / "stderr.txt"
                        logtext = log.read_text(errors="replace")
                        sample["stage_seconds"] = {
                            k: float(v)
                            for k, v in re.findall(
                                r"\[TIMING\] ([\w_]+): ([\d.]+)s", logtext
                            )
                        }
                        alignment = re.search(
                            r"\[TIMING\] traceback_alignment: .*\(dp=([\d.]+)s greedy=([\d.]+)s\)",
                            logtext,
                        )
                        if alignment:
                            sample["stage_seconds"].update(
                                {
                                    "traceback_alignment_dp": float(alignment[1]),
                                    "traceback_alignment_greedy": float(alignment[2]),
                                }
                            )
                        sample["thread_diagnostics"] = re.findall(
                            r"[^\n]*(?:effective_threads=|spawn tid=)[^\n]*", logtext
                        )
                        sample["spawned_helpers"] = (
                            len(re.findall(r"spawn tid=", logtext))
                            if kind == "threaded"
                            else None
                        )
                        stages = []
                        for line in sample["thread_diagnostics"]:
                            match = re.search(r"effective_threads=(\d+)", line)
                            if match:
                                stages.append(
                                    {
                                        "log": line,
                                        "effective_compute_threads": 1
                                        if "parallel=false" in line
                                        else int(match[1]),
                                    }
                                )
                        sample["thread_stages"] = stages
                        if stages:
                            sample["effective_compute_threads"] = max(
                                s["effective_compute_threads"] for s in stages
                            )
                    save_sample(sample)
                    if status != "PASS":
                        if oracle["status"] == "PASS" and result["raw_output_sha256"]:
                            diagnostic_diff(
                                oracle["output"],
                                result["output"],
                                Path(result["output"]).parent,
                            )
                        break  # Never benchmark a failed gate or failed repeat.
    groups = {}
    for s in samples:
        if s.get("timed") and s["status"] == "PASS":
            key = f"{s['program']}/{s['case_id']}/{s['kind']}/{s['threads_requested']}"
            groups.setdefault(key, []).append(s)
    medians = {}
    for key, rows in groups.items():
        group_failed = any(
            f"{s['program']}/{s['case_id']}/{s['kind']}/{s['threads_requested']}" == key
            and s["status"] != "PASS"
            for s in samples
        )
        if (
            len(rows) == manifest["measurement"]["timed_repetitions"]
            and not group_failed
        ):
            vals = [s["wall_seconds"] for s in rows]
            medians[key] = {
                "median_seconds": statistics.median(vals),
                "min_seconds": min(vals),
                "max_seconds": max(vals),
                "samples": len(vals),
                "peak_rss_bytes": max(s["peak_rss_bytes"] for s in rows)
                if all(s["peak_rss_bytes"] is not None for s in rows)
                else None,
            }
    failed = [s for s in samples if s["status"] != "PASS"]
    dump(
        dest / "summary.json",
        {
            "schema_version": SCHEMA,
            "status": "FAIL" if failed else "PASS_EXECUTED_SCOPE",
            "phase": args.phase,
            "selected_cases": len(selected),
            "sample_count": len(samples),
            "failures": failed,
            "medians": medians,
            "release_certification": "NOT_RUN",
            "final_broad_gate": "NOT_RUN",
        },
    )
    print(f"samples={len(samples)} failures={len(failed)}", flush=True)
    return 1 if failed else 0


def run_warm(args, manifest, dest):
    if args.targets and args.targets != ["serial"]:
        raise GateFailure("warm boundary currently supports serial only")
    if args.threads and args.threads != [1]:
        raise GateFailure("warm serial boundary requires one thread")
    selected = [
        c for c in manifest["cases"] if c["focus"] and c["focus"]["benchmark"] == "true"
    ]
    if args.case:
        selected = [c for c in selected if c["case_id"] in args.case]
        if {c["case_id"] for c in selected} != set(args.case):
            raise GateFailure("unknown warm case")
    rows, medians = [], {}
    warmups, repeats = (
        manifest["measurement"]["warmup_count"],
        manifest["measurement"]["timed_repetitions"],
    )
    if warmups < 1 or repeats < 5:
        raise GateFailure(
            "warm benchmark needs at least one warmup and five timed repeats"
        )
    for case in selected:
        directory = dest / case["program"] / case["case_id"]
        directory.mkdir(parents=True)
        print(f"warm: {case['program']}/{case['case_id']}", flush=True)
        oracle = execute(
            case["oracle_argv"],
            manifest["cwd"],
            directory / "oracle",
            environment(manifest),
            manifest["measurement"]["timeout_seconds"],
        )
        if raw_gate(oracle, case["losat_sha256"]) != "PASS":
            sample = {
                "program": case["program"],
                "case_id": case["case_id"],
                "status": "NOT_RUN",
                "reason": f"warm oracle gate {oracle['status']} failed",
                "timed": False,
            }
            rows.append(sample)
            with (dest / "samples.jsonl").open("a") as handle:
                handle.write(json.dumps(sample) + "\n")
            continue
        jobs = []
        for i in range(warmups + repeats):
            output = directory / f"repeat-{i}.txt"
            argv = thread_command(case, "serial", 1, manifest)[
                2 + len(manifest["runtime"]["flags"]) :
            ]
            # Strip host/runner, retaining Wasm argv[0] and ordered CLI args.
            argv = authority.replace_output(argv, output)
            jobs.append(
                {
                    "repeat": i,
                    "argv": argv,
                    "output": str(output),
                    "expected_sha256": case["losat_sha256"],
                }
            )
        dump(directory / "jobs.json", jobs)
        command = [
            manifest["runtime"]["node"],
            *manifest["runtime"]["flags"],
            manifest["runtime"]["runners"]["warm_serial"],
            manifest["artifacts"]["serial"]["path"],
            str(directory / "jobs.json"),
            "--out",
            "unused",
        ]
        execution = execute(
            command,
            manifest["cwd"],
            directory / "host",
            environment(manifest),
            manifest["measurement"]["timeout_seconds"] * len(jobs),
        )
        if execution["status"] != "PASS":
            raise GateFailure(f"warm host failed: {execution}")
        measured = json.loads(Path(execution["output"]).read_text())
        if len(measured["samples"]) != len(jobs):
            raise GateFailure("incomplete warm repeats")
        for sample in measured["samples"]:
            sample.update(
                {
                    "program": case["program"],
                    "case_id": case["case_id"],
                    "kind": "serial",
                    "threads_requested": 1,
                    "effective_compute_threads": 1,
                    "raw_output_sha256": digest(sample["output"]),
                    "raw_output_bytes": Path(sample["output"]).stat().st_size,
                    "status": "PASS" if sample["exit_status"] == 0 else "BAD_EXIT",
                    "timed": sample["repeat"] >= warmups,
                    "diagnostics_enabled": False,
                    "compile_seconds": measured["compile_seconds"],
                }
            )
            sample["status"] = raw_gate(sample, case["losat_sha256"])
            if (
                sample["process_lifetime_peak_rss_bytes"]
                > manifest["measurement"]["memory_budget_bytes"]
            ):
                sample["status"] = "MEMORY_BUDGET_FAIL"
            rows.append(sample)
            with (dest / "samples.jsonl").open("a") as handle:
                handle.write(json.dumps(sample, sort_keys=True) + "\n")
        current = rows[-len(jobs) :]
        if all(s["status"] == "PASS" for s in current):
            vals = [s["wall_seconds"] for s in current if s["timed"]]
            medians[f"{case['program']}/{case['case_id']}"] = {
                "median_seconds": statistics.median(vals),
                "min_seconds": min(vals),
                "max_seconds": max(vals),
                "samples": len(vals),
            }
    failed = [r for r in rows if r["status"] != "PASS"]
    dump(
        dest / "summary.json",
        {
            "status": "FAIL" if failed else "PASS_EXECUTED_SCOPE",
            "phase": "warm",
            "medians": medians,
            "failures": failed,
            "threaded_warm": "NOT_RUN",
            "resident": "NOT_RUN",
        },
    )
    return 1 if failed else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    f = sub.add_parser("freeze")
    for kind in ("native", "serial", "threaded"):
        f.add_argument(f"--{kind}", type=Path, required=True)
    f.add_argument("--oracle-dir", type=Path, required=True)
    f.add_argument("--source-archive", type=Path, required=True)
    f.add_argument("--node", default="node")
    f.set_defaults(cases=TESTS / "wasm_performance_cases.tsv")
    f.add_argument("--timeout", type=float, default=120)
    f.add_argument("--output", type=Path, required=True)
    r = sub.add_parser("run")
    r.add_argument("--manifest", type=Path, required=True)
    r.add_argument("--output", type=Path, required=True)
    r.add_argument(
        "--phase", choices=["audit", "benchmark", "profile", "warm"], required=True
    )
    r.add_argument("--targets", nargs="+", choices=["native", "serial", "threaded"])
    r.add_argument("--case", action="append")
    r.add_argument(
        "--threads",
        type=int,
        nargs="+",
        choices=[1, 2, 4, 8],
        help="explicit subset; omitted cells remain unexecuted",
    )
    args = parser.parse_args()
    try:
        return freeze(args) if args.action == "freeze" else run(args)
    except (
        GateFailure,
        authority.CertificationFailure,
        integrated.CertificationFailure,
        OSError,
        ValueError,
        KeyError,
        TypeError,
        subprocess.SubprocessError,
    ) as error:
        if args.output.is_dir():
            dump(args.output / "failure.json", {"status": "FAIL", "error": str(error)})
        print(f"FAIL: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
