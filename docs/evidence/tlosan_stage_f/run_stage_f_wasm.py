#!/usr/bin/env python3
"""Compare plain serial and real threaded WASI command outputs with native/NCBI."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SERIAL = ROOT / "LOSAT/target/serial-command/wasm32-wasip1/release/LOSAT.wasm"
THREADED = ROOT / "LOSAT/target/threaded-command/wasm32-wasip1-threads/release/LOSAT.wasm"
SERIAL_RUNNER = ROOT / "LOSAT/tests/run_losat_wasi.js"
THREADED_RUNNER = ROOT / "LOSAT/tests/run_losat_wasi_threads.js"
PHYSICAL = (
    "remaining_local_20260925__run_20260925__multi_hsp_20260924.mode2",
    "alternate_matrix_20260925__run_20260925__command_6",
    "natural_positive_20260925__run_20260925__heap_replacement",
    "natural_positive_20260925__containment_run_20260925__command_5",
    "kappa_heap_rejection_20260925__natural_c_d_early_20260925__command_5",
    "kappa_result_order_20260925__multi_query_20260924_default",
    "kappa_result_order_20260925__seg_hard_query_20260924_default",
)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def rows(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines()]


def command_at_threads(command: list[str], n: int) -> list[str]:
    result = command.copy()
    result[result.index("-num_threads") + 1] = str(n)
    return result


def run(command: list[str]) -> subprocess.CompletedProcess[bytes]:
    env = dict(os.environ)
    env["LOSAT_WASI_THREADS_DEBUG"] = "1"
    return subprocess.run(command, cwd=ROOT, capture_output=True, env=env)


# NCBI c++/src/algo/blast/api/prelim_stage.cpp:172-188:
# (*thread)->Run(); (*thread)->Join(&result);
# The WASI runners provide threads/IO only; native and NCBI bytes define output.
def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("physical_matrix", type=Path)
    parser.add_argument("codes_multi_matrix", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    assert SERIAL.is_file() and THREADED.is_file()
    physical = [row for row in rows(args.physical_matrix / "comparison.jsonl")
                if row["threads"] == 1 and row["repetition"] == 0
                and row["case"].removesuffix(f"__fmt{row['outfmt']}") in PHYSICAL]
    codes = [row for row in rows(args.codes_multi_matrix / "comparison.jsonl")
             if row["threads"] == 1 and row["repetition"] == 0]
    assert len(physical) == len(PHYSICAL) * 3 and len(codes) == 27 * 3
    expected = physical + codes
    env = {
        "node_version": subprocess.check_output(["node", "--version"]).decode().strip(),
        "serial_target": "wasm32-wasip1", "serial_feature": "no-default-features",
        "threaded_target": "wasm32-wasip1-threads", "threaded_feature": "wasm-threads",
        "serial_wasm_sha256": sha(SERIAL.read_bytes()),
        "threaded_wasm_sha256": sha(THREADED.read_bytes()),
        "serial_runner": str(SERIAL_RUNNER), "threaded_runner": str(THREADED_RUNNER),
        "serial_runner_sha256": sha(SERIAL_RUNNER.read_bytes()),
        "threaded_runner_sha256": sha(THREADED_RUNNER.read_bytes()),
        "native_executable_sha256": expected[0]["losat_executable_sha256"],
    }
    (out / "environment.json").write_text(json.dumps(env, indent=2) + "\n")
    output = out / "comparison.jsonl"
    with output.open("w") as log:
        for base in expected:
            name, fmt = base["case"], base["outfmt"]
            native_command = command_at_threads(base["losat_command"], 1)
            native = run(native_command)
            assert native.returncode == 0 and sha(native.stdout) == base["expected_sha256"], (name, fmt)
            assert sha(Path(base["query"]).read_bytes()) == base["query_sha256"]
            assert sha(Path(base["subject"]).read_bytes()) == base["subject_sha256"]
            for mode, n in (("serial", 1), ("threaded", 1), ("threaded", 2),
                            ("threaded", 4), ("threaded", 8)):
                repeats = (3 if mode == "threaded" and n > 1 and (
                    "natural_c_d_early" in name or "heap_replacement" in name or "code32_" in name
                ) else 1)
                for repetition in range(repeats):
                    prefix = ["node", "--no-warnings", "--experimental-wasi-unstable-preview1",
                              str(SERIAL_RUNNER if mode == "serial" else THREADED_RUNNER),
                              str(SERIAL if mode == "serial" else THREADED)]
                    command = prefix + command_at_threads(native_command[1:], n)
                    result = run(command)
                    activity = re.search(rb"\[losat-thread-activity\].*peak_active=(\d+)", result.stderr)
                    record = {
                        "case": name, "outfmt": fmt, "mode": mode, "threads": n,
                        "repetition": repetition, "command": command,
                        "query_sha256": base["query_sha256"],
                        "subject_sha256": base["subject_sha256"],
                        "native_command": native_command,
                        "ncbi_contract_sha256": base["expected_sha256"],
                        "native_sha256": sha(native.stdout),
                        "wasm_sha256": sha(result.stdout), "stdout_bytes": len(result.stdout),
                        "stderr_sha256": sha(result.stderr),
                        "peak_active": int(activity.group(1)) if activity else None,
                        "activity_log": activity.group(0).decode() if activity else None,
                        "exit": result.returncode,
                        "equal": result.returncode == 0 and result.stdout == native.stdout,
                    }
                    log.write(json.dumps(record, sort_keys=True) + "\n")
                    log.flush()
                    if not record["equal"]:
                        (out / "first_failure.json").write_text(json.dumps(record, indent=2) + "\n")
                        (out / "first_failure.native").write_bytes(native.stdout)
                        (out / "first_failure.wasm").write_bytes(result.stdout)
                        (out / "first_failure.stderr").write_bytes(result.stderr)
                        raise AssertionError(f"{name} fmt{fmt} {mode} n{n} repeat{repetition}")
            print("PASS", name, f"fmt{fmt}", flush=True)
    negative = run(["node", "--no-warnings", "--experimental-wasi-unstable-preview1",
                    str(SERIAL_RUNNER), str(SERIAL), *command_at_threads(expected[0]["losat_command"][1:], 2)])
    assert negative.returncode != 0 and not negative.stdout
    (out / "plain_wasi_threads_rejected.json").write_text(json.dumps({
        "command": ["node", "--no-warnings", "--experimental-wasi-unstable-preview1",
                    str(SERIAL_RUNNER), str(SERIAL), *command_at_threads(expected[0]["losat_command"][1:], 2)],
        "exit": negative.returncode, "stdout_sha256": sha(negative.stdout),
        "stderr_sha256": sha(negative.stderr), "stderr": negative.stderr.decode(),
    }, indent=2) + "\n")
    result = rows(output)
    (out / "summary.json").write_text(json.dumps({
        "comparison_count": len(result), "equal": sum(r["equal"] for r in result),
        "unique_case_formats": len({(r["case"], r["outfmt"]) for r in result}),
        "by_target_thread": {f"{mode}:{n}": sum(r["mode"] == mode and r["threads"] == n for r in result)
                             for mode, n in (("serial", 1), ("threaded", 1),
                                             ("threaded", 2), ("threaded", 4), ("threaded", 8))},
        "peak_active": {str(n): max((r["peak_active"] or 0 for r in result
                                    if r["mode"] == "threaded" and r["threads"] == n), default=0)
                        for n in (2, 4, 8)},
    }, indent=2) + "\n")
    print(f"{len(result)}/{len(result)} WASI outputs equal native/NCBI contract", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
