#!/usr/bin/env python3
"""Compare Stage G real and no-hit command-WASI bytes with native/NCBI."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
SERIAL = ROOT / "LOSAT/target/serial-command/wasm32-wasip1/release/LOSAT.wasm"
THREADED = ROOT / "LOSAT/target/threaded-command/wasm32-wasip1-threads/release/LOSAT.wasm"
SERIAL_RUNNER = ROOT / "LOSAT/tests/run_losat_wasi.js"
THREADED_RUNNER = ROOT / "LOSAT/tests/run_losat_wasi_threads.js"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> int:
    real, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    out.mkdir(parents=True, exist_ok=False)
    rows = [json.loads(line) for line in (real / "comparison.jsonl").read_text().splitlines()]
    native = [row for row in rows if row["threads"] == 1]
    assert len(native) == 18 and all(row["equal"] for row in native)
    environment = {
        "serial_wasm_sha256": sha(SERIAL.read_bytes()),
        "threaded_wasm_sha256": sha(THREADED.read_bytes()),
        "native_sha256": sha(Path(native[0]["losat_command"][0]).read_bytes()),
        "node_version": subprocess.check_output(["node", "--version"], text=True).strip(),
        "real_matrix": str(real),
        "real_matrix_sha256": sha((real / "comparison.jsonl").read_bytes()),
    }
    (out / "environment.json").write_text(json.dumps(environment, indent=2) + "\n")
    count = 0
    with (out / "comparison.jsonl").open("w") as log:
        for base in native:
            expected = (real / f"{base['case']}.fmt{base['outfmt']}.oracle").read_bytes()
            assert sha(expected) == base["expected_sha256"]
            assert sha(Path(base["query"]).read_bytes()) == base["query_sha256"]
            assert sha(Path(base["subject"]).read_bytes()) == base["subject_sha256"]
            assert base["actual_stderr_sha256"] == base["code1_cli_stderr_sha256"]
            # NCBI c++/src/algo/blast/api/prelim_stage.cpp:172-188:
            # (*thread)->Run(); (*thread)->Join(&result);
            # The hosts supply only command-WASI threads and I/O; compare raw bytes.
            for mode, n in (("serial", 1), ("threaded", 1), ("threaded", 2),
                            ("threaded", 4), ("threaded", 8)):
                args = base["losat_command"][1:]
                args[args.index("-num_threads") + 1] = str(n)
                command = ["node", "--no-warnings", "--experimental-wasi-unstable-preview1",
                           str(SERIAL_RUNNER if mode == "serial" else THREADED_RUNNER),
                           str(SERIAL if mode == "serial" else THREADED), *args]
                result = subprocess.run(command, cwd=ROOT, capture_output=True)
                row = {
                    "case": base["case"], "code": base["code"], "outfmt": base["outfmt"],
                    "mode": mode, "threads": n, "command": command,
                    "expected_sha256": sha(expected), "actual_sha256": sha(result.stdout),
                    "stdout_bytes": len(result.stdout), "stderr_sha256": sha(result.stderr),
                    "expected_stderr_sha256": base["code1_cli_stderr_sha256"],
                    "stderr_equal": sha(result.stderr) == base["code1_cli_stderr_sha256"],
                    "exit": result.returncode,
                    # NCBI c++/src/algo/blast/format/blast_format.cpp:1443-1451:
                    # formatter emits search warning strings on the error stream.
                    "equal": result.returncode == 0 and result.stdout == expected
                             and sha(result.stderr) == base["code1_cli_stderr_sha256"],
                }
                log.write(json.dumps(row, sort_keys=True) + "\n")
                log.flush()
                if not row["equal"]:
                    (out / "first_failure.json").write_text(json.dumps(row, indent=2) + "\n")
                    (out / "first_failure.expected").write_bytes(expected)
                    (out / "first_failure.actual").write_bytes(result.stdout)
                    (out / "first_failure.stderr").write_bytes(result.stderr)
                    raise AssertionError((base["case"], base["outfmt"], mode, n))
                count += 1
            print("PASS", base["case"], f"fmt{base['outfmt']}", flush=True)
    assert count == 90
    (out / "summary.json").write_text(json.dumps({"comparisons": count, "equal": count}, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
