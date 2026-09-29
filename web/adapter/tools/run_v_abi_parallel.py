#!/usr/bin/env python3
"""Runs web/adapter/tests/v_abi.js over a case list in parallel processes.

The full V-ABI suite contains genome-scale TBLASTX searches that take many minutes in
one serial reactor, so the list is split: one part per engine and program, and one
part per TBLASTX search; the longest parts start first. Each part is one `v_abi.js`
process (it registers its own subjects and reuses them across its searches). The
results of all parts are merged into DIR/v-abi-results.json and summarized.

Usage: run_v_abi_parallel.py --native LOSAT --serial SERIAL.wasm --threads THREADS.wasm
           --cases CASES.json --out DIR [--jobs 12]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ADAPTER = Path(__file__).resolve().parents[1]
SIZE_OPTIONS = ("-query", "-subject")


def input_size(search: dict) -> int:
    argv, total = search["argv"], 0
    for option in SIZE_OPTIONS:
        total += (Path(search["cwd"]) / argv[argv.index(option) + 1]).stat().st_size
    return total


def parts(searches: list[dict]) -> list[tuple[str, list[dict]]]:
    grouped: dict[str, list[dict]] = {}
    for index, search in enumerate(searches):
        key = f"{search['program']}-{index}" if search["program"] == "tblastx" else search["program"]
        grouped.setdefault(key, []).append(search)
    return sorted(grouped.items(), key=lambda item: -sum(map(input_size, item[1])))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ("native", "serial", "threads", "cases", "out"):
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument("--jobs", type=int, default=12)
    args = parser.parse_args()
    out = Path(args.out)
    work = [(engine, key, chunk)
            for key, chunk in parts(json.loads(Path(args.cases).read_text()))
            for engine in ("threads", "serial")]

    def run(item: tuple[str, str, list[dict]]) -> tuple[str, int]:
        engine, key, chunk = item
        directory = out / "parts" / f"{engine}-{key}"
        directory.mkdir(parents=True, exist_ok=True)
        # A result left by an earlier run must not be merged into this one.
        (directory / "v-abi-results.json").unlink(missing_ok=True)
        (directory / "cases.json").write_text(json.dumps(chunk, indent=2) + "\n")
        with (directory / "v-abi.log").open("w") as log:
            status = subprocess.run(
                ["node", str(ADAPTER / "tests/v_abi.js"), "--native", args.native, "--serial", args.serial,
                 "--threads", args.threads, "--cases", str(directory / "cases.json"), "--engines", engine,
                 "--out", str(directory)],
                stdout=log, stderr=subprocess.STDOUT).returncode
        print(f"{'ok' if status == 0 else 'FAILED'} {engine}-{key}", flush=True)
        return f"{engine}-{key}", status

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        statuses = list(pool.map(run, work))
    failed = [name for name, status in statuses if status != 0]
    results = []
    for engine, key, _ in work:
        path = out / "parts" / f"{engine}-{key}" / "v-abi-results.json"
        if path.exists():
            results.extend(json.loads(path.read_text()))
    (out / "v-abi-results.json").write_text(json.dumps(results, indent=2) + "\n")
    frozen = [(row, fmt, same) for row in results for fmt, same in row["frozen"].items()]
    differing = sorted({f"{','.join(row['cases'])} outfmt {fmt}" for row, fmt, same in frozen if not same})
    summary = {
        "artifacts": {name: hashlib.sha256(Path(getattr(args, name)).read_bytes()).hexdigest()
                      for name in ("native", "serial", "threads")},
        "runs": len(results),
        "searches": len({json.dumps(row["argv"]) for row in results}),
        "by_engine_and_threads": {f"{e} n{t}": sum(1 for row in results if row["engine"] == e and row["threads"] == t)
                                  for e, t in sorted({(row["engine"], row["threads"]) for row in results})},
        "frozen_equal": sum(1 for _, _, same in frozen if same),
        "frozen_compared": len(frozen),
        "frozen_differing": differing,
        "failed_parts": failed,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
