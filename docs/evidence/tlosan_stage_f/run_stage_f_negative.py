#!/usr/bin/env python3
"""Verify Stage E unsupported paths still fail after Stage F enables threads."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "docs/evidence/tlosan_stage_e"))
from run_stage_e_negative import CASES, LOSAT, QUERY, SUBJECT

# NCBI c++/src/algo/blast/blastinput/blast_args.cpp:3152-3187:
# arg_desc.SetConstraint(kArgNumThreads, new CArgAllowValuesGreaterThanOrEqual(1));
# Stage F replaces Stage E's threads=2 rejection with valid parallel search.
OPTIONS = {name: value for name, value in CASES.items() if name != "threads"}
OPTIONS["zero_threads"] = ["-num_threads", "0"]


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> int:
    out = Path(sys.argv[1]).resolve()
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    for name, option in OPTIONS.items():
        target = out / f"{name}.out"
        cmd = [str(LOSAT), "tblastn", "-query", str(QUERY), "-subject", str(SUBJECT),
               "-out", str(target), *option]
        result = subprocess.run(cmd, cwd=ROOT, capture_output=True)
        row = {
            "case": name, "command": cmd, "exit": result.returncode,
            "stdout_sha256": sha(result.stdout), "stderr_sha256": sha(result.stderr),
            "output_file_created": target.exists(),
            "pass": result.returncode != 0 and not result.stdout and not target.exists(),
        }
        rows.append(row)
    (out / "comparison.jsonl").write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows))
    assert len(rows) == 14 and all(row["pass"] for row in rows)
    print(f"{len(rows)}/{len(rows)} unsupported paths rejected")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
