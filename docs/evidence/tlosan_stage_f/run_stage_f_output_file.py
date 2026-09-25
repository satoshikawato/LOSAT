#!/usr/bin/env python3
"""Compare TBLASTN output-file bytes and empty stdout for all native thread counts."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "docs/evidence/tlosan_stage_e"))
from run_stage_e_cli_matrix import LOSAT, NCBI, read_commands


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


# NCBI c++/src/algo/blast/format/blast_format.cpp:1411-1458:
# formatter.PrintOneResultSet(...);
# Output-file bytes and empty stdout are checked after formatting.
def main() -> int:
    out = Path(sys.argv[1]).resolve()
    out.mkdir(parents=True, exist_ok=False)
    family, case, base = next((f, c, b) for f, c, b in read_commands()
                              if "natural_c_d_early" in f)
    query = Path(base[base.index("-query") + 1])
    subject = Path(base[base.index("-subject") + 1])
    rows = []
    for fmt in (0, 6, 7):
        expected_file = out / f"ncbi_fmt{fmt}.out"
        oracle_cmd = [str(NCBI), *base[1:], "-outfmt", str(fmt), "-out", str(expected_file)]
        expected = subprocess.run(oracle_cmd, cwd=ROOT, capture_output=True)
        assert expected.returncode == 0 and not expected.stdout
        expected_bytes = expected_file.read_bytes()
        for n in (1, 2, 4, 8):
            actual_file = out / f"losat_fmt{fmt}_n{n}.out"
            common = base[1:].copy()
            common[common.index("-num_threads") + 1] = str(n)
            command = [str(LOSAT), "tblastn", *common, "-outfmt", str(fmt),
                       "-out", str(actual_file)]
            actual = subprocess.run(command, cwd=ROOT, capture_output=True)
            row = {
                "case": f"{family}/{case}", "outfmt": fmt, "threads": n,
                "query_sha256": sha(query.read_bytes()),
                "subject_sha256": sha(subject.read_bytes()),
                "ncbi_command": oracle_cmd, "losat_command": command,
                "ncbi_file_sha256": sha(expected_bytes),
                "losat_file_sha256": sha(actual_file.read_bytes()) if actual_file.exists() else None,
                "ncbi_stdout_sha256": sha(expected.stdout), "losat_stdout_sha256": sha(actual.stdout),
                "ncbi_exit": expected.returncode, "losat_exit": actual.returncode,
                "equal": actual.returncode == 0 and not actual.stdout
                    and actual_file.exists() and actual_file.read_bytes() == expected_bytes,
            }
            rows.append(row)
            if not row["equal"]:
                raise AssertionError(f"outfmt {fmt} threads {n}")
    (out / "comparison.jsonl").write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows))
    print(f"{len(rows)}/{len(rows)} output files byte-identical and stdout empty")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
