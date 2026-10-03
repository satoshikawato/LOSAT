#!/usr/bin/env python3
"""Show that each TBLASTX regression case with an environment depends on it.

Runs `LOSAT tblastx` for every case of LOSAT/tests/fixtures/tblastx_regression/manifest.tsv
that sets `env`, once with the case's environment and once without it, and compares both
with the frozen NCBI output. A case discriminates its environment when the run without it
differs from the frozen output. The cases of NO_EFFECT show that a variable changes
nothing there, so the run without it must be the same.

Usage: env_discrimination.py --losat BIN [--out TSV]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "LOSAT/tests"))
import tblastx_regression_fixtures as fixtures  # noqa: E402

NO_EFFECT = {
    # tblastx_app.cpp:132-137: "Query is Empty!" ends the run before BATCH_SIZE is read.
    "env.batch_text_empty_query": "the empty query is checked before BATCH_SIZE is read",
    "env.batch0_empty_query_fmt0": "the empty query is checked before BATCH_SIZE is read",
    # showdefline.cpp kBits: only the outfmt 0 description table reads it.
    "ctoolkit.fmt7": "CTOOLKIT_COMPATIBLE changes only outfmt 0",
    "ncbirc.harmless_fmt0": "the .ncbirc entries change no output",
    # split_query_cxx.cpp:55-61: the splitter reads the sizes but an ungapped search is
    # never split; 9 and 2^64 - 1 (-1) are divisible by 3, and the overlap is any int.
    "env.chunk9_fmt7": "a chunk size divisible by 3 changes nothing in an ungapped search",
    "env.chunk_minus1_fmt6": "a chunk size divisible by 3 changes nothing in an ungapped search",
    "env.overlap_negative_fmt6": "the overlap changes nothing in an ungapped search",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--losat", required=True)
    parser.add_argument("--out")
    args = parser.parse_args()
    losat = str(Path(args.losat).resolve())
    lines = ["case_id\tenv\twith_env\twithout_env"]
    for row in fixtures.read_manifest():
        if not row["env"]:
            continue
        frozen = (fixtures.FIXTURES / f"{row['case_id']}.out").read_bytes()
        err_path = fixtures.FIXTURES / f"{row['case_id']}.err"
        frozen_err = err_path.read_bytes() if err_path.exists() else b""
        verdicts = []
        for case_env in (row["env"], ""):
            result = fixtures.run_case(losat, ["tblastx"], row["case_id"], row["argv"], case_env, row["losat_extra"])
            same = (result.stdout == frozen and result.stderr == frozen_err
                    and str(result.returncode) == row["exit"])
            verdicts.append("same" if same else "differs")
        lines.append(f"{row['case_id']}\t{row['env']}\t{verdicts[0]}\t{verdicts[1]}")
    text = "\n".join(lines) + "\n"
    if args.out:
        Path(args.out).write_text(text)
    print(text, end="")
    rows = [line.split("\t") for line in lines[1:]]
    bad = [r[0] for r in rows if r[2] != "same" or r[3] != ("same" if r[0] in NO_EFFECT else "differs")]
    print(f"{len(rows)} cases with an environment ({len(NO_EFFECT)} without effect by design), "
          f"{len(bad)} unexpected: {' '.join(bad)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
