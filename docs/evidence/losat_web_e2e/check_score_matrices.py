#!/usr/bin/env python3
"""Compare LOSAT's protein score matrices (LOSAT/src/utils/matrix.rs) with NCBI's util/tables/sm_*.c.

Each NCBI matrix is a 25x25 `TNCBIScore` array over "ARNDCQEGHILKMFPSTWYVBJZX*" plus its
`defscore`; LOSAT keeps BLOSUM62 as the `BLOSUM62` const and the others as `*_TEXT` blocks
(same symbol order) with the `defscore` in `standard_protein_matrix`.

Usage: check_score_matrices.py --ncbi-src C++DIR [--repo REPO]
Exit 0 when all nine matrices and default scores are equal.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

MATRICES = [("BLOSUM45", "blosum45", "Blosum45"), ("BLOSUM50", "blosum50", "Blosum50"),
            ("BLOSUM62", "blosum62", "Blosum62"), ("BLOSUM80", "blosum80", "Blosum80"),
            ("BLOSUM90", "blosum90", "Blosum90"), ("PAM30", "pam30", "Pam30"),
            ("PAM70", "pam70", "Pam70"), ("PAM250", "pam250", "Pam250"),
            ("IDENTITY", "identity", "Identity")]


def numbers(text: str) -> list[int]:
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    return [int(x) for x in re.findall(r"-?\d+", text)]


def ncbi_matrix(src: Path, cname: str) -> tuple[list[int], int, str]:
    text = (src / f"src/util/tables/sm_{cname}.c").read_text()
    body = re.search(r"static const TNCBIScore s_\w+PSM\[25 \* 25\] = \{(.*?)\};", text, re.S)
    tail = re.search(r"const SNCBIPackedScoreMatrix NCBISM_\w+ = \{\s*\"([^\"]+)\",\s*s_\w+PSM,\s*(-?\d+)", text)
    return numbers(body.group(1)), int(tail.group(2)), tail.group(1)


def losat_matrix(rs: str, name: str, variant: str) -> tuple[list[int], int]:
    if name == "BLOSUM62":
        body = re.search(r"pub const BLOSUM62: \[i8; BLOSUM62_SIZE \* BLOSUM62_SIZE\] = \[(.*?)\];", rs, re.S)
    else:
        body = re.search(rf"const {name}_TEXT: &str = r#\"(.*?)\"#;", rs, re.S)
    arm = re.search(rf"ScoringMatrix::{variant} => ProteinScoreMatrix \{{.*?defscore: (-?\d+),", rs, re.S)
    return numbers(body.group(1)), int(arm.group(1))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncbi-src", required=True, type=Path)
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[3])
    args = ap.parse_args()
    rs = (args.repo / "LOSAT/src/utils/matrix.rs").read_text()
    bad = 0
    for name, cname, variant in MATRICES:
        nscores, ndef, symbols = ncbi_matrix(args.ncbi_src, cname)
        lscores, ldef = losat_matrix(rs, name, variant)
        ok = symbols == "ARNDCQEGHILKMFPSTWYVBJZX*" and len(nscores) == 625 and nscores == lscores and ndef == ldef
        print(f"{name}\t{'same' if ok else 'DIFF'}\tscores={len(lscores)}\tdefscore={ldef}")
        bad += not ok
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
