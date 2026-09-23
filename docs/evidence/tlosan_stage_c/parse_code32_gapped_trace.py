#!/usr/bin/env python3
"""Extract code-32 gapped HSP trace from the Stage A comparison-only NCBI API."""
from __future__ import annotations

import hashlib
from pathlib import Path
import shutil
import sys

HERE = Path(__file__).resolve().parent
STAGE_A = HERE.parent / "tlosan_stage_a" / "api_20260923_verified"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: parse_code32_gapped_trace.py API_RUN_DIR NEW_OUTPUT_DIR")
    source = Path(sys.argv[1]).resolve()
    out = Path(sys.argv[2]).resolve()
    out.mkdir(parents=True, exist_ok=False)
    rows = [
        "code\tcall\tinit_index\tq_seed\ts_seed\tinit_q_start\tinit_s_start\t"
        "init_length\tinit_raw\thsp_index\thsp_raw\tcontext\tquery_frame\t"
        "q_start\tq_end\tq_gapped_start\tsubject_frame\ts_start\ts_end\ts_gapped_start\n"
    ]
    for code in (1, 32):
        output = source / f"code{code}_fmt6.out"
        expected = STAGE_A / f"code{code}_fmt6.out"
        assert output.read_bytes() == expected.read_bytes()
        stderr = source / f"code{code}_fmt6.stderr"
        shutil.copyfile(output, out / output.name)
        shutil.copyfile(stderr, out / stderr.name)
        lines = stderr.read_text().splitlines()
        inits = [line.split("\t") for line in lines if line.startswith("GAPPED_INIT\t")]
        hsps = [line.split("\t") for line in lines if line.startswith("GAPPED_HSP\t")]
        assert len(inits) == len(hsps)
        assert len(inits) == (1 if code == 1 else 2)
        for initial, hsp in zip(inits, hsps):
            assert initial[1] == hsp[1]
            assert initial[2] == hsp[2] == "0"
            rows.append("\t".join([str(code), *initial[1:], *hsp[2:]]) + "\n")
    (out / "gapped_code32_api.tsv").write_text("".join(rows))
    (out / "manifest.txt").write_text(
        "Comparison-only pinned NCBI C++ API -db oracle via FindGeneticCode(32).\n"
        "This establishes code-32 candidate-to-gapped HSP behavior, not local "
        "-subject statistics or headers.\n"
        "NCBI C/C++ source commit: 598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4\n"
        "NCBI CLI binary SHA256: e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0\n"
        f"Probe source SHA256: {sha(HERE / 'ncbi_gapped_trace.c')}\n"
        f"API run: {source}\n"
        "Rerun: gcc -shared -fPIC -std=c11 -O2 -o /tmp/tlosan-gapped-probe.so "
        "docs/evidence/tlosan_stage_c/ncbi_gapped_trace.c -ldl\n"
        "Rerun: LD_PRELOAD=/tmp/tlosan-gapped-probe.so bash "
        "docs/evidence/tlosan_stage_a/run_api_oracle.sh NEW_API_RUN_DIR\n"
        "Rerun: python3 docs/evidence/tlosan_stage_c/parse_code32_gapped_trace.py "
        "NEW_API_RUN_DIR NEW_OUTPUT_DIR\n"
    )
    paths = ["code1_fmt6.out", "code1_fmt6.stderr", "code32_fmt6.out",
             "code32_fmt6.stderr", "gapped_code32_api.tsv", "manifest.txt"]
    (out / "outputs.sha256").write_text(
        "".join(f"{sha(out / name)}  {name}\n" for name in paths))
    print("Code-1 control and code-32 C++ API output bytes unchanged; "
          "one and two gapped HSPs traced respectively.")


if __name__ == "__main__":
    main()
