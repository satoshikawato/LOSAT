#!/usr/bin/env python3
"""Run CLI-calibrated pinned NCBI local API for all 27 subject codes."""
import hashlib
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
STAGE_D = HERE.parent
STAGE_A = HERE.parent.parent / "tlosan_stage_a"
OUT, ORACLE, KAPPA, D_TRACE, CLI = map(Path, sys.argv[1:])
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
(OUT / "oracle_binary.sha256").write_text(f"{sha(ORACLE)}  tblastn_local_api_oracle\n")
manifest = [
    "NCBI source commit: 598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4",
    f"NCBI CLI SHA256: {sha(CLI)}",
    f"NCBI API oracle source SHA256: {sha(STAGE_D / 'tblastn_code32_local_oracle.cpp')}",
    f"NCBI Kappa probe SHA256: {sha(STAGE_D / 'ncbi_kappa_traceback_trace.c')}",
    f"NCBI D probe SHA256: {sha(STAGE_D / 'ncbi_d_call_trace.c')}",
]

def run(cmd, probe=None, prefix=None):
    env = os.environ.copy()
    if probe: env["LD_PRELOAD"] = str(probe)
    result = subprocess.run(cmd, capture_output=True, check=True, env=env)
    if prefix:
        selected = b"".join(line for line in result.stderr.splitlines(keepends=True)
                            if line.startswith(prefix))
        ordinary = b"".join(line for line in result.stderr.splitlines(keepends=True)
                            if not line.startswith(prefix))
        assert selected
        return result.stdout, ordinary, selected
    return result.stdout, result.stderr

# NCBI core/blast_engine.c:1407,1434-1443; blast_kappa.c:3577-3736:
# the comparison-only API must first reproduce the pinned code-1 CLI local
# stdout and exact Stage D call state on the retained code-32 subject.
cal_query = STAGE_A / "fixtures/query.faa"
cal_subject = STAGE_A / "fixtures/subject_code32.fna"
api = [str(ORACLE), str(cal_query), str(cal_subject), "1", "6"]
cli = [str(CLI), "-task", "tblastn", "-query", str(cal_query), "-subject",
       str(cal_subject), "-db_gencode", "1", "-num_threads", "1", "-outfmt", "6"]
api_out, api_err = run(api)
cli_out, cli_err = run(cli)
assert api_out == cli_out
assert api_out == (STAGE_D / "code32_local_api_20260925_cli_calibrated/code1_fmt6.out").read_bytes()
_, ordinary, d_state = run(api, D_TRACE, b"D_")
assert ordinary == api_err
assert d_state == (STAGE_D / "code32_local_api_20260925_cli_calibrated/code1_d.trace").read_bytes()
(OUT / "calibration.report.tsv").write_bytes(api_out)
(OUT / "calibration.d.trace").write_bytes(d_state)
manifest.append(f"Calibration API command: {(["$NCBI_LOCAL_API_ORACLE"] + api[1:])!r}")
manifest.append(f"Calibration CLI command: {cli!r}")
for row in (HERE / "fixtures.tsv").read_text().splitlines()[1:]:
    code, codon, standard, selected, positions, note = row.split("\t")
    q = HERE / "fixtures" / f"code{code}.faa"
    s = HERE / "fixtures" / f"code{code}.fna"
    for selected_code, label in ((1, "control"), (int(code), "selected")):
        cmd = [str(ORACLE), str(q), str(s), str(selected_code), "6"]
        plain, stderr = run(cmd)
        k_out, k_err, k_trace = run(cmd, KAPPA, b"K_TRACE_")
        d_out, d_err, d_trace = run(cmd, D_TRACE, b"D_")
        assert (k_out, d_out, k_err, d_err) == (plain, plain, stderr, stderr)
        stem = f"code{code}.{label}"
        (OUT / f"{stem}.report.tsv").write_bytes(plain)
        (OUT / f"{stem}.stderr").write_bytes(stderr)
        (OUT / f"{stem}.kappa.trace").write_bytes(k_trace)
        (OUT / f"{stem}.d.trace").write_bytes(d_trace)
        manifest.append(f"{stem} command: {(["$NCBI_LOCAL_API_ORACLE"] + cmd[1:])!r}")
    if note == "identical internal translation":
        assert (OUT / f"code{code}.control.report.tsv").read_bytes() == (OUT / f"code{code}.selected.report.tsv").read_bytes()
    else:
        assert (OUT / f"code{code}.control.report.tsv").read_bytes() != (OUT / f"code{code}.selected.report.tsv").read_bytes(), code
(OUT / "manifest.txt").write_text("\n".join(manifest) + "\n")
files = sorted(path for path in OUT.iterdir() if path.is_file() and path.suffix in {".tsv", ".trace", ".stderr", ".txt"})
(OUT / "outputs.sha256").write_text("".join(f"{sha(path)}  {path.name}\n" for path in files))
