#!/usr/bin/env python3
"""Save pinned NCBI Stage D traces for the remaining physical Stage C local FASTA cases."""
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
STAGE_D = HERE.parent
STAGE_C = HERE.parent.parent / "tlosan_stage_c"
sys.path.insert(0, str(STAGE_D))
from run_ncbi_parameter_trace import NCBI, NCBI_SHA256, NCBI_SOURCE, SOURCE_COMMIT, probe_run, sha

CASES = ("ambiguity_20260923", "lowercase_20260923", "multi_hsp_20260924",
         "multi_hsp_no_fence_20260924", "six_frame_merged_20260924")
D_PROBE = STAGE_D / "ncbi_d_call_trace.c"
K_PROBE = STAGE_D / "ncbi_kappa_traceback_trace.c"


def run_probe(command, plain, source, prefix, path):
    with tempfile.TemporaryDirectory(prefix="tlosan-d-remaining-") as tmp:
        run = probe_run(command, source, Path(tmp))
    selected = b''.join(line for line in run.stderr.splitlines(keepends=True) if line.startswith(prefix))
    ordinary = b''.join(line for line in run.stderr.splitlines(keepends=True) if not line.startswith(prefix))
    assert run.stdout == plain.stdout and ordinary == plain.stderr and selected
    path.write_bytes(selected)


def main():
    if len(sys.argv) != 2: raise SystemExit("usage: run_ncbi_remaining.py NEW_OUTPUT_DIRECTORY")
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(["git", "-C", str(NCBI_SOURCE), "rev-parse", "HEAD"], text=True).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256
    manifest = [f'NCBI source commit: {SOURCE_COMMIT}', f'NCBI executable SHA256: {NCBI_SHA256}',
                f'D probe SHA256: {sha(D_PROBE)}', f'Kappa probe SHA256: {sha(K_PROBE)}']
    for case in CASES:
        query = STAGE_C / case / 'query.faa'
        subject = STAGE_C / case / 'subjects.fna'
        manifest.extend([f'{case} query SHA256: {sha(query)}', f'{case} subject SHA256: {sha(subject)}'])
        for mode in (0, 2):
            name = f'{case}.mode{mode}'
            command = [str(NCBI), '-task', 'tblastn', '-query', str(query), '-subject', str(subject),
                       '-db_gencode', '1', '-num_threads', '1', '-evalue', '10000',
                       '-max_target_seqs', '500', '-comp_based_stats', str(mode),
                       '-sum_stats', 'false' if mode == 0 else 'true', '-seg', 'no',
                       '-outfmt', '6 qseqid sseqid score bitscore evalue nident positive length mismatch gaps gapopen qstart qend sstart send sframe']
            if case == 'lowercase_20260923': command.append('-lcase_masking')
            if case == 'multi_hsp_no_fence_20260924': command.extend(['-xdrop_gap', '5', '-xdrop_gap_final', '5'])
            plain = subprocess.run(command, capture_output=True, check=True)
            (output / f'{name}.report.tsv').write_bytes(plain.stdout)
            (output / f'{name}.stderr').write_bytes(plain.stderr)
            run_probe(command, plain, D_PROBE, b'D_', output / f'{name}.calls.trace')
            if mode == 2: run_probe(command, plain, K_PROBE, b'K_TRACE_', output / f'{name}.kappa.trace')
            manifest.append(f'{name} command: {command!r}')
    (output / 'manifest.txt').write_text('\n'.join(manifest) + '\n')
    (output / 'outputs.sha256').write_text(''.join(f'{sha(path)}  {path.name}\n' for path in sorted(output.iterdir()) if path.is_file()))

if __name__ == '__main__': main()
