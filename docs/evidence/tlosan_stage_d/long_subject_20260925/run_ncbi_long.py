#!/usr/bin/env python3
"""Save the fixed long Stage C subject's pinned NCBI local Stage D call trace."""
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
STAGE_D = HERE.parent
STAGE_C = HERE.parent.parent / "tlosan_stage_c"
sys.path.insert(0, str(STAGE_D))
from run_ncbi_parameter_trace import NCBI, NCBI_SHA256, NCBI_SOURCE, SOURCE_COMMIT, probe_run, sha


def main():
    if len(sys.argv) != 2: raise SystemExit("usage: run_ncbi_long.py NEW_OUTPUT_DIRECTORY")
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    probe = STAGE_D / "ncbi_d_call_trace.c"
    assert subprocess.check_output(["git", "-C", str(NCBI_SOURCE), "rev-parse", "HEAD"], text=True).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256
    query = STAGE_C / "long_multi_query_20260924/query.faa"
    source = STAGE_C / "run_20260923/subjects.fna"
    records = {}
    for part in source.read_text().split('>'):
        if part.strip():
            lines = part.splitlines()
            records[lines[0].split()[0]] = ''.join(lines[1:])
    plus1 = records['plus1'].encode()
    subject_bytes = b'>chunk_edge\n' + b'ATG'*4_999_950 + plus1 + b'ATG'*250 + b'\n'
    assert len(subject_bytes) == 15_000_962 + len(b'>chunk_edge\n\n')
    with tempfile.TemporaryDirectory(prefix="tlosan-d-long-") as temp:
        subject = Path(temp) / 'subjects.fna'
        subject.write_bytes(subject_bytes)
        command = [str(NCBI), '-task', 'tblastn', '-query', str(query), '-subject', str(subject),
                   '-db_gencode', '1', '-matrix', 'BLOSUM62', '-word_size', '3', '-threshold', '13',
                   '-window_size', '40', '-gapopen', '11', '-gapextend', '1', '-evalue', '10000',
                   '-num_threads', '1', '-comp_based_stats', '0', '-seg', 'no', '-sum_stats', 'false',
                   '-outfmt', '6 qseqid sseqid score bitscore evalue nident positive length mismatch gaps gapopen qstart qend sstart send sframe']
        plain = subprocess.run(command, capture_output=True, check=True)
        (output / 'report.tsv').write_bytes(plain.stdout)
        (output / 'stderr').write_bytes(plain.stderr)
        run = probe_run(command, probe, Path(temp))
        selected = b''.join(line for line in run.stderr.splitlines(keepends=True) if line.startswith(b'D_'))
        ordinary = b''.join(line for line in run.stderr.splitlines(keepends=True) if not line.startswith(b'D_'))
        assert run.stdout == plain.stdout and ordinary == plain.stderr and selected
        (output / 'calls.trace').write_bytes(selected)
        manifest = [f'NCBI source commit: {SOURCE_COMMIT}', f'NCBI executable SHA256: {NCBI_SHA256}',
                    f'D probe SHA256: {sha(probe)}', f'query SHA256: {sha(query)}',
                    f'insert subject source SHA256: {sha(source)}', f'generated subject SHA256: {sha(subject)}',
                    'subject: ATG x 4999950 codons + saved plus1 + ATG x 250 codons',
                    'command: ' + repr([arg if arg != str(subject) else '<generated-subject>' for arg in command])]
    (output / 'manifest.txt').write_text('\n'.join(manifest) + '\n')
    (output / 'outputs.sha256').write_text(''.join(f'{sha(path)}  {path.name}\n' for path in sorted(output.iterdir()) if path.is_file()))

if __name__ == '__main__': main()
