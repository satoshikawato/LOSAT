#!/usr/bin/env python3
"""Build a fixed subject with every NCBI4na ambiguity mask."""
from __future__ import annotations
import argparse
import hashlib
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
NCBI = Path('/home/kawato/micromamba/bin/tblastn')
NCBI_SHA256 = 'e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0'


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('output', type=Path, help='new directory')
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    out = args.output.resolve()
    assert hashlib.sha256(NCBI.read_bytes()).hexdigest() == NCBI_SHA256
    original = ROOT / 'docs/evidence/tlosan_stage_a/fixtures/subject_code1.fna'
    core = ''.join(line.strip() for line in original.read_text().splitlines()
                   if line and not line.startswith('>'))
    assert len(core) == 360
    letters = 'MRSVWYHKDBN'
    spectrum = core[:90] + letters + core[90 + len(letters):]
    (out / 'query.faa').write_bytes((HERE / 'run_20260923/query.faa').read_bytes())
    (out / 'subjects.fna').write_text(f'>spectrum\n{spectrum}\n')
    cmd = [str(NCBI), '-task', 'tblastn', '-query', str(out / 'query.faa'),
           '-subject', str(out / 'subjects.fna'), '-db_gencode', '1',
           '-matrix', 'BLOSUM62', '-word_size', '3', '-threshold', '13',
           '-window_size', '40', '-gapopen', '11', '-gapextend', '1',
           '-evalue', '10000', '-num_threads', '1', '-comp_based_stats', '0',
           '-seg', 'no', '-sum_stats', 'false', '-outfmt',
           '6 qseqid sseqid score qstart qend sstart send sframe qseq sseq']
    result = subprocess.run(cmd, capture_output=True, check=True)
    (out / 'raw_isolation_fields.out').write_bytes(result.stdout)
    (out / 'raw_isolation_fields.stderr').write_bytes(result.stderr)
    (out / 'oracle_manifest.txt').write_text(
        f'NCBI binary SHA256: {NCBI_SHA256}\n'
        'NCBI source commit: 598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4\n'
        f'Changed subject positions: {list(range(90, 90 + len(letters)))}\n'
        f'IUPAC ambiguity letters: {letters}\nCommand: {cmd!r}\n')
    files = ('query.faa', 'subjects.fna', 'raw_isolation_fields.out',
             'raw_isolation_fields.stderr', 'oracle_manifest.txt')
    (out / 'oracle_outputs.sha256').write_text(''.join(
        f'{hashlib.sha256((out / name).read_bytes()).hexdigest()}  {name}\n'
        for name in files))
    print('NCBI ambiguity spectrum rows:', len(result.stdout.splitlines()))


if __name__ == '__main__':
    main()
