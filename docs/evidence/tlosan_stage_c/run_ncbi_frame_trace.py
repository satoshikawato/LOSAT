#!/usr/bin/env python3
"""Record NCBI's preliminary translated subject bytes, comparison only."""
from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
NCBI = Path('/home/kawato/micromamba/bin/tblastn')
NCBI_SHA256 = 'e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0'
SOURCE_SHA = '598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4'
FRAMES = (1, 2, 3, -1, -2, -3)
FIELDS = '6 qseqid sseqid score qstart qend sstart send sframe qseq sseq'


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('fixture', type=Path)
    parser.add_argument('output', type=Path, help='new directory')
    args = parser.parse_args()
    fixture = args.fixture.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    output = args.output.resolve()
    assert hashlib.sha256(NCBI.read_bytes()).hexdigest() == NCBI_SHA256
    names = [line[1:].split()[0] for line in
             (fixture / 'subjects.fna').read_text().splitlines() if line.startswith('>')]
    with tempfile.TemporaryDirectory(prefix='tlosan-ncbi-frame-') as temp:
        probe = Path(temp) / 'ncbi_frame_trace.so'
        subprocess.run(['gcc', '-shared', '-fPIC', '-std=c11', '-O2', '-o',
                        str(probe), str(HERE / 'ncbi_frame_trace.c'), '-ldl'], check=True)
        command = [str(NCBI), '-task', 'tblastn', '-query', str(fixture / 'query.faa'),
                   '-subject', str(fixture / 'subjects.fna'), '-db_gencode', '1',
                   '-matrix', 'BLOSUM62', '-word_size', '3', '-threshold', '13',
                   '-window_size', '40', '-gapopen', '11', '-gapextend', '1',
                   '-evalue', '10000', '-num_threads', '1', '-comp_based_stats', '0',
                   '-seg', 'no', '-sum_stats', 'false', '-outfmt', FIELDS]
        env = os.environ.copy()
        env['LD_PRELOAD'] = str(probe)
        result = subprocess.run(command, env=env, capture_output=True, check=True)
        assert result.stdout == (fixture / 'raw_isolation_fields.out').read_bytes()
        rows = ['call\tsubject\tframe\tlength_aa\tsequence_ncbistdaa_hex\n']
        lines = result.stderr.decode().splitlines()
        assert len(lines) == 6 * len(names), len(lines)
        for call, line in enumerate(lines):
            tag, observed_call, observed_frame, length, aa_hex = line.split('\t')
            assert tag == 'FRAME' and int(observed_call) == call
            assert int(observed_frame) == FRAMES[call % 6]
            assert len(aa_hex) == 2 * int(length)
            rows.append(f'{call}\t{names[call // 6]}\t{observed_frame}\t{length}\t{aa_hex}\n')
        (output / 'ncbi_frame_trace.tsv').write_text(''.join(rows))
        (output / 'frame_trace_manifest.txt').write_text(
            f'NCBI source commit: {SOURCE_SHA}\nNCBI binary SHA256: {NCBI_SHA256}\n'
            f'Probe source SHA256: {hashlib.sha256((HERE / "ncbi_frame_trace.c").read_bytes()).hexdigest()}\n'
            'Comparison-only LD_PRELOAD; LOSAT runtime and build do not call NCBI.\n'
            'Final NCBI output bytes equal the unprobed raw_isolation_fields.out.\n'
            f'Command: {command!r}\n')
    files = ('ncbi_frame_trace.tsv', 'frame_trace_manifest.txt')
    (output / 'frame_trace_outputs.sha256').write_text(''.join(
        f'{hashlib.sha256((output / name).read_bytes()).hexdigest()}  {name}\n'
        for name in files))
    print(f'Captured {len(lines)} NCBI preliminary subject frames')


if __name__ == '__main__':
    main()
