#!/usr/bin/env python3
"""Capture every NCBI protein scan callback pair in its returned order."""
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
    parser.add_argument('--lcase-masking', action='store_true')
    args = parser.parse_args()
    fixture = args.fixture.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    output = args.output.resolve()
    assert hashlib.sha256(NCBI.read_bytes()).hexdigest() == NCBI_SHA256
    names = [line[1:].split()[0] for line in
             (fixture / 'subjects.fna').read_text().splitlines() if line.startswith('>')]
    with tempfile.TemporaryDirectory(prefix='tlosan-ncbi-candidates-') as temp:
        probe = Path(temp) / 'ncbi_candidate_trace.so'
        subprocess.run(['gcc', '-shared', '-fPIC', '-std=c11', '-O2', '-o',
                        str(probe), str(HERE / 'ncbi_candidate_trace.c'), '-ldl'], check=True)
        command = [str(NCBI), '-task', 'tblastn', '-query', str(fixture / 'query.faa'),
                   '-subject', str(fixture / 'subjects.fna'), '-db_gencode', '1',
                   '-matrix', 'BLOSUM62', '-word_size', '3', '-threshold', '13',
                   '-window_size', '40', '-gapopen', '11', '-gapextend', '1',
                   '-evalue', '10000', '-num_threads', '1', '-comp_based_stats', '0',
                   '-seg', 'no', '-sum_stats', 'false', '-outfmt', FIELDS]
        if args.lcase_masking:
            command.extend(['-lcase_masking'])
        env = os.environ.copy()
        env['LD_PRELOAD'] = str(probe)
        result = subprocess.run(command, env=env, capture_output=True, check=True)
        expected_name = 'lowercase_fields.out' if args.lcase_masking else 'raw_isolation_fields.out'
        assert result.stdout == (fixture / expected_name).read_bytes()
        rows = ['call\tsubject\tframe\tordinal\tquery_offset\tsubject_offset\n']
        ends = []
        parameters = []
        for line in result.stderr.decode().splitlines():
            parts = line.split('\t')
            call = int(parts[1])
            assert 0 <= call < len(names) * 6
            if parts[0] == 'CAND':
                assert len(parts) == 5
                rows.append(f'{call}\t{names[call // 6]}\t{FRAMES[call % 6]}\t'
                            + '\t'.join(parts[2:]) + '\n')
            elif parts[0] == 'PARAM':
                assert len(parts) == 5
                parameters.append((call, *(int(field) for field in parts[2:])))
            else:
                assert parts[0] == 'END' and len(parts) == 3
                ends.append((call, int(parts[2])))
        assert len(ends) == len(names) * 6, len(ends)
        assert len(parameters) == len(ends)
        assert all(call == i for i, (call, _) in enumerate(ends))
        assert sum(count for _, count in ends) == len(rows) - 1
        (output / 'ncbi_candidates.tsv').write_text(''.join(rows))
        (output / 'ncbi_word_params.tsv').write_text(
            'call\tx_dropoff_init\tx_dropoff\tcutoff_score\n'
            + ''.join('\t'.join(map(str, row)) + '\n' for row in parameters))
        (output / 'candidate_manifest.txt').write_text(
            f'NCBI source commit: {SOURCE_SHA}\nNCBI binary SHA256: {NCBI_SHA256}\n'
            f'Probe source SHA256: {hashlib.sha256((HERE / "ncbi_candidate_trace.c").read_bytes()).hexdigest()}\n'
            'Comparison-only LD_PRELOAD; LOSAT runtime and build do not call NCBI.\n'
            f'Final NCBI output bytes equal the unprobed {expected_name}.\n'
            f'WordFinder calls: {len(ends)}; candidate pairs: {len(rows)-1}\n'
            f'Command: {command!r}\n')
    files = ('ncbi_candidates.tsv', 'ncbi_word_params.tsv', 'candidate_manifest.txt')
    (output / 'candidate_outputs.sha256').write_text(''.join(
        f'{hashlib.sha256((output / name).read_bytes()).hexdigest()}  {name}\n'
        for name in files))
    print(f'Captured {len(rows)-1} candidate pairs from {len(ends)} NCBI WordFinder calls')


if __name__ == '__main__':
    main()
