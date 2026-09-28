#!/usr/bin/env python3
"""Capture pinned local-subject NCBI score/rank reference, without LOSAT runtime."""
from __future__ import annotations

import hashlib
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
NCBI = Path('/home/kawato/micromamba/bin/tblastn')
NCBI_SHA = 'e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0'


def main() -> None:
    fixture = Path(sys.argv[1]).resolve()
    output = Path(sys.argv[2]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    assert hashlib.sha256(NCBI.read_bytes()).hexdigest() == NCBI_SHA
    common = [str(NCBI), '-task', 'tblastn', '-query', str(fixture / 'query.faa'),
              '-subject', str(fixture / 'subjects.fna'), '-db_gencode', '1',
              '-matrix', 'BLOSUM62', '-word_size', '3', '-threshold', '13',
              '-window_size', '40', '-gapopen', '11', '-gapextend', '1',
              '-evalue', '10000', '-num_threads', '1']
    manifest = [f'NCBI SHA256: {NCBI_SHA}',
                'NCBI source commit: 598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4',
                'Local -subject, one thread; NCBI reference only.']
    for name, comp, seg, sums in (
            ('raw', '0', 'no', 'false'),
            ('default', '2', '12 2.2 2.5', 'true')):
        command = common + ['-comp_based_stats', comp, '-seg', seg,
                            '-sum_stats', sums,
                            '-outfmt', '6 qseqid sseqid score bitscore evalue']
        run = subprocess.run(command, capture_output=True, check=True)
        if run.stderr:
            (output / f'{name}.stderr').write_bytes(run.stderr)
        rows = run.stdout.decode().splitlines()
        old = (fixture / f'{"raw_isolation" if name == "raw" else "default_profile"}_fields.out').read_text().splitlines()
        assert [(x.split('\t')[1], x.split('\t')[2]) for x in rows] == [
            (x.split('\t')[1], x.split('\t')[2]) for x in old]
        (output / f'{name}_score_bits_evalue.tsv').write_text(
            'rank\tquery\tsubject\traw_score\tbit_score\tevalue\n'
            + ''.join(f'{i}\t{row}\n' for i, row in enumerate(rows, 1)))
        manifest.append('Command: ' + repr(command))
    (output / 'manifest.txt').write_text('\n'.join(manifest) + '\n')
    files = sorted(path for path in output.iterdir() if path.name != 'outputs.sha256')
    (output / 'outputs.sha256').write_text(''.join(
        f'{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n' for path in files))
    print('NCBI raw and default score/bit/E-value/rank reference captured')


if __name__ == '__main__':
    main()
