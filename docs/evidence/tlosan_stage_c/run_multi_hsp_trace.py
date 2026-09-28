#!/usr/bin/env python3
"""Reproduce a long local-subject TBLASTN comparison-only Stage C trace."""
from __future__ import annotations

import hashlib
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
NCBI_SOURCE = Path('/mnt/c/Users/genom/GitHub/ncbi-blast')
SOURCE_COMMIT = '598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4'
NCBI = Path('/home/kawato/micromamba/bin/tblastn')
NCBI_SHA256 = 'e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0'


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(['git', '-C', str(NCBI_SOURCE), 'rev-parse', 'HEAD'], text=True).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256
    original = HERE / 'run_20260923'
    query = (original / 'query.faa').read_bytes()
    records = {}
    for block in (original / 'subjects.fna').read_text().split('>'):
        if block.strip():
            name, *lines = block.splitlines()
            records[name] = ''.join(lines)
    randomizer = random.Random(20260924)
    def filler(length: int) -> str:
        return ''.join(randomizer.choice('ACGT') for _ in range(length))
    # Pinned NCBI c++/src/algo/blast/core/blast_engine.c:804-850:
    # for (context=first_context; context<=last_context; context++) {
    #     s_BlastSearchEngineOneContext(...);
    #     Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax);
    # }
    mode = sys.argv[2] if len(sys.argv) == 3 else ''
    assert mode in ('', '--six-frame', '--no-fence'), 'usage: run_multi_hsp_trace.py NEW_OUTPUT_DIR [--six-frame|--no-fence]'
    if mode == '--six-frame':
        names = ('plus1', 'plus2', 'plus3', 'minus1', 'minus2', 'minus3')
        subject = filler(600) + ''.join(records[name] + filler(601) for name in names[:-1]) + records[names[-1]] + filler(600)
        subject_id = 'six_frames'
    else:
        subject = filler(440) + records['plus1'] + filler(700) + records['plus1'] + filler(500)
        subject_id = 'two_copies'
    (output / 'query.faa').write_bytes(query)
    (output / 'subjects.fna').write_text('>' + subject_id + '\n' + subject + '\n')
    # Pinned NCBI c++/src/algo/blast/core/blast_traceback.c:1644-1684:
    # Blast_TracebackFromHSPList(...); if (fence_hit) retry on full subject.
    command = [str(NCBI), '-task', 'tblastn', '-query', str(output / 'query.faa'),
               '-subject', str(output / 'subjects.fna'), '-db_gencode', '1',
               '-matrix', 'BLOSUM62', '-word_size', '3', '-threshold', '13',
               '-window_size', '40', '-gapopen', '11', '-gapextend', '1',
               '-evalue', '10000', '-num_threads', '1', '-comp_based_stats', '0',
               '-seg', 'no', '-sum_stats', 'false', '-outfmt',
               '6 qseqid sseqid score qstart qend sstart send sframe qseq sseq']
    # NCBI c++/src/algo/blast/core/blast_parameters.c:455-463:
    # gap_x_dropoff = (Int4)(options->gap_x_dropoff*NCBIMATH_LN2/min_lambda);
    # gap_x_dropoff_final = (Int4)MAX(options->gap_x_dropoff_final*NCBIMATH_LN2/min_lambda, params->gap_x_dropoff);
    # The reduced pair yields raw 12 for preliminary and final extension,
    # exposing successful partial-window traceback without a fence.
    if mode == '--no-fence':
        command.extend(['-xdrop_gap', '5', '-xdrop_gap_final', '5'])
    reference = subprocess.run(command, capture_output=True, check=True)
    with tempfile.TemporaryDirectory(prefix='tlosan-c-multihsp-') as temp:
        probe = Path(temp) / 'ncbi_gapped_trace.so'
        subprocess.run(['gcc', '-shared', '-fPIC', '-std=c11', '-O2', '-Wall', '-Wextra',
                        '-o', str(probe), str(HERE / 'ncbi_gapped_trace.c'), '-ldl'], check=True)
        env = os.environ.copy()
        env['LD_PRELOAD'] = str(probe)
        traced = subprocess.run(command, capture_output=True, check=True, env=env)
    assert traced.stdout == reference.stdout, 'comparison probe changed NCBI output bytes'
    (output / 'ncbi_output.out').write_bytes(reference.stdout)
    (output / 'ncbi_gapped_trace.stderr').write_bytes(traced.stderr)
    lines = traced.stderr.decode().splitlines()
    selected = {
        'GAPPED_HSP': ('gapped_hsps.tsv',
                       'event\tcall\tindex\traw_score\tcontext\tquery_frame\tq_start\tq_end\tq_gapped_start\tsubject_frame\ts_start\ts_end\ts_gapped_start'),
        'TARGET_TRANSLATION': ('translation_windows.tsv',
                               'event\tframe\thsp_s_start\thsp_s_end\twindow_start\twindow_stop\treturned_length\tpartial'),
        'TRACEBACK_INPUT': ('traceback_events.tsv', 'event\ta\tb\tc\td\te\tf\tg\th'),
        'ENDPOINT_INPUT': ('endpoint_events.tsv', 'event\ta\tb\tc\td\te\tf\tg'),
    }
    for event, (filename, header) in selected.items():
        if event == 'TRACEBACK_INPUT':
            rows = [line for line in lines if line.startswith('TRACEBACK_')]
        elif event == 'ENDPOINT_INPUT':
            rows = [line for line in lines if line.startswith('ENDPOINT_')]
        else:
            rows = [line for line in lines if line.startswith(event + '\t')]
        (output / filename).write_text(header + '\n' + '\n'.join(rows) + '\n')
    counts = {event: sum(line.startswith(event + '\t') for line in lines)
              for event in ('GAPPED_HSP', 'TARGET_TRANSLATION', 'TRACEBACK_INPUT',
                            'TRACEBACK_OUTPUT', 'REEVAL_INPUT', 'REEVAL_OUTPUT',
                            'ENDPOINT_INPUT', 'ENDPOINT_OUTPUT')}
    assert counts['GAPPED_HSP'] >= 2
    expected_passes = 1 if mode == '--no-fence' else 2
    assert counts['TRACEBACK_INPUT'] == expected_passes
    assert counts['TRACEBACK_OUTPUT'] == expected_passes
    if mode == '--no-fence':
        assert any(line.startswith('TRACEBACK_OUTPUT\t0\t0\t') for line in lines)
    (output / 'manifest.txt').write_text(
        'Comparison-only NCBI trace. Probe never linked or called by LOSAT.\n'
        f'NCBI source commit: {SOURCE_COMMIT}\n'
        f'NCBI binary SHA256: {NCBI_SHA256}\n'
        f'Probe source SHA256: {sha(HERE / "ncbi_gapped_trace.c")}\n'
        f'Subject length: {len(subject)} nt\n'
        f'Counts: {counts}\n'
        'Final NCBI bytes unchanged by probe: yes\n'
        f'Command: {command!r}\n')
    files = ('query.faa', 'subjects.fna', 'ncbi_output.out',
             'ncbi_gapped_trace.stderr', 'gapped_hsps.tsv',
             'translation_windows.tsv', 'traceback_events.tsv',
             'endpoint_events.tsv', 'manifest.txt')
    (output / 'outputs.sha256').write_text(''.join(
        f'{sha(output / name)}  {name}\n' for name in files))
    print(counts)


if __name__ == '__main__':
    main()
