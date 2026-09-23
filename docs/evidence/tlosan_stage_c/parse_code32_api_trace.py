#!/usr/bin/env python3
"""Extract comparison-only code-32 DB WordFinder traces from pinned API runs."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
STAGE_A = HERE.parent / 'tlosan_stage_a' / 'api_20260923_verified'
FRAMES = (1, 2, 3, -1, -2, -3)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('candidate_run', type=Path)
    parser.add_argument('init_run', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    for run in (args.candidate_run, args.init_run):
        assert (run / 'manifest.txt').read_text().startswith(
            'NCBI source commit: 598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4')
        assert (run / 'code32_fmt6.out').read_bytes() == (STAGE_A / 'code32_fmt6.out').read_bytes()
        assert (run / 'subject32_g1_fmt6.out').read_bytes() == (STAGE_A / 'subject32_g1_fmt6.out').read_bytes()
    lines = (args.candidate_run / 'code32_fmt6.stderr').read_text().splitlines()
    candidates = ['frame\tordinal\tquery_offset\tsubject_offset\n']
    params = ['frame\tx_dropoff_init\tx_dropoff\tcutoff_score\n']
    ends = []
    for line in lines:
        fields = line.split('\t')
        if fields[0] == 'CAND':
            assert len(fields) == 5
            call = int(fields[1])
            candidates.append(f'{FRAMES[call]}\t' + '\t'.join(fields[2:]) + '\n')
        elif fields[0] == 'PARAM':
            assert len(fields) == 5
            params.append(f'{FRAMES[int(fields[1])]}\t' + '\t'.join(fields[2:]) + '\n')
        else:
            assert fields[0] == 'END' and len(fields) == 3
            ends.append((int(fields[1]), int(fields[2])))
    assert ends == [(i, count) for i, count in enumerate((134, 29, 8, 4, 2, 22))]
    assert sum(count for _, count in ends) == len(candidates) - 1 == 199
    (output / 'ncbi_candidates.tsv').write_text(''.join(candidates))
    (output / 'ncbi_word_params.tsv').write_text(''.join(params))

    lines = (args.init_run / 'code32_fmt6.stderr').read_text().splitlines()
    hsps = ['frame\tinit_index\tq_seed\ts_seed\tq_start\ts_start\tlength\traw_score\n']
    calls = []
    for line in lines:
        fields = line.split('\t')
        if fields[0] == 'WORD_FINDER':
            assert len(fields) == 3
            calls.append((int(fields[1]), int(fields[2])))
        elif fields[0] == 'INIT':
            assert len(fields) == 9
            hsps.append(f'{FRAMES[int(fields[1])]}\t' + '\t'.join(fields[2:]) + '\n')
    assert len(calls) == 6 and [call for call, _ in calls] == list(range(6))
    assert sum(count for _, count in calls) == len(hsps) - 1
    (output / 'ncbi_init_trace.tsv').write_text(''.join(hsps))
    (output / 'manifest.txt').write_text(
        'Pinned NCBI source: 598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4\n'
        'Stage A C++ API comparison-only DB search, code 32, one thread.\n'
        'Database statistics/headers are not local -subject oracles.\n'
        'Probes: ncbi_candidate_trace.c and ncbi_wordfinder_trace.c via LD_PRELOAD; '
        'never LOSAT runtime or build.\n'
        f'Candidate probe SHA256: {digest(HERE / "ncbi_candidate_trace.c")}\n'
        f'WordFinder probe SHA256: {digest(HERE / "ncbi_wordfinder_trace.c")}\n'
        f'code32_fmt6 output SHA256: {digest(args.candidate_run / "code32_fmt6.out")} '
        '(both probes byte-identical to Stage A unprobed output).\n'
        f'Six WordFinder calls, {len(candidates)-1} ordered candidate pairs, '
        f'{len(hsps)-1} saved initial HSPs.\n'
    )
    files = ('ncbi_candidates.tsv', 'ncbi_word_params.tsv', 'ncbi_init_trace.tsv', 'manifest.txt')
    (output / 'outputs.sha256').write_text(''.join(
        f'{digest(output / name)}  {name}\n' for name in files))
    print(f'code32: {len(candidates)-1} ordered candidates, {len(hsps)-1} init HSPs')


if __name__ == '__main__':
    main()
