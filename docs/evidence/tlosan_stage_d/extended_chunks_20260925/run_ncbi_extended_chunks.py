#!/usr/bin/env python3
"""Replay four generated Stage C local subjects through pinned NCBI Stage D."""
import hashlib
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
STAGE_D = HERE.parent
STAGE_C = HERE.parent.parent / 'tlosan_stage_c'
NCBI = Path('/home/kawato/micromamba/bin/tblastn')
NCBI_SHA = 'e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0'
SOURCE = Path('/mnt/c/Users/genom/GitHub/ncbi-blast')
COMMIT = '598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4'
CASES = ('long_chunk_20260924', 'masked_chunk_boundary_20260924',
         'no_range_middle_20260924', 'no_range_two_hits_20260924')
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()

def make_subject(case, insert):
    if case.startswith('no_range'):
        if case == 'no_range_two_hits_20260924':
            first, second = 4_998_450, 10_000_050
            subject = bytearray(b'ATG' * first + insert + b'A' +
                                b'ATG' * (second - first - 121) + insert + b'ATG' * 250)
        else:
            subject = bytearray(b'ATG' * 10_000_050 + insert + b'ATG' * 250)
        subject[4_999_000 * 3:10_000_000 * 3] = subject[4_999_000 * 3:10_000_000 * 3].lower()
        return b'>masked_middle\n' + subject + b'\n'
    subject = bytearray(b'ATG' * 4_999_950 + insert + b'ATG' * 250)
    if case == 'masked_chunk_boundary_20260924':
        subject[4_999_900 * 3:5_000_000 * 3] = subject[4_999_900 * 3:5_000_000 * 3].lower()
    return b'>chunk_edge\n' + subject + b'\n'

def main():
    out = Path(sys.argv[1]).resolve()
    out.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(['git', '-C', str(SOURCE), 'rev-parse', 'HEAD'], text=True).strip() == COMMIT
    assert sha(NCBI) == NCBI_SHA
    d_source = STAGE_D / 'ncbi_d_call_trace.c'
    k_source = STAGE_D / 'ncbi_kappa_traceback_trace.c'
    baseline = (STAGE_C / 'run_20260923/subjects.fna').read_text()
    records = {block.splitlines()[0]: ''.join(block.splitlines()[1:]) for block in baseline.split('>') if block.strip()}
    insert = records['plus1'].encode()
    assert len(insert) == 362
    manifest = [f'NCBI source commit: {COMMIT}', f'NCBI CLI SHA256: {NCBI_SHA}',
                f'D probe SHA256: {sha(d_source)}', f'Kappa probe SHA256: {sha(k_source)}',
                'Subjects reconstructed exactly from pinned Stage C scripts; query/subject SHA256 checked against each Stage C manifest.',
                'No NCBI artifact is used by LOSAT execution or build.']
    with tempfile.TemporaryDirectory(prefix='tlosan-d-extended-') as td:
        tmp = Path(td)
        probes = {}
        for kind, source in (('D', d_source), ('K', k_source)):
            lib = tmp / (kind + '.so')
            subprocess.run(['gcc', '-std=c11', '-O0', '-shared', '-fPIC', str(source), '-ldl', '-o', str(lib)], check=True)
            probes[kind] = lib
        for case in CASES:
            query = tmp / 'query.faa'
            subject = tmp / 'subjects.fna'
            query.write_bytes((STAGE_C / case / 'query.faa').read_bytes())
            subject.write_bytes(make_subject(case, insert))
            prior = (STAGE_C / case / 'manifest.txt').read_text()
            match = re.search(r'Input SHA256: query=([0-9a-f]+) subject=([0-9a-f]+)', prior)
            assert match and (sha(query), sha(subject)) == match.groups(), case
            manifest.append(f'{case}: query SHA256 {sha(query)}; subject SHA256 {sha(subject)}')
            for mode in (0, 2):
                stem = f'{case}.mode{mode}'
                cmd = [str(NCBI), '-task', 'tblastn', '-query', str(query), '-subject', str(subject),
                       '-db_gencode', '1', '-num_threads', '1', '-evalue', '10000',
                       '-max_target_seqs', '500', '-comp_based_stats', str(mode),
                       '-sum_stats', 'false' if mode == 0 else 'true', '-seg', 'no',
                       '-outfmt', '6 qseqid sseqid score bitscore evalue nident positive length mismatch gaps gapopen qstart qend sstart send sframe']
                if case != 'long_chunk_20260924': cmd.append('-lcase_masking')
                plain = subprocess.run(cmd, capture_output=True, check=True)
                (out / f'{stem}.report.tsv').write_bytes(plain.stdout)
                (out / f'{stem}.stderr').write_bytes(plain.stderr)
                for kind, prefix, suffix in (('D', b'D_', 'calls.trace'), ('K', b'K_TRACE_', 'kappa.trace')):
                    if kind == 'K' and mode == 0: continue
                    env = os.environ.copy(); env['LD_PRELOAD'] = str(probes[kind])
                    traced = subprocess.run(cmd, capture_output=True, check=True, env=env)
                    selected = b''.join(line for line in traced.stderr.splitlines(keepends=True) if line.startswith(prefix))
                    ordinary = b''.join(line for line in traced.stderr.splitlines(keepends=True) if not line.startswith(prefix))
                    assert selected and traced.stdout == plain.stdout and ordinary == plain.stderr, (case, mode, kind)
                    (out / f'{stem}.{suffix}').write_bytes(selected)
                manifest.append(f'{stem}: {len(plain.stdout.splitlines())} report rows; command {repr(cmd).replace(str(tmp), "$TMP")}')
    (out / 'manifest.txt').write_text('\n'.join(manifest) + '\n')
    (out / 'outputs.sha256').write_text(''.join(f'{sha(path)}  {path.name}\n' for path in sorted(out.iterdir()) if path.is_file()))

if __name__ == '__main__': main()
