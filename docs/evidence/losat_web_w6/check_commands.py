#!/usr/bin/env python3
"""Runs the reproduction commands of LOSAT Web under fixed conditions (S15 item 5, design §12.3).

The reproduction panel of a run (web/app/src/domain/reproduce.ts) gives, for each output format,
the LOSAT command and the NCBI BLAST+ command to compare it with, both from the run's argv. Design
§12.3 asks that the NCBI commands are not a mere change of the program's name but are checked
under fixed comparison conditions. This script takes the commands that reproduce.ts writes for a
fixed list of argvs (tests/unit/reproduce.test.ts, FIXED_CASES), puts small FASTA files of
LOSAT/tests/fasta/ (read only, named below) in one folder under the names of each argv, runs every
command there as the panel tells the user to (`bash -c '<command>'`), the LOSAT command with the
native CLI and the NCBI command with NCBI BLAST+ 2.17.0, and compares stdout, stderr and the exit
status byte for byte.

Status of each case and format:
  match               no exception applies, and both programs wrote the same bytes and status.
  differs             no exception applies, and the outputs differ: unexpected (exit 1).
  approved-exception  an approved exception applies (a non-default subject genetic code of
                      TBLASTN or TBLASTX, AGENTS.md); reported with "bytes equal" or "bytes
                      differ", never as a match.
  ncbi-refuses        reproduce.ts gives no NCBI command (NCBI cannot run the options as LOSAT
                      did); the NCBI command is run anyway and must fail, and LOSAT's must
                      succeed, or the status is unexpected (exit 1).

Usage (from the repository root; the machine-wide lock of the oracle runs, skill losat-oracle-runs):
  flock "$BUILD_ROOT/oracle.lock" python3 docs/evidence/losat_web_w6/check_commands.py \\
      --losat "$BUILD_ROOT/s14-gate-native/LOSAT" --ncbi-bin "$NCBI_BIN" \\
      --work "$BUILD_ROOT/s15-commands/run-<stamp>" --out <results.tsv> [--commands <json>] [--jobs 2]
Without --commands, the script writes them first with
  LOSAT_WEB_COMMANDS_OUT=<work>/commands.json npx vitest run tests/unit/reproduce.test.ts  (in web/app).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[3]
FASTA = REPOSITORY / 'LOSAT' / 'tests' / 'fasta'

# The FASTA files of each case of FIXED_CASES (query, subject), under LOSAT/tests/fasta/. They are
# fixtures of the certified outfmt 0/6/7 comparisons (LOSAT/tests/outfmt0_manifest.tsv: multi.blastn,
# e2d.blastn.dc, e2d.blastp.fast, width.tblastn, e2d.tblastn.both, e2d.tblastx.both) and, for
# BLASTX, the TBLASTN width fixture with the roles swapped.
FIXTURES = {
    'blastn.default': ('outfmt0/multi_query.fasta', 'outfmt0/multi_subject.fasta'),
    'blastn.options': ('outfmt0/multi_query.fasta', 'outfmt0/multi_subject.fasta'),
    'blastn.region': ('outfmt0/e2d_n_q.fa', 'outfmt0/e2d_n_s_plus.fa'),
    'blastp.default': ('outfmt0/e2d_p_q.faa', 'outfmt0/e2d_p_s3.faa'),
    'blastp.options': ('outfmt0/e2d_p_q.faa', 'outfmt0/e2d_p_s3.faa'),
    'tblastn.default': ('outfmt0/twidth_q.faa', 'outfmt0/twidth_s.fna'),
    'tblastn.options': ('outfmt0/e2d_t_pq.faa', 'outfmt0/e2d_t_ts.fa'),
    'tblastn.gencode4': ('outfmt0/twidth_q.faa', 'outfmt0/twidth_s.fna'),
    'tblastn.gencode32': ('outfmt0/twidth_q.faa', 'outfmt0/twidth_s.fna'),
    'tblastx.default': ('outfmt0/e2d_x_q.fa', 'outfmt0/e2d_t_ts.fa'),
    'tblastx.options': ('outfmt0/e2d_x_q.fa', 'outfmt0/e2d_t_ts.fa'),
    'tblastx.gencode5': ('outfmt0/e2d_x_q.fa', 'outfmt0/e2d_t_ts.fa'),
    'blastx.default': ('outfmt0/twidth_s.fna', 'outfmt0/twidth_q.faa'),
    'blastx.options': ('outfmt0/twidth_s.fna', 'outfmt0/twidth_q.faa'),
}
COLUMNS = [
    'case', 'program', 'outfmt', 'status', 'detail', 'losat_exit', 'ncbi_exit',
    'losat_stdout_sha256', 'ncbi_stdout_sha256', 'stdout_bytes', 'stderr_equal', 'losat_command', 'ncbi_command',
]
TIMEOUT_S = 600


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_commands(work: Path) -> Path:
    out = work / 'commands.json'
    env = {**os.environ, 'LOSAT_WEB_COMMANDS_OUT': str(out)}
    subprocess.run(['npx', 'vitest', 'run', 'tests/unit/reproduce.test.ts'], cwd=REPOSITORY / 'web' / 'app', env=env, check=True,
                   stdout=subprocess.DEVNULL)
    return out


def stage(case: dict, folder: Path) -> None:
    """Puts the case's FASTA files in its folder under the argv's names, as the panel says."""
    query, subject = FIXTURES[case['id']]
    folder.mkdir(parents=True, exist_ok=True)
    for name, source in ((case['query'], query), (case['subject'], subject)):
        target = folder / name
        if target.exists() and target.read_bytes() != (FASTA / source).read_bytes():
            raise SystemExit(f"{case['id']}: two different files would be named {name}")
        shutil.copyfile(FASTA / source, target)


def run(command: str, folder: Path, path: str, home: Path) -> tuple[int, bytes, bytes]:
    env = {'PATH': path, 'HOME': str(home), 'LC_ALL': 'C'}
    done = subprocess.run(['bash', '-c', command], cwd=folder, env=env, capture_output=True, timeout=TIMEOUT_S)
    return done.returncode, done.stdout, done.stderr


def check(case: dict, command: dict, folder: Path, path: str, home: Path) -> dict:
    losat = run(command['losat'], folder, path, home)
    ncbi = run(command['ncbi'], folder, path, home)
    same = losat == ncbi
    row = {
        'case': case['id'], 'program': case['program'], 'outfmt': command['format'],
        'losat_exit': losat[0], 'ncbi_exit': ncbi[0],
        'losat_stdout_sha256': sha256(losat[1]), 'ncbi_stdout_sha256': sha256(ncbi[1]),
        'stdout_bytes': len(losat[1]), 'stderr_equal': losat[2] == ncbi[2],
        'losat_command': command['losat'], 'ncbi_command': command['ncbi'],
    }
    if case['refused']:
        ok = losat[0] == 0 and ncbi[0] != 0
        first = ncbi[2].decode('utf-8', 'replace').strip().splitlines()
        row['status'] = 'ncbi-refuses' if ok else 'unexpected'
        row['detail'] = (f"NCBI exit {ncbi[0]}: {first[-1] if first else ''}" if ok
                         else f'expected LOSAT to search (exit 0) and NCBI to refuse; LOSAT exit {losat[0]}, NCBI exit {ncbi[0]}')
    elif case['exceptions']:
        row['status'] = 'approved-exception'
        row['detail'] = 'bytes equal' if same else 'bytes differ'
    else:
        row['status'] = 'match' if same else 'differs'
        row['detail'] = '' if same else ('stdout differs' if losat[1] != ncbi[1] else 'stderr or exit status differs')
    if losat[0] != 0 and not case['refused']:
        row['detail'] = (row['detail'] + f"; LOSAT stderr: {losat[2].decode('utf-8', 'replace').strip()[:200]}").strip('; ')
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--losat', required=True, type=Path, help='the native LOSAT CLI of the branch')
    parser.add_argument('--ncbi-bin', required=True, type=Path, help='the bin folder of NCBI BLAST+ 2.17.0')
    parser.add_argument('--work', required=True, type=Path, help='a new scratch folder')
    parser.add_argument('--out', required=True, type=Path, help='the TSV of results')
    parser.add_argument('--commands', type=Path, help='the JSON of tests/unit/reproduce.test.ts (written when absent)')
    parser.add_argument('--jobs', type=int, default=2, help='searches at once (each worker runs one at a time)')
    args = parser.parse_args()

    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    bin_dir = work / 'bin'
    bin_dir.mkdir()
    (bin_dir / 'LOSAT').symlink_to(args.losat.resolve())
    path = f'{bin_dir}:{args.ncbi_bin.resolve()}:/usr/bin:/bin'
    commands_path = args.commands.resolve() if args.commands else write_commands(work)
    commands = json.loads(commands_path.read_text())
    version = subprocess.run([str(args.ncbi_bin / 'blastn'), '-version'], capture_output=True, text=True).stdout.splitlines()[0]
    if '2.17.0' not in version:
        raise SystemExit(f'NCBI BLAST+ 2.17.0 expected, found {version!r}')

    rows: list[dict] = []
    with args.out.open('w') as out:
        out.write(f"# LOSAT {sha256(args.losat.read_bytes())} {args.losat}\n")
        out.write(f"# NCBI {version} {args.ncbi_bin}\n")
        out.write(f"# commands {sha256(commands_path.read_bytes())} {commands_path}\n")
        for case in commands['cases']:
            for name in FIXTURES[case['id']]:
                out.write(f"# fixture {case['id']} {sha256((FASTA / name).read_bytes())} LOSAT/tests/fasta/{name}\n")
        out.write('\t'.join(COLUMNS) + '\n')
        out.flush()
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            for case in commands['cases']:
                folder = work / case['id']
                stage(case, folder)
                for row in pool.map(lambda command: check(case, command, folder, path, work), case['commands']):
                    rows.append(row)
                    out.write('\t'.join(str(row[column]) for column in COLUMNS) + '\n')
                    out.flush()
                    print(f"{row['case']}\toutfmt {row['outfmt']}\t{row['status']}\t{row['detail']}", flush=True)
    unexpected = [row for row in rows if row['status'] in ('differs', 'unexpected')]
    counts = {status: sum(1 for row in rows if row['status'] == status) for status in sorted({row['status'] for row in rows})}
    print(f'{len(rows)} comparisons: ' + ', '.join(f'{status} {count}' for status, count in counts.items()))
    return 1 if unexpected else 0


if __name__ == '__main__':
    sys.exit(main())
