#!/usr/bin/env python3
"""The examples of LOSAT Web ("Try an example") and the self-test's expected outputs.

The examples are small cuts of public records in LOSAT/tests/fasta (docs/web/examples.md):
one query and one subject file per program, in web/app/src/examples/. The self-test
(web/app/src/application/self-test.ts) runs them through the engine in the browser and
compares the SHA-256 of outfmt 0, 6 and 7 with web/app/src/examples/self-test.json, which
this script writes from the native LOSAT CLI of the engine commit that the site's reactors
were built from, and says for each output whether NCBI BLAST+ writes the same bytes.

  python3 web/app/tools/examples.py cut        # (re)write the example FASTA files
  python3 web/app/tools/examples.py check      # the files are what `cut` writes
  flock "$BUILD_ROOT/oracle.lock" python3 web/app/tools/examples.py expected \\
      --losat <native LOSAT of the engine commit> --engine-commit <commit> \\
      --ncbi-bin "$NCBI_BIN"                   # write self-test.json

Remake self-test.json whenever the engine changes: tests/unit/self-test-engine.test.ts
(with LOSAT_WEB_REACTORS) fails until it agrees with the reactors again. The file is part
of the release, so it names no machine path (plan TD-11).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

APP = Path(__file__).resolve().parents[1]
REPOSITORY = APP.parents[1]
SOURCES = REPOSITORY / 'LOSAT' / 'tests' / 'fasta'
EXAMPLES = APP / 'src' / 'examples'
TABLE = EXAMPLES / 'self-test.json'
FORMATS = (0, 6, 7)
ROLES = ('query', 'subject')
# The programs that LOSAT Web lists; BLASTX's example is shipped before BLASTX runs (SX).
PROGRAMS = ('blastn', 'blastp', 'blastx', 'tblastn', 'tblastx')


@dataclass(frozen=True)
class Region:
    """Residues start..end (1-based, both included) of the only record of a nucleotide file."""

    file: str
    start: int
    end: int


@dataclass(frozen=True)
class Protein:
    """One whole record of a protein file, by its accession, as the file has it."""

    file: str
    accession: str


# The first 2,000 nt of five Candidatus Hepatoplasma and Tyloplasma genomes, which hold the
# gene of DnaA (the chromosomal replication initiator), and their DnaA proteins; see
# docs/web/examples.md.
REGIONS = {
    'crinochetorum_tokyo': Region('AP027132.fasta', 1, 2000),
    'crinochetorum_av': Region('NZ_CP006932.fasta', 1, 2000),
    'scabrum': Region('AP027133.fasta', 1, 2000),
    'vulgare': Region('AP027131.fasta', 1, 2000),
    'tyloplasma': Region('AP027078.fasta', 1, 2000),
}
DNAA = {
    'vulgare': Protein('AP027131.faa', 'BDV02112.1'),
    'crinochetorum_tokyo': Protein('AP027132.faa', 'BDV02707.1'),
    'scabrum': Protein('AP027133.faa', 'BDV03280.1'),
    'tyloplasma': Protein('AP027078.faa', 'BDU67348.1'),
}
DNAA_SUBJECTS = [DNAA['crinochetorum_tokyo'], DNAA['scabrum'], DNAA['tyloplasma']]
REGION_SUBJECTS = [REGIONS['crinochetorum_tokyo'], REGIONS['scabrum'], REGIONS['tyloplasma']]
CUTS: dict[str, dict[str, list[Region | Protein]]] = {
    'blastn': {
        'query': [REGIONS['crinochetorum_tokyo']],
        'subject': [REGIONS['crinochetorum_av'], REGIONS['scabrum'], REGIONS['vulgare']],
    },
    'blastp': {'query': [DNAA['vulgare']], 'subject': DNAA_SUBJECTS},
    'blastx': {'query': [REGIONS['vulgare']], 'subject': DNAA_SUBJECTS},
    'tblastn': {'query': [DNAA['vulgare']], 'subject': REGION_SUBJECTS},
    'tblastx': {'query': [REGIONS['vulgare']], 'subject': REGION_SUBJECTS},
}


def file_name(program: str, role: str) -> str:
    """As domain/examples.ts `exampleFileName`: the name the search form gives the file."""
    return f'example-{program}-{role}.fasta'


def records(path: Path) -> list[tuple[str, list[str]]]:
    """The records of a FASTA file: the definition line and the sequence lines, as written."""
    out: list[tuple[str, list[str]]] = []
    for line in path.read_text(encoding='ascii').splitlines():
        if line.startswith('>'):
            out.append((line, []))
        elif out:
            out[-1][1].append(line)
        elif line.strip():
            raise SystemExit(f'{path}: sequence before the first definition line')
    return out


def cut(piece: Region | Protein) -> str:
    path = SOURCES / piece.file
    found = records(path)
    if isinstance(piece, Protein):
        for definition, lines in found:
            if definition[1:].split()[0] == piece.accession:
                return '\n'.join([definition, *lines]) + '\n'
        raise SystemExit(f'{path}: no record {piece.accession}')
    if len(found) != 1:
        raise SystemExit(f'{path}: expected one record, found {len(found)}')
    definition, lines = found[0]
    residues = ''.join(lines)
    if not 1 <= piece.start <= piece.end <= len(residues):
        raise SystemExit(f'{path}: region {piece.start}-{piece.end} outside 1-{len(residues)}')
    width = len(lines[0])
    part = residues[piece.start - 1 : piece.end]
    body = [part[i : i + width] for i in range(0, len(part), width)]
    # The definition line stays as in the record; the region is said after it.
    return '\n'.join([f'{definition} [region {piece.start}-{piece.end}]', *body]) + '\n'


def example_text(program: str, role: str) -> str:
    return ''.join(cut(piece) for piece in CUTS[program][role])


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def command_cut(_: argparse.Namespace) -> None:
    EXAMPLES.mkdir(parents=True, exist_ok=True)
    for program in PROGRAMS:
        for role in ROLES:
            data = example_text(program, role).encode('ascii')
            (EXAMPLES / file_name(program, role)).write_bytes(data)
            print(f'{file_name(program, role)}\t{len(data)}\t{sha256(data)}')


def command_check(_: argparse.Namespace) -> None:
    wrong = [
        file_name(program, role)
        for program in PROGRAMS
        for role in ROLES
        if (EXAMPLES / file_name(program, role)).read_bytes() != example_text(program, role).encode('ascii')
    ]
    if wrong:
        raise SystemExit(f'not what `cut` writes: {", ".join(wrong)}')
    print('the example files are what `cut` writes')


def run(argv: list[str], cwd: Path) -> bytes:
    done = subprocess.run(argv, cwd=cwd, capture_output=True, check=False)
    if done.returncode != 0:
        raise SystemExit(f'{" ".join(argv)} failed ({done.returncode}): {done.stderr.decode(errors="replace")[:500]}')
    return done.stdout


def first_line(argv: list[str]) -> str:
    return subprocess.run(argv, capture_output=True, check=True, text=True).stdout.splitlines()[0].strip()


def command_expected(args: argparse.Namespace) -> None:
    losat = Path(args.losat).resolve()
    ncbi = Path(args.ncbi_bin).resolve()
    cases = []
    with tempfile.TemporaryDirectory(prefix='losat-web-examples-') as work:
        cwd = Path(work)
        for program in PROGRAMS:
            inputs = {}
            for role in ROLES:
                data = (EXAMPLES / file_name(program, role)).read_bytes()
                (cwd / file_name(program, role)).write_bytes(data)
                inputs[role] = {'file': file_name(program, role), 'sha256': sha256(data)}
            # The argv of the search form for an example: the program, the two inputs and
            # nothing else (default settings; LOSAT Web adds -num_threads, the CLI -outfmt).
            argv = [program, '-query', file_name(program, 'query'), '-subject', file_name(program, 'subject')]
            outputs = {}
            for fmt in FORMATS:
                ours = run([str(losat), *argv, '-outfmt', str(fmt)], cwd)
                theirs = run([str(ncbi / program), *argv[1:], '-outfmt', str(fmt)], cwd)
                outputs[str(fmt)] = {
                    'sha256': sha256(ours),
                    'bytes': len(ours),
                    'ncbi': 'same' if ours == theirs else 'differs',
                }
                print(f'{program}\toutfmt {fmt}\t{len(ours)} bytes\t{sha256(ours)[:16]}\tNCBI {outputs[str(fmt)]["ncbi"]}')
            cases.append({'program': program, 'argv': argv, 'inputs': inputs, 'outputs': outputs})
    table = {
        'format': 1,
        'about': 'Expected outputs of the examples (docs/web/examples.md) for the browser self-test; '
        'remake with web/app/tools/examples.py expected when the engine changes.',
        'engine': {
            'commit': args.engine_commit,
            'version': first_line([str(losat), '--version']),
            'nativeSha256': sha256(losat.read_bytes()),
        },
        'oracle': 'NCBI BLAST+ ' + first_line([str(ncbi / 'blastn'), '-version']).removeprefix('blastn: '),
        'cases': cases,
    }
    TABLE.write_text(json.dumps(table, indent=2) + '\n', encoding='utf-8')
    print(f'wrote {TABLE.relative_to(REPOSITORY)}')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest='command', required=True)
    commands.add_parser('cut', help='write the example FASTA files').set_defaults(handler=command_cut)
    commands.add_parser('check', help='check the example FASTA files').set_defaults(handler=command_check)
    expected = commands.add_parser('expected', help='write self-test.json')
    expected.add_argument('--losat', required=True, help='native LOSAT CLI built from the engine commit')
    expected.add_argument('--engine-commit', required=True, help='the engine commit of that CLI and of the reactors')
    expected.add_argument('--ncbi-bin', required=True, help='bin directory of NCBI BLAST+ 2.17.0')
    expected.set_defaults(handler=command_expected)
    args = parser.parse_args()
    args.handler(args)


if __name__ == '__main__':
    sys.exit(main())
