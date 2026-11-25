# cleanup-blastdb-volumes.py

`losat cleanup-blastdb-volumes.py` removes orphaned BLAST database volume files that are no longer listed in a database's `.pal` (protein) or `.nal` (nucleotide) alias file. It is a pure-Rust implementation of NCBI BLAST+ 2.17.0's `cleanup-blastdb-volumes.py` utility and needs no Python runtime.

## Synopsis

```bash
losat cleanup-blastdb-volumes.py -db DB -dbtype {prot,nucl} [-dry-run] [-verbose]
losat cleanup-blastdb-volumes.py -h
losat cleanup-blastdb-volumes.py -version
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `cleanup-blastdb-volumes.py ...`.

## Options

| Option | Meaning |
| --- | --- |
| `-h`, `--help` | Print the help text and exit with status `0`. |
| `-db DB` | Required. BLAST database name (searched directly, then in `BLASTDB`, then via `[BLAST] BLASTDB` in `.ncbirc` or `ncbi.ini` under the working directory, home directory, and `$NCBI`). |
| `-dbtype {prot,nucl}` | Required. Molecule type (`prot` for `.pal`/`.pin`/`.p??`, `nucl` for `.nal`/`.nin`/`.n??`). |
| `-dry-run` | Do not delete any files; print each extra volume (`Will remove extra volume ...`). |
| `-verbose` | With `-dry-run`, also print each file that would be removed (`Will remove ...`). |
| `-version` | Print the version line (`cleanup-blastdb-volumes.py LOSAT 0.2.0`) and exit with status `0`. |

## What it removes

1. The database is located by checking `{db}.pal`/`.pin` (or `.nal`/`.nin`) directly, under the `BLASTDB` environment variable, and under `[BLAST] BLASTDB` from `.ncbirc` or `ncbi.ini` in the current working directory, the user's home directory, and `$NCBI`. If no alias or index file is found, it writes `Cannot find protein|nucleotide <db> BLAST database` to standard error and exits with status `1`.
2. If the database has no `.pal` (protein) or `.nal` (nucleotide) alias file (for example, a single-volume database), it exits silently with status `1`, matching NCBI's script.
3. For each `DBLIST` line in the alias file, every volume matching `{db}.*.pin` (or `.nin`) whose base name is not listed on that line is treated as an extra volume.
4. For each extra volume, it removes every matching `{vol}.p??` (or `.n??`) file and `{vol}.tar.gz.md5`, printing `Removed <path>` for each removed file (or, with `-dry-run`, printing `Will remove extra volume <vol>` and, when `-verbose` is also given, `Will remove <path>`).

## Differences from NCBI BLAST+

- **`-version`**: prints `cleanup-blastdb-volumes.py LOSAT 0.2.0` where NCBI's script prints `cleanup-blastdb-volumes.py 1.0`.
- **Python version**: LOSAT matches NCBI BLAST+ 2.17.0's `cleanup-blastdb-volumes.py` running under Python 3.13 (including `argparse` wording and terminal-width wrapping). Under Python 3.12, `argparse` quotes invalid `-dbtype` choices (`(choose from 'prot', 'nucl')`), whereas Python 3.13 and LOSAT print `(choose from prot, nucl)`.
- **Not a C++ Toolkit application**: like NCBI's `cleanup-blastdb-volumes.py`, this subcommand uses Python `argparse`-style options (exit status `2` on argument errors) and does not take NCBI C++ Toolkit flags such as `-logfile`, `-conffile`, `-dryrun`, `-xmlhelp`, or `-version-full`.
- **Error inputs**: inputs where NCBI's Python script exits with a Python traceback (status `1`) are listed in [What LOSAT refuses](../refusals.md#cleanup-blastdb-volumespy).
