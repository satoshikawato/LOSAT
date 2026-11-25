# Coming from NCBI BLAST+

LOSAT is made to replace NCBI BLAST+ 2.17.0 without changing your commands. You can call it as `losat <program> ...`, or make the NCBI program names (`blastn`, `makeblastdb` and so on) run LOSAT, so that scripts, workflow tools and other programs that look for `blastn` on your `PATH` need no changes.

## Call a program through `losat`

Put the program name after `losat` and keep the rest of the command line:

| NCBI BLAST+ | LOSAT |
|:---|:---|
| `blastn -query q.fa -db nt_db -outfmt 6` | `losat blastn -query q.fa -db nt_db -outfmt 6` |
| `makeblastdb -in p.faa -dbtype prot` | `losat makeblastdb -in p.faa -dbtype prot` |
| `blastdbcmd -db p_db -info` | `losat blastdbcmd -db p_db -info` |

The programs are `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, `psiblast`, `rpsblast`, `rpstblastn`, `makeblastdb`, `blastdbcmd`, `dustmasker`, `segmasker`, `windowmasker`, `convert2blastmask`, `makembindex`, `makeprofiledb`, `blastdb_aliastool`, `blastdbcheck`, `cleanup-blastdb-volumes.py` and `rpsbproc`.

## Use the NCBI names: `losat --install-shims`

`losat --install-shims [DIR]` makes the 20 program names above in `DIR`. Each name is a link to `losat`, and an executable called by one of these names runs that program. Without `DIR`, the names are made in the folder that holds `losat`.

<!-- example: from-ncbi-blast-1 -->
```bash
losat --install-shims ./ncbi-names
```

The command prints the path of each name it made. A few details:

- **Links and copies.** On Linux and macOS each name is a symbolic link. In the folder of `losat` the link is relative, so you can move the folder. In another folder the link holds the full path of `losat`. On Windows, where links need a privilege, each name is a copy of `losat.exe`.
- **It never overwrites a program that is not LOSAT.** If a file of one of the 19 names is already in `DIR`, for example NCBI's own `blastn`, the command stops, changes nothing and exits with status 1. It lists the files in the way.
- **It is safe to run again.** Names that already link to this `losat` are kept. Names made by an earlier LOSAT, for example after an upgrade or a move, are replaced.
- **`-install-shims` with one dash** is accepted too, because `losat -help` shows it that way.

<!-- example: from-ncbi-blast-2 -->
```bash
mkdir -p tools
touch tools/blastn
losat --install-shims tools
```

The example above stops with an error like this, and leaves `tools` as it was:

```text
losat --install-shims: these files exist and are not this executable; nothing was installed: /path/to/tools/blastn
```

To use the names, put `DIR` on your `PATH`. Put it before NCBI's directory if both are installed and you want LOSAT to answer:

```bash
export PATH="$PWD/ncbi-names:$PATH"
blastn -query query.fna -subject subject.fna -outfmt 6
makeblastdb -in subject.faa -dbtype prot -out subject_db
```

Check which program answers with `-version`. LOSAT prints its own name and version in the first line, then a `Package:` line that names the BLAST+ release it follows:

```text
blastn: LOSAT 0.2.0
 Package: LOSAT 0.2.0, output compatible with blast 2.17.0
```

The link points at the `losat` file itself. If you move or replace `losat`, run `losat --install-shims DIR` again.

## What stays the same

- **Options.** LOSAT reads the options of each program as NCBI BLAST+ 2.17.0 does: single-dash names such as `-query`, `-subject`, `-db`, `-evalue`, `-outfmt` and `-num_threads`, the same values, and the same defaults. `-help` and `-h` open with NCBI's `USAGE` and `DESCRIPTION` headings. `-version` answers in NCBI's two lines, with LOSAT's version where NCBI writes its own.
- **Output.** On the release test set, the output is NCBI BLAST+ 2.17.0's, byte for byte, in every output format a program writes, apart from the line that names the program's version (`BLASTN (LOSAT 0.2.0)` where NCBI writes `BLASTN 2.17.0+`; see [Differences from NCBI BLAST+](differences.md#the-version-line-and-the-error-prefixes-name-losat)). Scores, E-values, coordinates and ordering are the same, and so are NCBI's error messages and exit codes for errors in input files and for errors raised after the arguments are read, apart from the prefix (`LOSAT engine error:` where NCBI writes `BLAST engine error:`). A mistake in the argument syntax (an unknown option, a missing value, a value out of range) prints NCBI's USAGE text and NCBI's error lines and exits with status 1, as NCBI does, apart from LOSAT's words (the second error line has no `(CArgException::eX)`); see [Differences from NCBI BLAST+](differences.md#argument-errors-and-help). Database files written by `losat makeblastdb` have the same sequence and header bytes as NCBI's; the index and metadata files differ only in their creation date (compared with NCBI's `makeblastdb` on Linux x86-64).
- **Databases.** NCBI databases work with LOSAT, and LOSAT databases work with NCBI BLAST+. LOSAT finds a database as NCBI does: the current directory, then `BLASTDB`, then `.ncbirc`. See [BLAST databases](databases.md).
- **Threads.** `-num_threads` changes the speed, not the output. LOSAT writes what NCBI writes with `-num_threads 1` (apart from the version line, which names LOSAT), at every thread count. `-mt_mode` is accepted and does not change the output: `blastn`, `blastp`, `blastx` and `tblastn` take 0, 1 or 2, `rpsblast` and `rpstblastn` take 0 or 1, and `tblastx` and `psiblast` have no `-mt_mode`, as in NCBI BLAST+ 2.17.0.

## What differs

LOSAT keeps to NCBI's results. A small set of differences is on purpose, and each has a written decision. [Differences from NCBI BLAST+](differences.md) lists them. In short:

- A local `-subject` search of `tblastx` and `tblastn` translates the subject with the genetic code that `-db_gencode` selects. NCBI's local search uses code 1. `tblastn` also accepts genetic code 32.
- On macOS, the last (15th) digit of a full-precision E value can differ by 1.
- Where NCBI's result depends on undefined behaviour in its code, such as a signed overflow, LOSAT computes the value without it, or refuses the input.
- In SAM output (`-outfmt 17`), the `CL:` tag of the `@PG` line holds LOSAT's own command line.
- Where NCBI writes its version, LOSAT writes its own name and version (`BLASTN (LOSAT 0.2.0)` where NCBI writes `BLASTN 2.17.0+`), and its error prefixes say `LOSAT` (`LOSAT engine error:` where NCBI writes `BLAST engine error:`).

LOSAT also refuses what it does not implement. It stops with an error instead of running something else. [What LOSAT refuses](refusals.md) lists every option and value, for example `-remote`, `-entrez_query`, `-gilist` and `-ipglist`. The NCBI programs and helper scripts that LOSAT does not provide are listed in the [Introduction](introduction.md).

If a command that works with NCBI BLAST+ fails or gives different output with LOSAT, and the case is not in those two pages, please report it; see [Troubleshooting](troubleshooting.md).
