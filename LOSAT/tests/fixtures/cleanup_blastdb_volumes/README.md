# cleanup-blastdb-volumes.py fixtures (row RG-D2c)

What they are: NCBI BLAST+ 2.17.0 `cleanup-blastdb-volumes.py` runs that
`tests/cleanup_blastdb_volumes_fixtures.rs` repeats with `losat cleanup-blastdb-volumes.py`, each
row in a fresh copy of `in/`.

## Origin

- Program: official NCBI BLAST+ 2.17.0 `cleanup-blastdb-volumes.py` (run under Python 3.13,
  specifically Python 3.13.3; LOSAT matches NCBI's script under Python 3.13) and `makeblastdb`
  (`ncbi-blast-2.17.0+-x64-linux`), run on 2026-11-25 in a fresh directory holding `in/`, with
  `LC_ALL=C TZ=UTC`, standard input `/dev/null`, and a clean `HOME` outside the directory unless
  the row sets `HOME` or `NCBI` explicitly.
- Inputs (`in/`):
  - `p5`: version 5 multi-volume protein database (`makeblastdb -blastdb_version 5`) with extra
    volumes `p5.09` (covering `.p??` extensions, `.tar.gz.md5`, and preserved non-matching files
    `.pa`, `.pinx`, `.pdb-lock`, `.tar.gz`, `.txt`, `.nin`) and `p5.99`, plus unrelated
    `other.09.pin`.
  - `p4`: version 4 multi-volume protein database with no extra volumes.
  - `n4`, `n5`: version 4 and version 5 multi-volume nucleotide databases with extra volumes
    `n4.09` and `n5.09` (and preserved `.na`, `.ninx`, `.pin`).
  - `sv`: single-volume protein database (`.pin` only, no `.pal`).
  - `quoted`, `multi_dl`, `no_dl`, `sub/sdb`, `envdb/edb`, `dots`, `globmeta/a[0]`: alias and path
    layouts covering quoted `DBLIST`, multiple `DBLIST` lines, alias without `DBLIST`,
    subdirectories, `BLASTDB`, `.ncbirc` / `ncbi.ini` (`cfg_*`, with `<RUN>` expanded to the run
    directory), dot-only / multi-dot volume names, and glob metacharacters.
- Rows (`cases.tsv`, tab-separated `name`, `env`, `cwd`, `args`):
  - `H01`–`H09`: `-h`, `--help`, `-help`, `--hel`, help precedence with unknown flags, and
    repeated-`h` ignored explicit arguments (`-hh=x`, `-hh=`, `-hh-x`).
  - `P01`–`P06`, `N01`–`N05`: protein and nucleotide v4/v5 dry-run, verbose dry-run, and deletion.
  - `S01`–`S04`: single-volume without alias, wrong molecule type, missing databases.
  - `D01`–`D09`: `DBLIST` quoting, multiple `DBLIST` lines, missing `DBLIST`, subdirs, `./`, dots,
    and glob metacharacters.
  - `E01`–`E10`: `BLASTDB` (single, trailing `/`, colon-separated), `.ncbirc` and `ncbi.ini` in
    working directory, `HOME`, `$NCBI`, lowercase `[blast]`, and `[DEFAULT]` + `[BLAST]`.
  - `A01`–`A21`: missing required arguments, invalid choices, option-like values, `=` form,
    negative numbers, `-`, unrecognized arguments, `--`, ambiguous prefixes (including `-=x` and
    `-=`), explicit arguments to flags, unambiguous abbreviations, and duplicate `-dbtype`.
  - `W01`–`W06`: terminal-width wrapping of `-h` and argument-error usage with `COLUMNS=40`, `60`,
    `120`, and `200`.
- `-version` is tested directly in `tests/cleanup_blastdb_volumes_fixtures.rs` and
  `tests/entry_names.rs` because LOSAT prints `cleanup-blastdb-volumes.py LOSAT 0.2.0` where NCBI
  prints `cleanup-blastdb-volumes.py 1.0`.
- Per-volume file order: both NCBI's Python script (`glob.glob` / `os.scandir`) and LOSAT emit
  `.p??`/`.n??` files within an extra volume in filesystem `read_dir` order before `.tar.gz.md5`.
  Because each test copies `in/` to a fresh directory whose `read_dir` order can differ from the
  recording machine's, the fixture test checks exact volume order, exact match to the test copy's
  `fs::read_dir` order followed by `.tar.gz.md5`, exact match to NCBI's frozen per-volume line
  sets, exact stderr, exact exit status, and exact file additions / deletions (`ncbi/<row>.files`).
