# A mask file whose Seq-locs nest deeper than NCBI's stack holds (PD-LOSAT-NCBI-DEFECTS exception 28)

Line NP item M (b), 2026-11-23. NCBI BLAST+ 2.17.0 (the bioconda build and the official Linux x86-64 release) and
LOSAT, `makeblastdb -mask_data`.

## What NCBI does

makeblastdb reads a mask file (`Blast-db-mask-info`, ASN.1 text or binary) with NCBI's serial library
(`CObjectIStreamAsn`, `CObjectIStreamAsnBinary`). A Seq-loc is the one type of a mask file that holds values of its
own type: `mix` (SEQUENCE OF Seq-loc), `equiv` (SET OF Seq-loc), and an element written as `:Seq-loc` in NCBI's
pointer syntax. NCBI reads it by recursion on the main thread's stack. Each level opens the choice, its variant, the
list and the list's element (`c++/src/serial/objistr.cpp:1482-1500`, `c++/src/serial/objistrasn.cpp:1245-1274`),
which takes about 1.3 KB of stack.

With the default stack limit of 8 MiB, NCBI ends with SIGSEGV on a nest of about 6,500 to 6,900 levels:
- exit status 139;
- stdout holds only the log's first lines ("Building a new DB ... Maximum file size");
- stderr is empty;
- no database is written.

The table gives each nesting path (5 runs per depth, minimal environment): the deepest nest read in all 5 runs and
the shallowest from which every run crashed.

| nest (`gen_mask_nests.py` kind) | bioconda, binary | bioconda, text | official, binary | official, text |
|---|---|---|---|---|
| `mix` | 6,449 / 6,473 | 6,622 / 6,634 | 6,785 / 6,809 | 6,874 / 6,898 |
| `equiv` | 6,450 / 6,474 | 6,610 / 6,634 | 6,785 / 6,809 | 6,873 / 6,897 |
| `mix` and `equiv` alternating (`alt`) | 6,449 / 6,473 | 6,613 / 6,637 | 6,784 / 6,808 | 6,875 / 6,899 |
| `mix` of `:Seq-loc` elements (`ptr`) | 6,452 / 6,464 | 6,612 / 6,636 | 6,785 / 6,809 | 6,872 / 6,896 |

A bisect found a first crash depth for each path; the depths up to 60 levels either side of it, in steps of 12, ran
5 times each. The depths in between crash in some runs only (for example bioconda binary `mix` at 6,461 levels in 4
of 5 runs, official text `mix` at 6,886 in 3 of 5). The official build's table is the third of four surveys
(2026-11-23 10:49, 20:53 and 22:07, 2026-11-24 02:32 JST); the other three put each boundary within 11 levels of
it. Several things move the crash depth: the build, the stack limit, the lengths of the command line and the
environment, and the stack's random start at each run.

With no stack limit (`ulimit -s unlimited`), NCBI reads every nest tried, up to 100,000 levels, and writes its
normal database (at 100,000 levels in 0.2-0.55 s, with 205-239 MB of memory for the bioconda build and 170-200 MB
for the official one). For a file damaged deep inside a nest, NCBI writes its report, which names every frame of the
nest. These reports grow with the depth; at 100,000 levels they take about 25 s and are 2 MB long.

A number too large for `int` in a text file gives a different report: one line per frame. It took 0.5 s at 401
levels, 245 s at 5,000 levels, and 4,933 s (253 MB of memory) at 20,000 levels, for a report of 80,014 lines and
17 MB.

## What LOSAT does

LOSAT keeps the levels and frames on the heap (`TypedReader::nest`, `LOSAT/src/report/serial/asn_types/typed.rs`),
and the Seq-id walks of the levels (`Build::seq_loc_parts`, `LOSAT/src/writedb/mask_info/reader.rs`). It reads every depth as far as memory allows, with an 8 MiB or a 2 MiB stack, and writes
what NCBI writes with no stack limit: the exit status, the log, the database and, for a damaged file, the report.
Its reports take time linear in the depth.

- In all 90 rows of `results.tsv` LOSAT's exit status, standard output, database and the nest frames of its report
  equal those of NCBI's bioconda run with no stack limit; the 60 rows run with the official release too equal
  that build's.
- At 100,000 levels LOSAT takes at most 0.12 s and 34 MB, where NCBI takes 0.2-0.55 s and 170-239 MB; a damaged
  file at most 0.13 s and 24 MB, where NCBI takes 24-28 s and 205-210 MB.
- The number too large for `int` at 20,000 levels takes LOSAT 0.23-0.38 s and 51-53 MB. Its report has the same
  message and the same 80,012 frame lines in the same order as NCBI's (`compare_reports.py`: NCBI's location
  prefixes removed, LOSAT's words through `crate::wording`).

Dev before this change refused every nest of 97 levels and more ("objects nested more than 400 levels deep").

## The runs

`results.tsv` has one row per kind, format and depth (LOSAT's columns: the build of commit bb678eff6). Each row
gives:
- the exit statuses of NCBI's 5 runs with 8 MiB, for the bioconda build and, for the kinds that read, the official
  release;
- NCBI's run with no stack limit, dev's run, and LOSAT's runs with 8 MiB and with 2 MiB, each as
  `rc/signal/wall s/max RSS KB`;
- whether LOSAT's exit status, standard output, database and nest frames in standard error equal those of NCBI's
  run with no stack limit.

The kinds: `mix`, `equiv`, `alt` and `ptr` as above; `wide` (two intervals per level); `wideid` (a second Seq-id,
so that the walk stops); and damaged files: `badvar` (an unknown variant innermost), `badstrand` (strand 7),
`bignum` (a `from` past `int`) and `cut` (the file ends inside).

LOSAT's report is NCBI's in LOSAT's words (`crate::wording`). The byte comparison through that mapping is the fixture
rows `n_*` of `LOSAT/tests/seqdb/mask_parse/expected.tsv`, frozen from the official build's runs with no stack limit.

Inputs: `python3 -I gen_mask_nests.py KIND FORMAT DEPTH > m.asn`. Each run used:

```sh
env -i PATH=/usr/bin:/bin HOME=<empty dir> BLAST_USAGE_REPORT=false TZ=UTC \
  prlimit --as=8589934592 --core=0 --stack=<8388608 | unlimited | 2097152> --cpu=1200 \
  makeblastdb -in mask_nucl.fna -dbtype nucl -out db -mask_data m.asn
```

`mask_nucl.fna` is `LOSAT/tests/seqdb/mask_nucl.fna`. NCBI's run with no stack limit of `bignum` text at 20,000
levels had `--cpu=14400` (it did not end within 1,200 s). The database holds its build time as text, and the
comparison masks the time but not its length: LOSAT's runs used a time zone (`TZ=UTC+6`) whose hour had one digit,
as in NCBI's runs.
