# rpsbproc fixtures (PS16a, PS16b)

Inputs and frozen NCBI outputs for `losat rpsbproc` (NCBI rpsbproc 0.5.1, the 2025-05-13 package).

## Provenance

- Archives (`arch/*.asn`): NCBI BLAST+ 2.17.0 `rpsblast -outfmt 11`, one thread, against CDD v3.22
  (`Cdd_LE.tar.gz` of 2026-10-27, run in its directory as `-db Cdd`), queries from E. coli K-12 MG1655
  (`GCF_000005845.2_ASM584v2_protein.faa`). The commands are in `ARCH_COMMANDS.tsv`. `lcl.asn` searches a profile
  database made by `makeprofiledb` whose ids are `lcl|` (LOSAT's rpsblast fixture `prof_rps`). `bin.asn` is the binary
  archive (`ARCHIVE_FORMAT=bin`), `one.xml` is `-outfmt 5`, `garbage.txt` is not an archive.
- `arch/bad/*.asn`: `one.asn` edited by hand (truncated, an unknown member, a duplicate member, a trailing comma, no
  `results`, no `alignments`, comments, CRLF line ends, two objects without a line end between them, and so on).
  Keep their bytes as they are (CRLF included).
- `data/`: the CDD v3.22 data files of NCBI's FTP site (taken 2026-11-11; `family_superfamily_links` as rewritten on
  2026-11-04), cut to the rows that the archives reach. `cdtrack.txt` is empty, as NCBI publishes it for v3.22.
  For every case, NCBI's output with this cut data equals its output with the full files.
- `oracle/<case>.{out,err,rc}` (and `.outfile` for `-o`): NCBI's rpsbproc 0.5.1 Linux x64 binary of the 2025-05-13
  package, run in a directory that holds `arch/` and `data/`, with the arguments of `CASES.tsv` and stdin from
  `/dev/null` unless the case names a file. Bioconda's `rpsbproc 0.5.1` (built on the BLAST+ 2.17.0 toolkit) gave the
  same bytes for every case.
- The cases of the PS16a review (2026-11-11; the same two NCBI binaries, same bytes), run as above:
  - `root/*.asn`, which the tests copy to the run directory itself (the report's header names the `-i` file as given):
    `bracket_member.asn`, `neg0_length.asn`, `org_descr.asn` and `nolength.asn` are `arch/one.asn`, an edited copy of
    it and `arch/one.asn` again (a member id in brackets, `length -0`, an `org` descriptor, no `inst.length`; oracle
    requests O3-2, O3-1, O3-6, O3-9 of the review's angle 3). `q_tab_utf8.asn` is NCBI rpsblast 2.17.0's own archive
    (`-db Cdd -evalue 0.01 -outfmt 11`, CDD v3.22) of two copies of `NP_414783.1` whose definition lines hold a TAB and
    a UTF-8 letter (O3-10).
  - `data_with/<dir>/`: the data files of `data/` that the case `d_<name>` edits; the tests make `<dir>/` from `data/`
    and these files (oracle request A4-1 of the review's angle 4): `d_acc0` adds a `cddid.tbl` row with an empty
    accession, `d_wsline` a line of one space, `d_annot10` a `cddannot.dat` line whose coordinates are its last
    field, `d_bs12` writes the threshold of `smart00421` as `0x10`, `d_bsnul` writes the bit score of `CHL00113` as
    `1x`, a 0 byte and `y` (NCBI quotes the string up to the 0 byte, exit 2), and `d_ct3`'s `cdtrack.txt` is the line
    `smart00421` alone.
  - `xml_missing_input`: `-i` names a missing file and the standard input is XML (oracle request M4 of the review).

## CASES.tsv

`case`, `use` (`compare`: LOSAT must give the same stdout, stderr and exit code; `linux`: the same, on Linux only; `refuse`: LOSAT refuses, exit 1, with
the bytes of `refused/<case>.{out,err,rc}`: the XML input until the XML reader of a later item, PS16d, and the
refusals of PS16b below), `stdin` (`-` = empty), `setup`
(`data-without:<file>`: a directory `data_missing/` holding the data files except that one; `data-with:<dir>`: a
directory `<dir>/` holding the data files with those of `data_with/<dir>/` in their place; `data-with:<a>+<dir>`: the
same with those of `data_with/<a>/` first, then those of `data_with/<dir>/`), `arguments` (split on
spaces). PS16a's cases of `-t feats`, `-t both` and the default `-t` (named `refuse_*` until PS16b) are `t_feats`,
`t_both` and `t_default`.

## Sites and motifs (PS16b)

- Archives added for PS16b (the same NCBI rpsblast 2.17.0 run, `-evalue 0.01`, commands in `ARCH_COMMANDS.tsv`):
  - `mstate.asn` (`BATCH_SIZE=1`: three sessions): NP_414781.1 and NP_416506.1 (cd08071, "Zn binding site", motif
    `HHD`), then NP_416542.1 (cd00438, "catalytic site", motif `HKY`). NCBI's `MotifCheck` keeps the positions of the
    first motif site it checks in the process (a function-static list) and checks every later motif site against them:
    the second query's `HHD` site prints, the third query's `HKY` site does not, across sessions.
  - `mfirst.asn`: NP_416542.1 alone; its `HKY` site prints.
  - `gen2.asn`: NP_416215.1 and NP_415978.1, whose hits have generic sites of hierarchy roots (with `ctprev`).
  - `motdom.asn`: NP_416435.1, a `MOTIFS` line of a curated domain (a source PSSM-ID other than 0).
  - `sdsite.asn`: NP_416698.1, sites and motifs of a structural domain (`sd`, source 0).
  - `sd2.asn` (PS16a): its second query NP_418053.1 has structural-domain hits and prints no `SITES` block.
- `data/` stays PS16a's. The PS16b archives read data sets of their own under `data_with/`, cut by the same rule
  (the rows the named archives reach), so that no file passes 64 KB:
  - `sites/`: the PS16b archives; `cdtrack.txt` empty (from `data/`), as NCBI publishes it for v3.22.
  - `ctprev/`: the PS16b archives and `three.asn`, with `cdtrack.txt` = their rows of the CDD v3.23 preview
    `cdtrack.txt` (`pub/mmdb/cdd/preview/cdtrack.txt`, 2026-11-04) and the `cddid.tbl` rows of the parents and roots
    those rows name: the `Generic` sites.
  - `rootpid/`: `gen2.asn` as `ctprev/`, but the parents and roots that `gen2.asn` does not hit (cd04410) are
    renumbered to PSSM-ID 100,000,000 + their own: a Generic site whose root has a PSSM-ID of 100,000,000 or more
    prints the root's accession instead of the number.
  - `badmotif/`: `mfirst.asn` as `sites/`, with the motif of cd00438 "catalytic site" written `K#R`, which NCBI's
    PROSITE parser refuses (exit 3, the output ends inside `SITES`).
- For every PS16b case, NCBI's output with the cut data equals its output with the full files (`ctprev/`: the full
  v3.22 files with the preview `cdtrack.txt`), and Bioconda's rpsbproc 0.5.1 gave the same bytes as NCBI's binary.
- The cases of the PS16b review (2026-11-12; oracle requests R1, R2, R5, R6, R20, R21; NCBI's binary and Bioconda's
  gave the same bytes), run as above. `data_with/rd1/`, `rd5/`, `rd4/`, `rd13/` and `rd14/` hold only `cddannot.dat`,
  that of `sites/` with the motif of cd00438 "catalytic site" (`HKY`) replaced; the cases lay them over `sites/`
  (`data-with:sites+rd1`):
  - `rd1` `Hx(0,10)K` and `rd5` `H{P}(0,10)K`: NCBI's match reads the residues `HKY` of the site beyond their end
    (the query has 185), but every byte there gives the same answer: the site prints.
  - `rd4` `K` 0xE9 `R`: a parse error whose message quotes the byte as it is (exit 3).
  - `rd13` `Hx(99999999999999999999)K` (`atol` gives `LONG_MAX`) and `rd14` `Hx(200000000000000)K`: the minimal map
    passes the string's `max_size()` (`basic_string::_M_replace_aux`) or needs 2e14 bytes (`std::bad_alloc`); exit 3,
    the report cut as for a parse error (NCBI run under `prlimit --as=4294967296`).
  - `arch/mfirst_delta.asn`: `mfirst.asn` with its query as `repr delta`, the residues in two `ext delta` literals
    (NCBIeaa, 100 and 85 residues); case `both.mfirst_delta` with `sites/`.
- The cases of the second PS16b review (2026-11-12, oracle O2 on ws: NCBI's binary and Bioconda's 0.5.1, the same
  bytes), as above: `rd15` `Hx(4611686018427387903)K` and `rd18` `Hx(9223372036854775806)K`, minimal maps of 2^62 and
  2^63 - 1 letters, which NCBI's C++20 builds allocate (`max_size()` is 2^63 - 1) and fail with `std::bad_alloc`
  (exit 3, the report cut as for a parse error).
- `refused/<case>.{out,err,rc}`: LOSAT's own bytes (written from LOSAT's output and checked by hand: standard output is
  NCBI's up to the line before the refused site or query, standard error ends with the refusal). The refusals of
  PS16b (PD-LOSAT-NCBI-DEFECTS 1.15, explicit rejections):
  - `refuse_rd2.mfirst`: `data_with/rd2/` `Hx(0,10)W` (review oracle R4): NCBI's match reads beyond the residues
    `HKY` and its answer depends on the bytes there (NCBI's stack: the review's runs dropped the site, but no line of
    NCBI's code decides it).
  - `refuse_rd19.mfirst`: `data_with/rd19/` `H-{<}(0,67108863)-{<}(0,67108863)-K` (review 2 oracle O3): NCBI holds
    1.5 GiB of match records at once (exit 0 under a 4 GiB address space, `std::bad_alloc` under 1.5 GiB).
  - `refuse_rd20.mfirst`: `data_with/rd20/` `Hx(536870910)K` (O3): the minimal map's string of 512 MB grows to
    1 GiB while the old buffer lives (as rd19 on NCBI).
  - `refuse_nodata`: `arch/q_nodata.asn`, `mfirst.asn` without the query's `seq-data` (review oracle R19), with
    `sites/`: the site's positions are outside the (empty) residues, where NCBI reads memory outside the string
    (its two builds print different bytes).

## The NCBI C++ Toolkit options (line RG, item E1c)

- The rows `e1c_*` (2026-11-13): `-dryrun` (rpsbproc's own `DryRun` after its `Init`), `-logfile log.txt` and
  `-conffile c.ini`. `oracle/<case>.{out,err,rc}` and `.outfile` / `.log` (the `out.txt` and `log.txt` the run
  made): NCBI's rpsbproc 0.5.1, the binary of the 2025-05-13 package and bioconda's (on the BLAST+ 2.17.0 toolkit),
  which gave the same bytes, run as the test runs LOSAT (the run directory above, no `LOSAT_` or `NCBI_` variable,
  `NCBI_CONFIG_PATH` = the run directory, stdin from `/dev/null` or the named file).
- `conf/`: the registry files of the setups `conf:<file>` (copied as `c.ini`), `ini:<file>` (as `rpsbproc.ini`) and
  `ncbirc:<file>` (as `.ncbirc`): `blast.ini` is `[BLAST]` alone; `garbage.ini` a line that is not an entry,
  `sec.ini` an unclosed section, `quote.ini` a value with a quote inside, `cont.ini` a continued value then a line
  that is not an entry (NCBI's report, exit 2); `bom2.ini` a UTF-8 byte order mark and one byte more (NCBI's
  `Error: (110.4) Error reading the registry after line 1: `, the run goes on); `bom1.ini` the first byte of one (as a
  `.ncbirc`: the same post; as a `.ini` NCBI posts nothing).
- `arch/bad/subj_sequences.asn`, `matrix_integer.asn`, `evalue_rawscore.asn` (the review of item E1c): `arch/one.asn`
  edited by hand: `subject sequences { }`, `MatrixName` `value integer 5`, `EvalueThreshold` `value cutoff raw-score
  5`. NCBI writes the `SESSION` line's tab before each getter throws; with `-o log.txt -logfile log.txt` its stream
  keeps the part of that line until the program ends, and writes it over the start of the exception's report.
  `matrix_integer_longdb.asn` is `matrix_integer.asn` with a database name of 5,000 bytes, which libstdc++'s buffer of
  8,192 bytes keeps (an insertion is written directly only when it does not fit in what is left);
  `matrix_integer_db8165.asn` one of 8,165 bytes, after which the buffer holds 8,190 bytes and takes the `char` of the
  next tab (`put`) without writing (the review 2 of item E1c).
- `use` `linux` (`linux_cout_longdb`): run only on Linux (`rpsbproc_gives_ncbis_bytes_through_a_pipe_on_linux`):
  `-logfile /dev/stdout` with the standard output a pipe, where glibc's buffer of `cout` (4,096 bytes) decides the
  order of the report and the diagnostics.
- `refused/refuse_xml.*`: LOSAT's bytes for the XML input with a successful `Init`, refused when its run would read the
  archive (NCBI's standard error up to its warning `Warning: XML format support is deprecated!`, then the refusal).
  With `-dryrun` NCBI reads no archive and passes, and LOSAT gives its bytes (`e1c_dry_xml`).
