# makeprofiledb fixtures (row PS10a)

What they are: NCBI BLAST+ 2.17.0 `makeprofiledb` runs that `tests/makeprofiledb_fixtures.rs` repeats with
`losat makeprofiledb`, and NCBI `rpsblast` 2.17.0 output on two of the databases it built.

## Origin
- Program: NCBI BLAST+ 2.17.0 `makeprofiledb` and `rpsblast` (bioconda build, `makeprofiledb: 2.17.0+`, build Aug 11
  2025), run on the local development host (WSL2) by the PS10 survey oracle on 2026-11-21: `run_cases.sh` of
  `psiblast-20261025/PR/ps10/oracle/` (outside the repository), each case in a fresh directory holding `in/`, with
  `env -i PATH=/usr/bin:/bin TZ=UTC LC_ALL=C BLAST_USAGE_REPORT=false HOME=<run dir>`, under `prlimit --core=1`, the
  wall clock frozen at 2026-11-21 05:07:00 UTC by an `LD_PRELOAD` time shim.
- Inputs (`in/`): `pc.smp`, `pe.smp`, `pf.smp` are `psiblast 2.17.0 -num_iterations 2 -out_pssm` of 32-36 residue
  protein cuts (scores and frequency ratios, scale 1); `c1.smp`-`c3.smp` are CDD v3.22 cd00001-cd00003 cut to their
  first 22/18/20 columns (scale 100, word thresholds); `fc.smp`, `fe.smp` are frequency-ratio-only forms of the
  psiblast PSSMs (the IMPALA path); `t_source.smp`, `t_orgnotax.smp` are `pc.smp` with a `source` taxid and an
  organism without a taxid; `e_trunc.smp` and `e_unknownfield.smp` are `pc.smp` cut in its frequency ratios and with
  an unknown member, `e_notpssm.smp` a `Seq-entry`, `e_empty.smp` empty (the reader's errors). Added on
  2026-11-22 (the rows `R*`, run the same way that day): `nv00.smp`-`nv03.smp`, generated score PSSMs of 3-11
  columns (scale 100), and their variants written by `PR/ps10/scripts/fix_inputs.py` (same place): `z_bh0`,
  `z_bh1`, `z_bh2`, `z_bhorg`, `z_bhnoid` carry binary definition lines (a `user` object `ASN1_BlastDefLine`: an
  empty set, one line, two lines whose first has no Seq-id, beside an organism, on a query without Seq-ids);
  `z_noid_d`, `z_noid_n` have no Seq-id (with and without descriptors); `z_l2p31`, `z_lneg1`, `z_lneg3`,
  `z_lneg1r`, `z_lneg1f` give the query a length past `int` (2^31; 2^32-1 and 2^32-3 with a matching negative
  `numColumns`, by row, with frequency ratios only). Lists (`*.lst`) name the PSSMs relative to the working
  directory; `nope.smp` is a file that does not exist.
- Then `freeze_fixtures.py` (same place) kept the 65 rows of `cases.tsv` and, per row, NCBI's standard output,
  standard error and exit status as written (`ncbi/<row>.{out,err,rc}`; the test maps NCBI's words to LOSAT's with
  `losat::wording`) and `ncbi/<row>.files`: for each file the run wrote (outside `in/`), its path, size and FNV-1a 64
  digest; in a `.pin` or `.pal` the build's date (its 4-byte length and padding in a `.pin`) is replaced by `<DATE>`
  first and the size is `*`.
- `rps/q3.faa` is the survey's three-protein query; `rps/R01.{6,0}` and `rps/R05.{6,0}` are NCBI rpsblast
  `-query q3.faa -db db/db -outfmt 6|0` on NCBI's databases of rows `O00` (score PSSMs) and `S_f2`
  (frequency-ratio-only PSSMs).
- Rows `QC01`-`QC07` (row PS10c, run the same way on 2026-11-23 by `psiblast-20261025/PR/ps10c/fix/run_cases.sh`,
  inputs from `PR/ps10c/scripts/make_fixture_inputs.py`, rows added by `PR/ps10c/scripts/freeze_ps10c.py`):
  `fc.smp` cut to its first 12 columns with its query written as letters or codes, the query of a frequency-ratio
  profile, which NCBI converts with `CSeqConvert`'s tables. `q_il.smp`: IUPACaa in lower case (the capitals'
  codes; `rps/QC01.{6,0}` is NCBI rpsblast on its database); `q_mp.smp`, `q_ml.smp`, `q_ms30.smp`: NCBIeaa `+`,
  NCBIeaa in lower case, the NCBIstdaa code 30 (each no residue), one letter short of the query's length (NCBI's
  `Query sequence lengths mismatch` comes first); `q_xp.smp`: NCBIeaa `+` with the matrix `FOO` (the engine's
  unknown matrix comes first). The same queries of the right length make NCBI's engine write outside its table
  of residues (then a crash, an abort, or a database, by the code), which LOSAT refuses. `q_ng.smp`, `q_ngi.smp`:
  a negative frequency ratio (exit status 1), with BLOSUM62 and with IDENTITY 15/2 (NCBI's error comes before
  IDENTITY's missing table is read).
- Rows `QD01`-`QD46` (row PS10d, run the same way on 2026-11-23 by `psiblast-20261025/PR/ps10d/fix/run_cases.sh`,
  twice, with no difference; inputs from `PR/ps10d/scripts/make_fixture_inputs.py`, rows added by
  `PR/ps10d/scripts/freeze_ps10d.py`): a profile's query in a nucleotide coding, which NCBI's writer stores as
  the sequence's bytes. `d1.smp`, `d2.smp`, `d3.smp` are `c1.smp`-`c3.smp` cut to 6, 5 and 8 columns (without
  `weightedResFreqsPerPos`), `e1`, `e2` the first two without `finalData` (frequency ratios only); a suffix
  `_na`, `_n2`, `_n4` writes the query as `iupacna`, `ncbi2na` or `ncbi4na` (bases `ATGGCG`, `GATTA`,
  `CCGATAGT`). `QD01`-`QD15`: such a query alone, first, in the middle, last, and all three; `QD16`-`QD27`: in a
  profile with frequency ratios only, alone, in the middle, last, and with `-max_smp_vol 2` (NCBI's `Unknown
  conversion.`, exit status 255, after the files it writes); `QD28`, `QD29`: data longer than the length
  (`ib=..,n4sz=..` lines on the standard output, also with `-logfile`); `QD30`-`QD37`: the unused bits of the
  last byte set, a pad nibble, a letter more for an odd length, iupacna two and three letters short (NCBI reads
  the string's terminating zero), a length of 1 with empty iupacna (with frequency ratios only: `Query sequence
  lengths mismatch` after the line), ncbi4na with `N`; `QD38`: `mol na` in the middle of a list (the writer's
  molecule error first); `QD39`: iupacna with a two-byte UTF-8 letter (the reader's error for each byte, which it
  stores as `#`); `QD40`: `-max_smp_vol 2` with longer data in both volumes (a size line each, the alias file);
  `QD41`-`QD46`: the codings that the reader keeps as generic values, in the middle of a list (from the PS10d
  review's cases O1-O5): with frequency ratios only `ncbi8na` (`Unknown conversion.`, after the writer's error
  posted as it closes), `ncbi8aa` (its bytes are the codes; the writer's error at the next query), `ncbipna` and
  `ncbipaa` (NCBI's unknown exception), empty `ncbi8na` (the lengths' error); with scores `ncbi8na` (the writer's
  error at the next query). `l_QR01.lst`-`l_QR05.lst` and their inputs (`d1_n2_less1`, `d1_n4_less1`,
  `d1_na_less4`, `d1_n2_empty`, `e1_na_empty`) are queries shorter than NCBI's writer reads for their length:
  NCBI's two runs gave different `.psq` bytes for the first two, and SIGSEGV for empty ncbi2na;
  `l_QR06.lst` (`e2_n8` alone: NCBI's process ends as its writer closes, review O6) and `l_QR07.lst` (a gap,
  where NCBI's conversion reads past its table; not run) have no row either. The test checks LOSAT's refusal of
  each.
- Rows `NE00`-`NE27` (line RG item D3c-3c, run the same way on 2026-11-24 by
  `remaining-gaps-20261112/D/d3/c3c/tools/run_probeB3.sh`, inputs from `c3c/tools/gen_probeB3.py`, rows added by
  `c3c/tools/freeze_ne.py`, outside the repository): an element stated as `NULL` in a list of a profile's query.
  `zs.smp` is `z_bh1.smp` without its binary definition lines; `zs_d1`, `zs_d2`, `zs_d0` hold `descr { NULL }`,
  `{ title, NULL }` and `{ NULL, title }`: NCBI's loops over the descriptors dereference the element (`Error: ERROR:
  Attempt to access NULL pointer.`, exit status 255) alone, as the second profile (`NE04`: the files of the volume
  begun), as the first, with `-taxid`, `-blastdb_version 4`, `-index false`, and with `-max_smp_vol 1` (`NE09`:
  `Empty database`, exit status 1, no file). `n_bh`, `n_bh0` (`NE10`-`NE12`) are `z_bh1.smp` with `NULL` after and
  before its binary definition lines (after them only the taxid loop reaches it: `-taxid 5` gives the database);
  `e_oss_empty_null` (`NE13`, `NE14`) has an empty octet string, so the definition lines are built from the
  descriptors; `e_bh_org_null`, `e_bh_null_org` (`NE15`, `NE16`) put an organism before or after the `NULL` (the
  taxid loop stops at the first organism). `NE17`-`NE26` hold a `NULL` element in lists makeprofiledb does not
  read (Seq-hist ids and assembly, Seq-ext seg and delta, BioSource subtype, a feature's dbxref and qual, a
  Variation-ref set, population-data and clinical-test, whose deprecation NCBI posts): the database. `NE27`
  (D3c-3c review R2, `c3c/probeF/mpdb`, run the same way) is `e_noid_oss_null`: a Bioseq with an empty id
  list whose descriptors hold a title, binary definition lines of an empty octet string and a `NULL`; NCBI
  builds the definition lines from the descriptors before it looks at the ids, so it dereferences the element
  there too.

- Rows `B01`-`B22`, `BE01`-`BE12`, `RB04`, `RB04m`, `RB04b`, `RB08`-`RB29` (row PS10b, frozen on
  2026-11-25 with the official NCBI BLAST+ 2.17.0 `makeprofiledb` and `rpsblast` binaries):
  - `B01`-`B09`: `-binary` on ASN.1 BER binary PSSMs (`pc`/`pe`/`pf` scored, `c1`/`c2`/`c3` CDD,
    `fc`/`fe` frequency-ratio-only IMPALA path), multi-volume (`-max_smp_vol 2`), `-index false`,
    `-blastdb_version 4`, `-scale`/`-threshold`/`-gapopen`/`-gapextend`, `-logfile`, `-dryrun`;
    `rps/B01.{6,0}` and `rps/B03.{6,0}` are NCBI `rpsblast` outputs on the `-binary`-built DBs.
  - `B10`-`B22`: short-form definite length (`B10`, valid), long-form definite length (`B11`,
    reproducing NCBI's `error while reading freqRatios`), `byRow TRUE` (`B12`), binary defline
    inside binary `.smp` (`B13`), `UTF8String` tag (`B14`, suppressed warning), bad visible chars
    in binary title (`B15`, 803.7 posts), `StringStore` `ncbieaa`/`iupacaa` bytes `>= 0x80` and
    `< 0x20` (`B16`, `B17`), odd-length `iupacna` overhang byte `>= 0x80` (`B18`), `iupacna`
    control byte `< 0x20` (`B19`), `mol na` with `iupacna` `>= 0x80` (`B20`, NCBI molecule error
    first), `params NULL` (`B21`), BER special REAL `0` / `-0.0` (`B22`).
  - `BE01`-`BE12`: `-binary` errors (text with `-binary`, binary without `-binary`, mixed binary
    and text with leftovers, empty, truncated `30 80`, missing mandatory `pssm`, mid-file
    truncation with leftovers, unknown member `[9]`, bad tag inside `numRows`, huge string 2^63,
    huge octet string 2^62 with 803.4 post, `@1` pointer overflow).
  - `RB04`, `RB04m`, `RB04b`, `RB08`-`RB29`: damaged and edge-case binary definition lines
    (`ASN1_BlastDefLine`) decoded through the shared table-driven BER reader (`Blast-def-line-set`):
    `0xDEADBEEF` (`RB04`, mid-list `RB04m`, in binary `.smp` `RB04b`), truncated `30` (`RB08`),
    missing `seqid` (`RB09`), overflowing `taxid` (`RB10`), bad visible chars in title and Seq-id
    (`RB11`, 803.7 posts), bad visible char followed by syntax error (`RB12`), `NULL`
    `Blast-def-line` (`RB13`), `UTF8String` tag (`RB14`), huge string 2^63 (`RB15`), `:Blast-def-line`
    and `:Seq-id` pointer forms with definite length and trailing bytes (`RB16`), and `NULL`
    `Seq-id` elements carried to NCBI's exact failing point: `-blastdb_version 4` (`RB17`, `RB18`
    with `-taxid 9606`: 802.5 post + `byte 15` error, rc 255), non-empty `-taxid_map` in v4 and v5
    (`RB19`, `RB20`: `Attempt to access NULL pointer.`, rc 255), varied `NULL` positions (1st
    defline 2nd Seq-id `RB21` byte 28; 2nd defline 1st Seq-id `RB22` byte 56; 2nd defline 2nd
    Seq-id `RB23` byte 69), `-index false` replacing the first slot with `gnl|BL_ORD_ID|0` in v4
    and v5 (`RB24`, `RB25`: rc 0; 2nd Seq-id `RB26` byte 49), `org` descriptor populating the map
    in v5 (`RB27`: `Attempt to access NULL pointer.`, rc 255), v5 after a valid profile (`RB28`:
    802.5 + `byte 15`, rc 255, 9 files), and bad visible char before `NULL` `Seq-id` (`RB29`:
    803.7 + 802.5 + `byte 17`, rc 255).
  - `l_QR08.lst`-`l_QR15.lst` (no `cases.tsv` row; tested as refusals): `iupacna` byte `>= 0x80`
    in pair positions (`QR08`-`QR10`), `:Seq-interval` other-type pointer in binary (`QR11`) and
    ASN.1 text (`QR14`, where NCBI 2.17.0 prints `class not found: Seq-interval`, rc 255), `NULL`
    `Seq-id` inside `ASN1_BlastDefLine` in v5 with 0 prior OIDs (`QR12`: none, `-taxid 9606`, and
    empty `-taxid_map`, where NCBI aborts rc 134), 1 GiB string length past EOF (`QR13`), and
    `NULL` `Seq-id` followed by a missing input file in v5 and v4 (`QR15`, where `Close` throws in
    the unwinding of the input error and NCBI aborts rc 134).

- Rows `D01`-`D32` (row PS14, frozen on 2026-11-26 with the official NCBI BLAST+ 2.17.0
  `makeprofiledb` and `rpsblast` binaries; each row run twice with zero difference):
  - `D01`-`D07`: `-dbtype delta` on CDD scoremats (`c1`, `c2`, `c3`) in v5 (`D01`), v4 (`D02`),
    `-index false` (`D03`), `-binary` (`D04`, `c1.smp.bin`-`c3.smp.bin`), multi-volume `-max_smp_vol
    2` (`D05`), `-taxid 9606` (`D06`), and `-taxid_map` (`D07`); `rps/D01.{6,0}` and
    `rps/D04.{6,0}` are NCBI `rpsblast` outputs with `rps/qcdd.faa` on the built delta databases.
  - `D08`-`D09`: missing `weightedResFreqsPerPos` (`D08`, `pc.smp`, rc 1) and missing
    `numIndeptObsr` (`D09`, rc 1).
  - `D10`-`D17`: `-obsr_threshold` exact boundary `5.0 + 1.0 == 6.0` vs `4.9999 + 1.0 < 6.0`
    (`D10`, and still excluded with `-exclude_invalid false` in `D11`), `x_ValidateCd` column-sum
    epsilon `1.0001` vs `1.00011` (`D12` excluded with NCBI's `it conatins an invalid CD `, `D13`
    kept with `-exclude_invalid false`), `obsr < 1.0` at `-obsr_threshold 0` (`D14` excluded, `D15`
    kept), `-obsr_threshold -10 -exclude_invalid false` with `obsr = -0.5` (`D16`, valid `Uint4`
    `500`), and `.obsr` run-length compression boundaries relative to run head with
    `obsr.len() != seq_size` (`D17`).
  - `D18`-`D20`: all profiles excluded via `-obsr_threshold 10000` in v5 (`D18`), v4 (`D19`), and
    `-max_smp_vol -1` (`D20`): `Empty database`, rc 1, no files left.
  - `D21`-`D23`: column count mismatch (`D21`, and still rc 255 with `-exclude_invalid false` in
    `D22`), and mismatch skipped when excluded earlier by `-obsr_threshold 100` (`D23`).
  - `D24`-`D27`: frequency-ratio-only delta profile (`D24`), bad visible characters (`803.7`: 1x
    for excluded, 3x for included across NCBI's three read passes, `D25`), Pass-2 scaling mismatch
    cleaning up `.wcounts`/`.obsr`/RPS files (`D26`), and Pass-1 error after exclusion (`D27`).
  - `D28`-`D32`: `-logfile` with exclusions (`D28`), multi-volume alias (`D29`), and Pass-1 error
    (`D30`, matching `std::ofstream` buffered flushing); multi-volume split after exclusion (`D31`,
    4 inputs with 1 excluded -> 1+2 volumes); and Pass-2 duplicate `Seq-id` cleanup (`D32`).
  - `l_QR16.lst`-`l_QR22.lst` (no `cases.tsv` row; tested as refusals alongside `-dbtype cobalt`):
    `byRow TRUE` (`QR16`), `numRows 26` with weighted frequencies (`QR17`), empty `numIndeptObsr`
    (`QR18`), partial trailing column of 29 entries (`QR19`), and out-of-range / NaN floats
    (`QR20`-`QR22`: `obsr = -2`, `weightedResFreqsPerPos = -0.5`, `obsr = NaN`).

The whole oracle (248 cases) is compared by `compare.py` of the survey; these rows are the part kept in the
repository.
