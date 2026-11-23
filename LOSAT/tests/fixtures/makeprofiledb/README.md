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

The whole oracle (248 cases) is compared by `compare.py` of the survey; these rows are the part kept in the
repository (about 1.05 MB).
