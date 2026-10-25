# Product Decision: TBLASTX local-subject genetic code

- Decision ID: `PD-TBLASTX-LOCAL-SUBJECT-GENCODE`
- Version: 1.0
- Date: 2026-10-25 (this record). The deviation itself was approved by the Owner at the
  TBLASTX v0.1.0 certification (`docs/release/tblastx_v0.1.0_certification.md`, recorded
  2026-08-30; "Approved local-subject `db-gencode` deviation") and is kept for v0.2.0 by the
  Owner on 2026-10-25.
- Status: Accepted (approved deviation, kept for v0.2.0). This record only writes down what
  `PD-PURE-RUST-RUNTIME-AUTHORITY` (1.1, "Existing compatibility decisions"), the root
  [`AGENTS.md`](../../AGENTS.md) (bit-perfect rule, approved project exception) and
  `docs/release/v0.1.0.md` ("TBLASTX approved deviation") already name; it widens nothing.

## Scope

`tblastx` only. The TBLASTN counterpart is
[`PD-TLOSAN-LOCAL-GENCODE-32`](PD-TLOSAN-LOCAL-GENCODE-32.md), which is confined to TBLASTN
and does not change TBLASTX.

## Decision

With `-subject` and a `-db_gencode` other than 1, LOSAT `tblastx` translates the subject with
the selected genetic code and uses that translation for the search, the scores and the
report. NCBI BLAST+ 2.17.0 translates a local `-subject` with genetic code 1 whatever
`-db_gencode` says, and only hands the selected code to the report formatter. LOSAT does not
reproduce that, because the NCBI behaviour searches a sequence with a code the user did not
ask for.

### Exact condition

All of these must hold; otherwise LOSAT equals NCBI byte for byte:

1. the program is `tblastx`;
2. the subject is a local FASTA file given with `-subject` (not `-db`);
3. `-db_gencode` is a value other than 1 (the accepted values are those of NCBI's CLI,
   `genetic_code` in `LOSAT/src/blastinput/value_parsers.rs:728`: 1-6, 9-16, 21-31, 33).

`-db` searches with a non-default `-db_gencode` are identical to NCBI's (frozen rows
`db.tblastx.gc4.6` and `db.tblastx.gc11.6` in `LOSAT/tests/db_search_manifest.tsv`, lines 125
and 126, and `db.sa.tblastx.gc.15` / `.16`, lines 332 and 333). The default code 1 and every
`-query_gencode` value are identical to NCBI's (case `p14_ap027131_ap027133_query4`,
`-query_gencode 4 -db_gencode 1`: raw-byte identical, 14,871 rows;
`docs/release/tblastx_v0.1.0_certification.md`, "Approved local-subject `db-gencode`
deviation").

### NCBI behaviour and where it comes from

NCBI 2.17.0 source (the tree at `ncbi-blast-2.17.0+-src/c++/src/`):

- `algo/blast/blastinput/blast_args.cpp:2547-2557` reads a local subject with
  `ReadSequencesToBlast` into a `CObjMgr_QueryFactory`; nothing there sets a genetic code.
- `algo/blast/api/blast_setup_cxx.cpp:800-810` takes the subject's genetic code from
  `subjects.GetGeneticCodeId(i)`, and `include/algo/blast/api/sseqloc.hpp:157-173,203-208`
  initialise that value to `BLAST_GENETIC_CODE` (1) for a sequence read this way.
- `app/blast/tblastx_app.cpp:155-156` gives the formatter `opt.GetQueryGeneticCode()` and
  `opt.GetDbGeneticCode()`, so the report is formatted with the selected code while the
  search used code 1 (the same boundary that `PD-TLOSAN-LOCAL-GENCODE-32` records for
  TBLASTN).

### LOSAT behaviour

`GeneticCode::from_id(args.db_gencode)` builds the subject translation in the search
(`LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:1735`, `:2956`, with the comment at
`:2957-2959`: "LOSAT intentionally treats `--db-gencode` as the subject genetic code for local
`-s/--subject` searches. Unlike BLAST+, search/scoring and reporting both honor the explicit
subject code") and in the reports (`LOSAT/src/algorithm/tblastx/blast_xml.rs:259`,
`LOSAT/src/algorithm/tblastx/report.rs:226,553`). The exception does not touch timing,
ordering, scoring formulas, statistics, pruning, coordinates or formatting: they follow NCBI
for whatever translation the subject has.

## Observable effect

The HSP set differs from NCBI's. LOSAT's HSPs are those of the correct translation of the
subject; NCBI's come from a code-1 translation and are displayed with the selected code. Of
the 20 cases of the gbdraw P1-P2 profile, the six that use a non-default subject code differ
(14 are raw-byte identical). The largest probe, `d06_ap027131_ap027133_db4`
(`-query_gencode 1 -db_gencode 4`), gives 12,672 rows in NCBI and 13,644 in LOSAT
(certification document, same section; LOSAT SHA-256
`719339294d0f22d0f415dd33106a0fd93c25ac372b1456d19d92390fce081df1`, also in
`LOSAT/tests/platform_native_v010_canonical.tsv:47`).

## Verification

- `LOSAT/tests/tblastx_v010_parity_manifest.tsv`: the `contract` column is
  `approved_db_gencode_deviation` for exactly `d01_nz_self_code4`, `d02_ap027132_nz_code4`,
  `d03_ap027078_ap027131_code4`, `d04_ap027131_ap027133_code4`, `d05_ap027133_ap027132_code4`
  and `d06_ap027131_ap027133_db4` (lines 13-17 and 21), and `parity` for the other 14.
  `LOSAT/tests/audit_tblastx_v010.py:41,97,114` accepts that contract only with explicit
  genetic-code arguments and a `db_gencode` other than 1, and requires one such case that is
  repeatability-checked (`:149`).
- Gate A: the six cases are classified `HSP_SET_DIFF` in
  `LOSAT/tests/platform_native_v010_canonical.tsv` (lines 39-43, 47) and
  `LOSAT/tests/platform_native_v011_canonical.tsv` (lines 42-46 and 50), and
  `LOSAT/tests/certify_platform_native_v010.py:2323-2335` fails the certification if the
  number of `approved_db_gencode_deviation` TBLASTX rows is not 6 or `d06` changes class.
  The class is declared, never promoted to exact parity.
- Outfmt 0: rows `tblastx.code4.0` and `tblastx.code4.7` of `LOSAT/tests/outfmt0_manifest.tsv`
  (lines 88 and 89) carry the same contract; LOSAT must equal the NCBI `-db` output outside
  the database lines (header comment at line 10).
- The decision is capped at those six cases (`PD-PURE-RUST-RUNTIME-AUTHORITY`); a seventh
  case with another contract value, or any other TBLASTX difference, is a defect, not part of
  this decision.

## Release-notes sentence

`tblastx -subject` with a `-db_gencode` other than 1 translates the subject with that genetic
code for the search, the scores and the report, where NCBI BLAST+ 2.17.0 searches a local
subject with genetic code 1, so the HSP sets differ by design; `-db`, the default code and
`-query_gencode` match NCBI.
