# PD-TLOSAN-LOCAL-GENCODE-32

Status: approved TBLASTN v0.2.0 product requirement; this Stage A record is
neither an implementation nor a certification claim.

## Decision

LOSAT `tblastn -subject` applies the selected `-db_gencode` to subject
translation, candidate search, HSP re-evaluation, scoring, statistics,
coordinates, and `outfmt 0/6/7` reporting. The accepted IDs are exactly the
27 IDs in the pinned NCBI `gc.prt`: 1–6, 9–16, 21–33 (including 32).
The TBLASTN CLI accepts 32 and explicitly rejects invalid IDs without
substituting code 1. This is confined to TBLASTN and does not change TBLASTX.

The `-subject` exception is confined to effects of the selected non-default
**subject** code. NCBI remains authoritative for call order, candidate rules,
HSP construction, linking, pruning, scoring and statistical formulas,
coordinates, sorting, and output format. Code 1 local output must match NCBI
bytes. For other codes, NCBI `-db` supplies a translation/search/reporting
reference, but its database statistics and headers are a separate contract;
do not compare `-subject` and `-db` raw bytes as if they were equivalent.

## Evidence and source boundary

- `blast_args.cpp:997-1056` restricts CLI `-db_gencode` to 26 IDs, omitting
  32, while `gc.prt:340-347` defines 32. `blast_aux.cpp:588-613` implements
  `FindGeneticCode(id)` through `CGen_code_table::GetNcbieaa`.
- `blast_args.cpp:2538-2557` builds local subjects through
  `CObjMgr_QueryFactory`; `blast_setup_cxx.cpp:800-810` takes the genetic code
  from the subject source. `tblastn_app.cpp:242-265` passes the option's DB
  code to the formatter. The pinned CLI results are in
  `docs/evidence/tlosan_stage_a/cli_20260923_verified/`.
- On the code-4 fixture, NCBI `-subject -db_gencode 4` kept the code-1
  score of 231 bits and E-value 1.75e-85 while showing 100% identity;
  `-db` with code 4 produced 251 bits, E-value 3.39e-93, and an
  additional HSP. This is a local-subject code application boundary, not
  permission to alter scoring or statistics independently.

Code 32 is not certified until a comparison-only NCBI C++ API harness runs
`FindGeneticCode(32)` and checks translation plus the resulting TBLASTN HSP,
statistics, and format 0/6/7 output. That harness must not enter LOSAT's
runtime, build, or distributed artifacts. The Stage A evidence record tracks
its build and validation state.
