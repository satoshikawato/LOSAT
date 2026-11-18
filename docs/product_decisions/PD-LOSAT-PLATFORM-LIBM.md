# Product Decision: a platform's libm and the last digit of a full-precision value

- Decision ID: `PD-LOSAT-PLATFORM-LIBM`
- Version: 1.1
- Date: 2026-10-25 (1.0, macOS: Owner, release board question 2, option (a); for v0.2.0);
  2026-11-18 (1.1, WASI: Owner, see Amendments)
- Status: Accepted for v0.2.0. A bit-identical `exp` / `log` port is a later line, not part of
  this decision.

## Scope

Two platforms where the libm that LOSAT calls for `exp` and `log` can differ by 1 ULP from glibc's
on an x86-64 CPU with FMA, where NCBI BLAST+ 2.17.0 wrote the frozen bytes of the tests (Linux x64):

- **macOS** (arm64 and x86_64), since 1.0: Apple's libm, when LOSAT computes an E value. The
  condition is an output format that prints the E value with all its digits. Of the formats LOSAT
  writes, these are the ASN.1 text formats: outfmt 11 (BLAST archive, `real { M, 10, E }` with 15
  digits, `LOSAT/src/report/serial/real.rs:882`) and outfmt 8 (Seqalign, same writer). The XML and
  JSON formats (5, 13, 14, 15, 16) round the E value to 6 significant digits first (`round6`,
  `LOSAT/src/report/blast_xml/hsp.rs:18`), and the tabular formats (6, 7, outfmt 0 tables) print 2
  or 3 significant digits; neither is affected in practice. This record names outfmt 11 and 8 only.
- **WASI** (`wasm32-wasip1` serial and `wasm32-wasip1-threads`), since 1.1: Rust's `f64::ln` there is
  wasi-libc's `log` (musl's, the same algorithm as glibc's without FMA; WebAssembly has no FMA
  instruction). The condition is psiblast's ASN.1 PSSM (`-out_pssm`) written with its diagnostics,
  which psiblast adds when `-out_ascii_pssm` is given too (`PSIDiagnosticsRequestNewEx`,
  `LOSAT/src/core/psi/mod.rs`): its `informationContent` values are sums of `log(qOverPEstimate)`
  terms (`blast_psi_priv.c:2360-2378`; LOSAT `save_diagnostics`, `LOSAT/src/core/psi/stages.rs`),
  printed as ASN.1 REAL values `{ M, 10, E }` with 15 digits. The ASCII PSSM prints the same values
  with 2 decimals and is not affected in the cases tested. This record names, for WASI, that file
  only: it is the one output of the release test set in which the difference was seen (Observation,
  "What was compared on WASI").

It is distinct from [`PD-NCBI-PLATFORM-VARIANCE`](PD-NCBI-PLATFORM-VARIANCE.md), which governs
NCBI's own platform variance in the certification gates. It changes no Gate A manifest and no
production code.

## Decision

On macOS, the 15-digit mantissa of a full-precision E value may differ from the Linux and
Windows value by 1 in its last digit, with the same exponent. On WASI, the 15-digit mantissa of a
value of psiblast's ASN.1 PSSM may differ from the Linux x64 value by 1 in its last digit, with the
same exponent. LOSAT does not work around either: no platform-specific code path, no normalization
of the output. The tests accept the difference on these platforms and nowhere else: the macOS
platform checks, and the WASI guard for the PSSM files.

## Observation

macOS (1.0):

- `blastx` case `lcbatch_sub` (`LOSAT/tests/outfmt11_blastx.rs`): LOSAT on macOS prints
  `real { 693693758698358, 10, -80 }`; NCBI's frozen bytes (Linux x64, as are LOSAT's on Linux
  x64, Linux arm64 and Windows) have `real { 693693758698357, 10, -80 }`.
- Cause: Apple's `exp` / `log` differ from the other platforms' by 1 ULP, which moves the last
  of the 15 printed digits.
- `platform-checks` runs `40171992341` and `40175184944` failed only in `cargo test
  (macos-x64)` and `cargo test (macos-arm64)`, in `archives_equal_ncbis_at_one_and_four_threads`
  and `the_flow_does_not_change_the_archives`, on this case.
- `platform-checks` run `41062695644` (dev) failed only in `cargo test (macos-x64)` and `cargo
  test (macos-arm64)`, in `stats::spouge::tests::pinned_double_rounding_erfc_and_spouge_calls_match_ncbi_double_bits`
  (`LOSAT/src/stats/spouge.rs`). The pinned ErfC rows passed (the port uses no libm); the pinned
  `blast_spouge_stoe` rows, whose E value multiplies three libm `exp` results, differ from the
  Linux bits by 2 ULP on arm64 and 4 ULP on x64 (first failing row each). Same accepted
  difference, no change of decision: on macOS that test accepts up to 8 ULP and is exact elsewhere.
- NCBI BLAST+ built on macOS has the same property: it calls the same libm (see
  `PD-NCBI-PLATFORM-VARIANCE` for NCBI's own platform variance).

WASI (1.1):

- In 14 psiblast `-db` rows of `LOSAT/tests/db_search_manifest/` (`db.ps.asn.tk.b.6`,
  `db.ps.inpssm.f.each.tk`, `db.ps.inpssm.f.last.tk`, `db.ps.it.f.both.tk`, `db.ps.it.f.last.tk`,
  `db.ps.it.f.eachlast.tk`, `db.ps.it.f.each.pk5`, `db.ps.it.f.each.three`, `db.ps.sp.f.each.tk`,
  `db.ps.q1.lc.tk.it0.6`, `db.ps.q2.b80.each`, `db.ps.q3.cull1.each`,
  `db.ps.q3.xall.inpssm.it2.each`, `db.ps.q3b.paja.each`), the serial and the threaded WASI
  commands write the report, standard error, exit status and ASCII PSSM of NCBI's frozen run, and
  an ASN.1 PSSM whose `informationContent` differs in the last of its 15 digits, for example
  `{ 347923230575625, 10, -15 }` where the frozen file and the native commands have
  `{ 347923230575626, 10, -15 }`. Every other byte of these files is the same. (Found by the
  psiblast WASI work of 2026-11-17, #1281.)
- NCBI's own output depends on the CPU: on x86-64, glibc picks an FMA variant of `log()` on a CPU
  with FMA. NCBI BLAST+ 2.17.0 on a Xeon 6338 (glibc 2.35) writes the frozen bytes by default; with
  glibc's FMA variants off (`GLIBC_TUNABLES=glibc.cpu.hwcaps=-AVX2,-FMA,-FMA4`) it writes the WASI
  bytes: on each of the 14 rows its report, standard error, exit status and every PSSM file (60
  files) are the bytes of the serial WASI command (2026-11-18, the command built from dev
  `e87bb4538`).
- Rust's `f64::ln` on `wasm32-wasip1` differs from glibc 2.39's `log` by 1 ULP in 7 of 200,000
  random inputs, and in none with glibc's FMA variant off.

What was compared on WASI (1.1), to bound the scope: the serial command of dev `e87bb4538` ran
every row of the NCBI fixture sets of the release test plan, D (`outfmt0_manifest.tsv`, 2,602
`-subject` rows) and E (`db_search_manifest/`, 6,571 `-db` rows), at one thread, against the frozen
bytes (2026-11-18). Outside psiblast's ASN.1 PSSM no value differs in its digits: every `-outfmt 11`
and `-outfmt 8` row of the five search programs gives its expected bytes (483 `-subject` rows; 563
`-db` rows, 18 of them LOSAT's refusals), and so does every psiblast `-subject` row (95, the
`-out_pssm` row among them). The 14 rows above are the only ones of this record's kind: their 24
differing PSSM files hold 89 values that are 1 apart in the last digit and no other differing byte,
and the threaded command gives the same 14 rows and bytes (set G, one thread). The same run found
two differences of other causes, which this record does not cover: blastx with an XML or JSON format
(`-outfmt 5`, 13 to 16) stops before the search on WASI (88 `-subject` and 210 `-db` rows), and 8
`-db` rows whose alias databases total 2^32 letters or more are refused on wasm32. Those rows show
nothing about their values; the XML and JSON formats round the E value to 6 digits (the macOS scope
above).

## Verification

One rule, in two copies: `archives_match_allowing_libm_last_digit` in
`LOSAT/tests/outfmt11_blastx.rs` (the `-outfmt 11` blastx tests, macOS) and
`libm_last_digit_differences` in `LOSAT/tests/xinclude_fixtures.py` (the Python tools, WASI). A
difference is accepted only when every byte outside the ASN.1 REAL values is identical and each
differing value has the same exponent and a mantissa of 15 digits, as the expected one, that
differs from it by exactly 1. A mantissa 2 away, another exponent, a mantissa of another or fewer
digits (NCBI drops trailing zeros, so a shorter mantissa 1 away is not a last-digit change), or any
other byte fails. The Rust copy reads the values of a Seq-align (`real { M, 10, E }`); the Python
copy reads every REAL value `{ M, 10, E }` of ASN.1 text, which in an archive is the same token and
in an ASN.1 PSSM is how the values are written (`lambda { M, 10, E }`, the values of
`informationContent { ... }`). Their unit tests run on every OS (`cargo test`, and
`LOSAT/tests/test_ci_fast_regressions.py` in CI).

- macOS: the `-outfmt 11` blastx tests apply the rule with `cfg!(target_os = "macos")`; on every
  other OS the comparison is exact.
- WASI: `LOSAT/tests/check_wasm_threading.py` (set G of the release test plan, `wasm-threading.yml`)
  runs every psiblast `-db` row (`db_search_fixtures.py check --program psiblast --threads 1`)
  through each WASI command it is given, with `--libm-last-digit` naming its native command. The
  rule applies to the PSSM files only; the report, standard error and exit status are compared byte
  for byte. A PSSM file over 64 KB is not frozen (it is compared by its SHA-256); its expected bytes
  are those of the native command's run of the row, which must have the frozen SHA-256. A row taken
  by the rule is listed as `libm-last-digit` with its files. Native runs of the same checker (sets
  D, E, F) pass no option and stay exact.

## Release-notes sentence

On macOS, an E value printed with full precision (`-outfmt 11`, `-outfmt 8`) can differ from the
other platforms' output by 1 in its last (15th) digit because Apple's libm differs from the
other platforms' by 1 ULP; NCBI BLAST+ built on macOS has the same property, and the formats
that print 2 to 6 digits are not affected in the cases tested. On WASI, a value of psiblast's
ASN.1 PSSM (`-out_pssm` with `-out_ascii_pssm`) can differ in the same way, because WebAssembly's
`log()` has no FMA; NCBI BLAST+ writes the same digit on an x86-64 CPU without FMA.

## Amendments

- 1.1 (2026-11-18, Owner): WASI is added, with the same rule. The Owner answered 「WASI=b」 on
  2026-11-18 at 12:2x JST, in session losat-baselines-04, to this question (as shown to the Owner):
  on WASI, psiblast with `-db` writing a PSSM gives, in 14 of 954 rows, a last digit of
  `informationContent` that differs from Linux (glibc's FMA variant of `log()` against WASI's
  `log()` without FMA; NCBI itself gives the WASI value on a CPU without FMA); (a) port the FMA
  variant of `log` to every platform / (b) extend the exception record of 10-25 for macOS,
  PD-LOSAT-PLATFORM-LIBM, to WASI (recommended) / (c) README only. Not taken: (a) and (c). The
  guard (`check_wasm_threading.py`) runs the psiblast `-db` rows under the rule; the rule's two
  copies now also require 15-digit mantissas, as the Decision says (1.0's test compared digit
  counts only, which let a shorter mantissa 1 away pass).
