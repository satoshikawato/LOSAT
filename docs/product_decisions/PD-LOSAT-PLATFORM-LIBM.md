# Product Decision: a platform's libm and the last digit of a full-precision value

- Decision ID: `PD-LOSAT-PLATFORM-LIBM`
- Version: 1.2 (2026-11-20, Owner: the WASI clause of 1.1 is retired)
- Date: 2026-10-25 (1.0, macOS: Owner, release board question 2, option (a); for v0.2.0);
  2026-11-18 (1.1, WASI: Owner, see Amendments); 2026-11-20 (1.2, WASI retired: Owner, see
  Amendments)
- Status: macOS clause (1.0): Accepted for v0.2.0; it stays, and is unverified on macOS (no Apple
  device). WASI clause (1.1): retired by 1.2 (Owner, 2026-11-20) because `exp`, `log`, `log10` and
  `pow` give the bits of glibc 2.35's FMA versions on one Rust path on every target (#1344: the code
  of Arm Optimized Routines, MIT, and Sun's fdlibm for `log10`).

## Scope

Two platforms where the libm that LOSAT calls for `exp` and `log` can differ by 1 ULP from glibc's
on an x86-64 CPU with FMA, where NCBI BLAST+ 2.17.0 wrote the frozen bytes of the tests (Linux x64).
Since 1.2 only macOS is in force; the WASI paragraph is kept as it stood in 1.1:

- **macOS** (arm64 and x86_64), since 1.0: Apple's libm, when LOSAT computes an E value. The
  condition is an output format that prints the E value with all its digits. Of the formats LOSAT
  writes, these are the ASN.1 text formats: outfmt 11 (BLAST archive, `real { M, 10, E }` with 15
  digits, `LOSAT/src/report/serial/real.rs:882`) and outfmt 8 (Seqalign, same writer). The XML and
  JSON formats (5, 13, 14, 15, 16) round the E value to 6 significant digits first (`round6`,
  `LOSAT/src/report/blast_xml/hsp.rs:18`), and the tabular formats (6, 7, outfmt 0 tables) print 2
  or 3 significant digits; neither is affected in practice. This record names outfmt 11 and 8 only.
- **WASI** (`wasm32-wasip1` serial and `wasm32-wasip1-threads`), since 1.1, no longer in force
  under 1.2 (see Amendments): Rust's `f64::ln` there is
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
Windows value by 1 in its last digit, with the same exponent. LOSAT does not work around it: no
platform-specific code path, no normalization of the output. The tests accept the difference on
macOS and nowhere else: the macOS platform checks.

Under 1.1, the same held on WASI for the 15-digit mantissa of a value of psiblast's ASN.1 PSSM.
Under 1.2 WASI is no longer in force: the port gives the Linux x64 bytes there (Amendments,
Verification), and the WASI guard compares the PSSM files byte for byte.

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

WASI (1.1; describes the builds before the port of #1344):

- In 14 psiblast `-db` rows of `LOSAT/tests/db_search_manifest/` (`db.ps.asn.tk.b.6`,
  `db.ps.inpssm.f.each.tk`, `db.ps.inpssm.f.last.tk`, `db.ps.it.f.both.tk`, `db.ps.it.f.last.tk`,
  `db.ps.it.f.eachlast.tk`, `db.ps.it.f.each.pk5`, `db.ps.it.f.each.three`, `db.ps.sp.f.each.tk`,
  `db.ps.q1.lc.tk.it0.6`, `db.ps.q2.b80.each`, `db.ps.q3.cull1.each`,
  `db.ps.q3.xall.inpssm.it2.each`, `db.ps.q3b.paja.each`), the serial and the threaded WASI
  commands write the report, standard error, exit status and ASCII PSSM of NCBI's frozen run, and
  an ASN.1 PSSM whose `informationContent` differs in the last of its 15 digits, for example
  `{ 347923230575625, 10, -15 }` where the frozen file and the native commands have
  `{ 347923230575626, 10, -15 }`. Every other byte of these files is the same. (Found by the
  psiblast WASI work of 2026-11-17, #1281.) Correction (1.2, the v0.2.0 release line, 2026-11-20): in these
  14 rows the differing value comes from `pow` in the pseudocounts (relative entropy to the power
  K, `LOSAT/src/core/psi/stages.rs`), not from `log` as this record said in 1.1.
- NCBI's own output depends on the CPU: on x86-64, glibc picks an FMA variant of `log()` on a CPU
  with FMA. NCBI BLAST+ 2.17.0 on a Xeon 6338 (glibc 2.35) writes the frozen bytes by default; with
  glibc's FMA variants off (`GLIBC_TUNABLES=glibc.cpu.hwcaps=-AVX2,-FMA,-FMA4`) it writes the WASI
  bytes: on each of the 14 rows its report, standard error, exit status and every PSSM file (60
  files) are the bytes of the serial WASI command (2026-11-18, the command built from dev
  `e87bb4538`).
- Rust's `f64::ln` on `wasm32-wasip1` differs from glibc 2.39's `log` by 1 ULP in 7 of 200,000
  random inputs, and in none with glibc's FMA variant off.

What was compared on WASI (1.1, before the port), to bound the scope: the serial command of dev
`e87bb4538` ran every row of the NCBI fixture sets of the release test plan, D (`outfmt0_manifest.tsv`, 2,602
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

The rule for macOS: `archives_match_allowing_libm_last_digit` in `LOSAT/tests/outfmt11_blastx.rs`
(the `-outfmt 11` blastx tests). A difference is accepted only when every byte outside the ASN.1
REAL values is identical and each differing value has the same exponent and a mantissa of 15
digits, as the expected one, that differs from it by exactly 1. A mantissa 2 away, another exponent,
a mantissa of another or fewer digits (NCBI drops trailing zeros, so a shorter mantissa 1 away is
not a last-digit change), or any other byte fails. The copy reads the values of a Seq-align
(`real { M, 10, E }`). Its unit tests run on every OS (`cargo test`).

- macOS: the `-outfmt 11` blastx tests apply the rule with `cfg!(target_os = "macos")`; on every
  other OS the comparison is exact. This copy stays under 1.2.
- WASI, since 1.2: `LOSAT/tests/check_wasm_threading.py` (set G of the release test plan,
  `wasm-threading.yml`) runs every psiblast `-db` row through each WASI command it is given and
  compares the PSSM files byte for byte, as the report, standard error and exit status (1,030 rows
  equal on the serial WASI command of the port's branch, 2026-11-20). The Python copy of the rule
  (`libm_last_digit_differences` in `LOSAT/tests/xinclude_fixtures.py`, and `db_search_fixtures.py
  check --libm-last-digit`) is removed. A PSSM file over 64 KB is not frozen (it is compared by its
  SHA-256); its expected bytes are those of the native command's run of the row, which must have the
  frozen SHA-256. Native runs of the same checker (sets D, E, F) stay exact.
- WASI, 1.1 (before the port): the same guard applied the rule to the PSSM files with
  `--libm-last-digit` naming its native command; a row taken by the rule was listed as
  `libm-last-digit` with its files.

### macOS

The port of #1344 runs on macOS too (`f64::mul_add` where the CPU has an FMA instruction, the
exact emulation elsewhere). It is unverified on macOS (no Apple device): nobody has checked it on a
Mac. The macOS clause (1.0) stays until the Owner retires it.

## Release-notes sentence

On macOS, an E value printed with full precision (`-outfmt 11`, `-outfmt 8`) can differ from the
other platforms' output by 1 in its last (15th) digit because Apple's libm differs from the
other platforms' by 1 ULP; NCBI BLAST+ built on macOS has the same property, and the formats
that print 2 to 6 digits are not affected in the cases tested. WASI is no longer named (1.2): its
psiblast ASN.1 PSSM files give the Linux x64 bytes.

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
- 1.2 (2026-11-20, Owner; proposed the same day by line LIBM): the WASI clause of 1.1 is
  retired; the macOS clause (1.0) stays and is unverified on macOS (no Apple device). Why: the branch
  `fix/libm-glibc-fma` (#1344) gives `exp`, `log`, `log10` and `pow` the bits of glibc 2.35's
  x86-64 FMA versions on one Rust path on every target (`LOSAT/src/utils/glibc_fma_libm.rs`, the
  code of Arm Optimized Routines, MIT, and Sun's fdlibm for `log10`): `f64::mul_add`
  where there is an FMA instruction (aarch64; x86-64 CPUs that have FMA, detected at run time), an
  exact emulation of the FMA elsewhere (WebAssembly, x86-64 without FMA). Every call site uses it,
  and clippy refuses the libm methods of `f64` elsewhere. Measured 2026-11-20 on the serial WASI
  command of that branch: all 1,030 psiblast `-db` rows of `LOSAT/tests/db_search_manifest/` give the
  frozen bytes byte for byte, with no rule (the 14 rows of 1.1 among them); every `-outfmt 8` and
  `-outfmt 11` row of D (483) and E (587) is equal; the 20 WQA web-page cases give the native bytes
  (dev's WASI command differed on 4 of them, the four tblastx cases, in E values printed with
  all 15 digits). On that
  branch the WASI guard compares the psiblast PSSM files byte for byte and the Python copy of the
  rule is removed; the macOS Rust copy stays. Correction to 1.1: in the 14 psiblast rows the
  differing value came from `pow` in the pseudocounts, not from `log` (the v0.2.0 release line, 2026-11-20).
  The Owner answered with no words: on 2026-11-20 at 05:39 JST the Owner pressed "Ready for review"
  on pull request #1345, which its description named as the way to accept option A of this question
  (Q-LIBM-PD, as shown to the Owner): retire the WASI clause of PD-LOSAT-PLATFORM-LIBM (1.1, Owner
  2026-11-18)? (A) retire it (recommended) / (B) keep 1.1; the macOS clause stays in both. Not
  taken: (B). #1345 merged into dev the same minute.
