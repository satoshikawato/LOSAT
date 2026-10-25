# Product Decision: Apple libm and the last digit of a full-precision E value

- Decision ID: `PD-LOSAT-PLATFORM-LIBM`
- Version: 1.0
- Date: 2026-10-25 (Owner, release board question 2, option (a); for v0.2.0)
- Status: Accepted for v0.2.0. A bit-identical `exp` / `log` port is a later line, not part of
  this decision.

## Scope

The platform is macOS (arm64 and x86_64), where LOSAT calls Apple's libm for `exp` and `log`
when it computes an E value. The condition is an output format that prints the E value with all
its digits. Of the formats LOSAT writes, these are the ASN.1 text formats: outfmt 11 (BLAST
archive, `real { M, 10, E }` with 15 digits, `LOSAT/src/report/serial/real.rs:882`) and
outfmt 8 (Seqalign, same writer). The XML and JSON formats (5, 13, 14, 15, 16) round the E value
to 6 significant digits first (`round6`, `LOSAT/src/report/blast_xml/hsp.rs:18`), and the
tabular formats (6, 7, outfmt 0 tables) print 2 or 3 significant digits; neither is affected in
practice. This record names outfmt 11 and 8 only.

It is distinct from [`PD-NCBI-PLATFORM-VARIANCE`](PD-NCBI-PLATFORM-VARIANCE.md), which governs
NCBI's own platform variance in the certification gates. It changes no Gate A manifest and no
production code.

## Decision

On macOS, the 15-digit mantissa of a full-precision E value may differ from the Linux and
Windows value by 1 in its last digit, with the same exponent. LOSAT does not work around it: no
platform-specific code path, no normalization of the output. The platform-checks tests accept
the difference on macOS and nowhere else.

## Observation

- `blastx` case `lcbatch_sub` (`LOSAT/tests/outfmt11_blastx.rs`): LOSAT on macOS prints
  `real { 693693758698358, 10, -80 }`; NCBI's frozen bytes (Linux x64, as are LOSAT's on Linux
  x64, Linux arm64 and Windows) have `real { 693693758698357, 10, -80 }`.
- Cause: Apple's `exp` / `log` differ from the other platforms' by 1 ULP, which moves the last
  of the 15 printed digits.
- `platform-checks` runs `40171992341` and `40175184944` failed only in `cargo test
  (macos-x64)` and `cargo test (macos-arm64)`, in `archives_equal_ncbis_at_one_and_four_threads`
  and `the_flow_does_not_change_the_archives`, on this case.
- NCBI BLAST+ built on macOS has the same property: it calls the same libm (see
  `PD-NCBI-PLATFORM-VARIANCE` for NCBI's own platform variance).

## Verification

`archives_match_allowing_libm_last_digit` in `LOSAT/tests/outfmt11_blastx.rs` is the one
comparison of the `-outfmt 11` blastx tests. With the macOS rule on (`cfg!(target_os =
"macos")`), a difference is accepted only when every byte outside `real { M, 10, E }` tokens is
identical and each differing token has the same exponent, the same digit count and mantissas
that differ by exactly 1. A mantissa 2 away, another exponent, a changed digit count, or any
other byte fails. On every other OS the comparison is exact. Its unit test runs on every OS.

## Release-notes sentence

On macOS, an E value printed with full precision (`-outfmt 11`, `-outfmt 8`) can differ from the
other platforms' output by 1 in its last (15th) digit because Apple's libm differs from the
other platforms' by 1 ULP; NCBI BLAST+ built on macOS has the same property, and the formats
that print 2 to 6 digits are not affected in the cases tested.
