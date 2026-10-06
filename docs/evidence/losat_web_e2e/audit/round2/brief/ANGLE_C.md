# Angle (c): TBLASTX — options and application flow

Read `COMMON.md` first. Work dir: `/home/kawato/.cache/losat-web-gui-target/s08pb-audit/c/`. Program: `tblastx` (NCBI) vs `$FINAL tblastx`.

1. Re-run round 1: every reproducing command and harness case in `src/docs/evidence/losat_web_e2e/audit/round1/tblastx.md` (TX-1..TX-10 and the ~1000 argv) — S08+a's copy is `/home/kawato/.cache/losat-web-gui-target/s08pa/audit-rerun/tblastx/` (`cmp.sh` and argv lists); copy it into `c/` and point LOSAT at `$FINAL`. Give a verdict per round 1 finding and totals per class. Skip NCBI runs that do not end (TX-2 window values); they are rejected by D8.
2. New:
   - `-threshold` extremes (`+inf`, `1e300`, `2147483648`, `0.5`, `1`) with `-word_size 2/3/4` on small inputs: in an earlier aborted gate run LOSAT ended with signal 9 on `-threshold +inf` once; measure LOSAT's time and max RSS (`/usr/bin/time -v`, `timeout 300`) and compare output with NCBI. Report any run where LOSAT uses far more memory or time than NCBI.
   - `-db_gencode` 2, 5, 6, 9, 11, 12 with `-query_gencode` variants and `-evalue` values (approved exception: compare with NCBI `-db` built by `makeblastdb -dbtype nucl`, outside database lines; see COMMON.md).
   - `-culling_limit` (ported in S08+, `LOSAT/src/algorithm/tblastx/hsp_culling.rs`): values 1, 2, 3, 10, 0x2, with `-max_target_seqs` 1..5 and `-evalue` 10/1000, on multi-subject inputs.
   - NCBI toolkit words in an option's value position (`-out -version`, `-out -dryrun`, `-outfmt -version-full`, ...): explicit rejection (pending the maintainer).
   - The CLI changes of S08+a in `LOSAT/src/cli.rs` (`ncbi_preparsed_toolkit_word`, `unknown_option_error`) must not change any accepted TBLASTX argv: re-run the accepted argv lists of round 1.
3. Report per COMMON.md.
