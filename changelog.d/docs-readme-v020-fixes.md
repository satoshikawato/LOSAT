- docs: README's program table names `-outfmt 18` for the five search programs (the output-format table
  and the programs had it); the refusal list no longer says that `blastx` refuses `-query_loc` and
  `-subject_loc` with `-outfmt` 5, 8, 9 and 11 to 16 (it writes NCBI BLAST+ 2.17.0's bytes there) and
  names the refusal of a protein subject over 5,000,000 residues in `blastx`.
  `docs/losatx_blastx_native_cli.md` drops its old status (HOLD) and its list of 30 fields.
- docs: `docs/release/v0.1.0.md` says v0.1.0 was published (GitHub Release, 2026-09-17); the "Configuration
  boundary" of `docs/cli_v2_migration.md` says that stdin and `-db` were missing at v0.1.0 and work since
  v0.2.0; README says that `--install-shims` replaces a link or copy an earlier LOSAT made.
- docs: README's output-format table says "all five programs" for outfmt 11 (it said every program) and
  points to the program table for the formats of `psiblast`, `rpsblast` and `rpstblastn`; the
  `convert2blastmask` paragraph no longer says that `losat makeblastdb` lacks `-mask_data`.
- docs: README says that the Web GUI follows the page layout of NCBI BLAST (it said "the layout and terms";
  its words are LOSAT's own).
