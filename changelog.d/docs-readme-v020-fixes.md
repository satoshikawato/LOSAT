- docs: README's program table names `-outfmt 18` for the five search programs (the output-format table
  and the programs had it); the refusal list no longer says that `blastx` refuses `-query_loc` and
  `-subject_loc` with `-outfmt` 5, 8, 9 and 11 to 16 (it writes NCBI BLAST+ 2.17.0's bytes there) and
  names the refusal of a protein subject over 5,000,000 residues in `blastx`.
  `docs/losatx_blastx_native_cli.md` drops its old status (HOLD) and its list of 30 fields.
