- `psiblast` through the library (`losat::api`) no longer takes the `delim=` of `-outfmt 6`, `7`, `10` and `20`: its
  reports keep the comma or the tab, as NCBI BLAST+ 2.17.0's `psiblast` and LOSAT's command line write them, and
  `rpsblast -outfmt "6 delim=X"` (or `7`) counts the identities of masked queries as its default report does. One rule
  for every program and both paths: blastn, blastp, blastx, tblastn and tblastx take `delim=`; psiblast, rpsblast and
  rpstblastn check it and keep their own delimiter (line RG, series F, D-36).
- `blastp` and `psiblast` through the library accept a format number written as NCBI reads it (`-outfmt 06`,
  `"06 std"`, `+7`, `00`), as the command line already did; `rpsblast` and `rpstblastn -outfmt "20 delim=; qseqid
  staxids"` stop with NCBI's error for a delimiter that the fields use, not with LOSAT's refusal of outfmt 20 (line RG,
  series F, D-36 review 1).
- `blastp` and `psiblast` through the library run `-outfmt 0` to `4` with words after the number (`"0 std"`,
  `"0 delim=@"`, `"00 qseqid"`, `"1 qseqid"`) and drop the words, as NCBI BLAST+ 2.17.0 and LOSAT's command line do,
  instead of refusing them. `psiblast`, `rpsblast` and PSI-TBLASTN (`tblastn -in_pssm`) with `-outfmt 13` or `14` and
  no `-out` file now stop with NCBI's `Please provide a file name for outfmt 14.`, in NCBI's order, instead of LOSAT's
  refusal of the format (line RG, series F, D-36 review 1b).
