- `psiblast` through the library (`losat::api`) no longer takes the `delim=` of `-outfmt 6`, `7`, `10` and `20`: its
  reports keep the comma or the tab, as NCBI BLAST+ 2.17.0's `psiblast` and LOSAT's command line write them, and
  `rpsblast -outfmt "6 delim=X"` (or `7`) counts the identities of masked queries as its default report does. One rule
  for every program and both paths: blastn, blastp, blastx, tblastn and tblastx take `delim=`; psiblast, rpsblast and
  rpstblastn check it and keep their own delimiter (line RG, series F, D-36).
- `blastp` and `psiblast` through the library accept a format number written as NCBI reads it (`-outfmt 06`,
  `"06 std"`, `+7`, `00`), as the command line already did; `rpsblast` and `rpstblastn -outfmt "20 delim=; qseqid
  staxids"` stop with NCBI's error for a delimiter that the fields use, not with LOSAT's refusal of outfmt 20 (line RG,
  series F, D-36 review 1).
- Through the library, `blastp` runs `-outfmt 0` to `4` and `psiblast` runs `-outfmt 0` with words after the number
  (`"0 std"`, `"0 delim=@"`, `"00 qseqid"`, `"1 qseqid"`) and drop the words, as NCBI BLAST+ 2.17.0 and LOSAT's command
  line do, instead of refusing them (`psiblast`'s `1` to `4` stay refused, with or without words). `psiblast`,
  `rpsblast` and PSI-TBLASTN (`tblastn -in_pssm`) with `-outfmt 13` or `14` and no `-out` file now stop with NCBI's
  `Please provide a file name for outfmt 14.`, in NCBI's order, instead of LOSAT's refusal of the format (line RG,
  series F, D-36 review 1b).
- With `-outfmt 13` or `14` and no `-out` file, or another `-outfmt` value that NCBI BLAST+ 2.17.0 stops at (`17`
  outside `blastn`, `19`, `21`, a delimiter that the fields use, a number it cannot read), every search program stops
  with NCBI's error even when the run also sends its search strategy (`-export_search_strategy -` or `/dev/stdout`) or
  `psiblast`'s PSSM (`-out_pssm`, `-out_ascii_pssm`) to the standard output, instead of LOSAT's refusal of that file:
  NCBI never writes them. `blastx` through the library gives NCBI's `Please provide a file name for outfmt 13.` for
  13 and 14 without a receiver of their side files, as its command line and the other programs' library do, instead
  of LOSAT's own refusal. With an `-outfmt` value that NCBI cannot read (`99`), `psiblast`, `rpsblast` and
  `rpstblastn` create the `-export_search_strategy` file before NCBI's error, as NCBI does (line RG, series F, D-36
  review 2).
- Every search program's refusal of a search strategy on a report stream (`-export_search_strategy -`, or the file of
  `-out` or of the standard output) now comes where NCBI BLAST+ 2.17.0 would start its report. A run that NCBI ends
  before that point now ends as NCBI's: an empty or white-space query gives `Query is Empty!` with exit status 0, a
  missing database NCBI's database error, an option error NCBI's error (before: LOSAT's refusal, exit status 1). With
  a strategy file and an `-outfmt` value that NCBI cannot read (`99`), every program also creates the `-out` file
  before NCBI's error, as NCBI does (line RG, series F, D-36 review 2b).
