- README: "What LOSAT refuses" names every option the parser refuses by name, per program (`tblastn -in_pssm`
  and `-use_sw_tback`, `rpsblast`'s strategy options, the NCBI C++ Toolkit options `-xmlhelp`, `-logfile`,
  `-conffile`, `-dryrun`, `-help-full` and the `-version-full` forms, `tblastx -max_intron_length`) and some value
  refusals that were missing (`tblastn -evalue` of DBL_MAX, `-threshold` with the
  compressed lookup table, a hit list size that overflows NCBI's 32-bit size, `blastn -ungapped` with `-outfmt` 8, 9, 11,
  12 and 17), and "Programs LOSAT does not provide" lists all 28 executables of the NCBI BLAST+ 2.17.0 `bin` directory
  that LOSAT lacks. A unit test (`cli::readme_refusals_tests`) fails when a refusal list or an NCBI executable is not
  in the README again.
