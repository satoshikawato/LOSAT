# E2g angle C findings (running)

(working notes; final text below at end)
- opt coverage: NCBI 2.17.0 `-help` lists 71 option names; LOSAT implements 18 + -help and lists 52 as unported (cli.rs:205-262). comm shows no NCBI option is neither implemented nor listed. All 52 x 3 syntaxes exit 2 with a message naming LOSAT (opt_rejection_test.tsv).
- NCBI CNcbiApplication standard flags -version-full, -xmlhelp, -conffile, -logfile, -dryrun (accepted by NCBI, absent from -help) are not in INVENTORY (grep) and LOSAT gives "unknown option or argument" exit 2 (explicit, but not a mapped rejection).
- blastn-short is cheap: NCBI -task blastn-short == LOSAT -task blastn -reward 1 -penalty -3 -evalue 1000 -word_size 7 -dust no, byte-identical outfmt 0/6/7 on a 5-query set.
