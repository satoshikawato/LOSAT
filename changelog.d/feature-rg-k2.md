- `losat segmasker`: what NCBI BLAST+ 2.17.0 `segmasker` 1.0.0 writes, byte for byte. It finds the
  low-complexity intervals of the protein sequences of a FASTA file (`-in`, the standard input by
  default) or of a protein BLAST database (`-infmt blastdb`) with SEG and writes them in the eight
  `-outfmt` forms (`interval`, `fasta`, `seqloc_asn1_bin`/`_text`/`_xml`,
  `maskinfo_asn1_bin`/`_text`/`_xml`), with `-parse_seqids`, `-window`, `-locut`, `-hicut`, `-out`,
  `-logfile`, `-conffile`, `-dryrun`, and NCBI's error messages and exit codes. NCBI's
  `makeblastdb -mask_data` reads its `maskinfo_*` files. The library API runs it as `run_segmasker`
  (pylosat: `segmasker`), and `losat --install-shims` links `segmasker` too. The SEG of the
  searches' masking is the same code; their output does not change.
- `losat blastdbcmd` writes the byte 255 for a protein residue code above 27 (only in a hand-made
  or damaged `.psq`) in its FASTA and `%s` letters, as NCBI does; it wrote a letter of the code's low
  bits before.
