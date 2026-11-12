- `losat dustmasker`: what NCBI BLAST+ 2.17.0 `dustmasker` 1.0.0 writes, byte for byte. It finds the
  low-complexity intervals of the sequences of a FASTA file (`-in`, the standard input by default) or
  of a BLAST database (`-infmt blastdb`) and writes them in the nine `-outfmt` forms (`interval`,
  `acclist`, `fasta` with `-hard_masking`, `seqloc_asn1_bin`/`_text`/`_xml`,
  `maskinfo_asn1_bin`/`_text`/`_xml`), with `-parse_seqids`, `-window`, `-level`, `-linker`, `-out`,
  `-logfile`, `-conffile`, and NCBI's error messages and exit codes. NCBI's `makeblastdb -mask_data` reads its
  `maskinfo_*` files. `losat --install-shims` links `dustmasker` too. blastn's query masking uses the
  same DUST code; its output does not change and it runs no slower.
