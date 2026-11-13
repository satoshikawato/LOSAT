- `losat convert2blastmask`: what NCBI BLAST+ 2.17.0 `convert2blastmask` writes, byte for byte. It
  turns the lower case letters of FASTA input (`-in`, the standard input by default), and its lines
  `from - to`, into the masking information that NCBI's `makeblastdb -mask_data` reads, in the four
  `-outfmt` forms (`maskinfo_asn1_text`, the default, `maskinfo_asn1_bin`, `maskinfo_xml` and
  `interval`), with `-parse_seqids`, `-masking_algorithm`, `-masking_options`, `-out`, and NCBI's
  error messages and exit codes. The library API runs it as `run_convert2blastmask` (pylosat:
  `convert2blastmask`), and `losat --install-shims` links `convert2blastmask` too.
