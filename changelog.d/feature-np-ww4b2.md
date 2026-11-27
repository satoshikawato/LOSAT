- LOSAT Web: a makeblastdb tab, next to WindowMasker, builds a BLAST database in the browser with every
  option of `losat makeblastdb`, grouped as its help groups them: FASTA pasted or opened (`.gz`, `.bz2`
  and `.zst` too), ASN.1 files, or a BLAST database as `-in` (`-input_type`), `-dbtype`, `-out`,
  `-mask_data`, `-taxid_map`, `-gi_mask`, `-logfile` and the rest. The command beside Run LOSAT is the
  command line that runs, and a refused option is marked in its field. A finished run lists the files it
  made with their sizes; each one downloads, or all of them as one `.tar`. "Use as Database (-db)" opens
  the Search tab with the new database chosen (blastn for nucleotide, blastp for protein), and "Use as
  Repeats database (-filtering_db)" does the same for a nucleotide database. A WindowMasker run that
  writes mask data that makeblastdb reads (`-outfmt maskinfo_asn1_bin` or `maskinfo_asn1_text`) has "Use as
  -mask_data in makeblastdb", which adds its file to the makeblastdb form.
