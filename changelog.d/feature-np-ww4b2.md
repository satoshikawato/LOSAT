- LOSAT Web: a makeblastdb tab, next to WindowMasker, builds a BLAST database in the browser with every
  option of `losat makeblastdb`, grouped as its help groups them: FASTA pasted or opened (`.gz`, `.bz2`
  and `.zst` too), ASN.1 files, or a BLAST database as `-in` (`-input_type`), `-dbtype`, `-out`,
  `-mask_data`, `-taxid_map`, `-gi_mask`, `-logfile` and the rest. The command beside Run LOSAT is the
  command line that runs, and a refused option is marked in its field. A finished run lists the files it
  made with their sizes; each one downloads, or all of them as one `.tar`. "Use as Database (-db)" opens
  the Search tab with the new database chosen (blastn for nucleotide, blastp for protein), and "Use as
  Repeats database (-filtering_db)" does the same for a nucleotide database. A WindowMasker run that
  writes mask data (`-outfmt maskinfo_asn1_bin`, `maskinfo_asn1_text` or `maskinfo_xml`) has "Use as
  -mask_data in makeblastdb", which adds its file to the makeblastdb form.
  A makeblastdb run is marked Supported when its options match a run that LOSAT's tests compare byte for
  byte with NCBI BLAST+ 2.17.0: masks from `-mask_data` (with or without `-parse_seqids`, `-title`,
  `-mask_id` and `-mask_desc`), ASN.1 and BLAST database inputs, `-gi_mask` across several volumes,
  `-oid_masks`, `-randomize`, `-metadata_output_prefix`, `-hash_index`, `-taxid` and `-taxid_map`;
  compressed FASTA (`.gz`, `.bz2`, `.zst`) is among the compared inputs.
