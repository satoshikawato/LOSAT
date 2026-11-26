- `losat makeblastdb`: added `-input_type asn1_txt` and `-input_type asn1_bin` support for ASN.1
  text and binary `Seq-entry` and `Bioseq-set` inputs (including `.gz`, `.bz2`, `.zst`, and
  `stdin`), matching NCBI BLAST+ 2.17.0 byte-for-byte across `ASN1_BlastDefLine` binary headers,
  `CFastaReader` user objects, `s_GetTaxId`, nucleotide delta sequences, `-mask_data`,
  `-hash_index`, `-oid_masks`, `-logfile`, `NCBI_MAKEBLASTDB_SKIP_UNVERIFIED_BIOSEQS`, and
  `NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS`, within the explicit,
  source-taxname, and default protein title scope (explicitly refusing title paths planned for
  follow-up PRs and NCBI buffer-overread / null-dereference inputs).
