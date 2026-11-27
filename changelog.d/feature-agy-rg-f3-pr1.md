- `losat blast_formatter`: writes the reports of a search again from its BLAST archive (the file of
  `-outfmt 11` from LOSAT or NCBI BLAST+, a file or the standard input with `-archive -`) without
  running the search again, as NCBI BLAST+ 2.17.0 `blast_formatter` writes them: outfmt `8`, `9` and
  `12` (the Seq-aligns in ASN.1 text, binary ASN.1 and JSON) for each query batch of the archive, for
  archives of every search program with `-subject` or `-db`. Its command line is NCBI's (`-archive`,
  `-outfmt`, `-show_gis`, `-num_descriptions`, `-num_alignments`, `-line_length`, `-html`, `-sorthits`,
  `-sorthsps`, `-max_target_seqs`, `-out`, `-parse_deflines`, `-h`, `-help`, `-help-full`, `-version`,
  `-xmlhelp`, `-dryrun`, `-logfile` and `-conffile`), with NCBI's help, argument errors, messages and
  exit statuses. It reads archives in ASN.1 text, the default of `-outfmt 11`; it refuses `-outfmt`
  0 to 7, 10, 11, 13 to 18 and 20, an archive in ASN.1 binary or XML, a database of the archive that is
  not found as NCBI finds a database, `-rid`, `PRE_FETCH_SEQS_LIMIT=2147483647`,
  archives with a stored message of severity Fatal and archives on which NCBI BLAST+ has no defined
  result. The library API runs it as `run_blast_formatter`
  (pylosat: `blast_formatter`), and `losat --install-shims` links `blast_formatter` too.
