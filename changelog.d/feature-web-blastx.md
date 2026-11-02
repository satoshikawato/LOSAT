- BLASTX has the shared entry of LOSAT Web (web ABI v2), `run_local_blastx`: one search over the registered
  records written as outfmt 0, 6 and 7 in one call, the warnings once (the command line's standard error), the
  final hit list and the byte range of each HSP in each format. The query batches are searched by the command
  line's loop, so the bytes equal the command line's of each format. Other formats, a custom field list,
  `-db` and `-parse_deflines` are refused. The adapter and the app do not offer BLASTX yet.
