- PSIBLAST and TBLASTN with `-in_pssm`, `-db` and `-export_search_strategy`: a PSSM file in ASN.1
  binary or XML whose author names hold a number too large for its type is now refused after the
  start of the outfmt 0 to 4 report (the program line, the references and the `Database:` lines),
  where BLAST+ starts to write the strategy file and stops with its exception, as `-out_pssm` with
  such a file already was. Before, the refusal came before any output.
