- `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, `rpsblast` and `rpstblastn` with
  `-parse_deflines`: a query batch with a repeated Seq-id in which every searched query fails query
  setup (its searched interval passes the end of the earlier record of its id or is empty, or a
  translated query is shorter than one codon) now ends with NCBI's exit-3 setup error
  (`CSeqVector_CI: cannot locate segment at ...` / `Sequence contains no data`) in every output
  format, as NCBI BLAST+ does; with `blastn -dust no -filtering_db`, in the form of NCBI's repeat
  search. In `blastn` with DUST (the default) or WindowMasker filtering, a batch in which any
  searched query passes the end of its earlier record now ends the same way, with that query's
  message once for each searched query. Not with `-export_search_strategy`; `psiblast` stays
  refused there because NCBI BLAST+ crashes.
- `blastp`, `tblastn`, `blastx` and `tblastx` with `-parse_deflines` and an XML, XML2 or JSON
  output format (`-outfmt 5`, `13`, `14`, `15`, `16`) whose `-matrix` name NCBI's XML writer does
  not know: a non-split batch whose queries have no valid Karlin-Altschul context (such as only `X`
  or `N` letters) now also ends with NCBI's `unsupported score matrix` error when the run's
  subjects pass setup, matching NCBI BLAST+ output and exit status.
