- `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, `rpsblast` and `rpstblastn` with `-parse_deflines`: when every
  searched query of a batch with a repeated Seq-id fails query setup (its searched interval passes
  the earlier record's end or is empty, alongside any translated query shorter than one codon),
  LOSAT now ends with NCBI's exit-3 setup error (`CSeqVector_CI: cannot locate segment at ...` /
  `Sequence contains no data`) in every output format (except `-export_search_strategy`, and
  `psiblast` stays refused where NCBI BLAST+ crashes).
- `blastp`, `tblastn`, `blastx` and `tblastx` with `-parse_deflines` and an XML, XML2 or JSON
  output format (`-outfmt 5`, `13`, `14`, `15`, `16`) whose `-matrix` name NCBI's XML writer does
  not know: a non-split batch whose queries have no valid Karlin-Altschul context (such as only `X`
  or `N` letters) now also ends with NCBI's `unsupported score matrix` error when the run's
  subjects pass setup, matching NCBI BLAST+ output and exit status.
