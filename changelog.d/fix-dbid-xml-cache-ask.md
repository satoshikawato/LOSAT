- `-db` searches in the XML2 and JSON formats (outfmt 13 to 16; blastn, blastp, blastx, tblastn and
  tblastx) now refuse a run for a taxonomy database whose damaged last entry shares NCBI's taxonomy
  cache slot with another taxid of a database volume, where a reported subject has a taxid, as the
  other formats do where NCBI fetches a subject. NCBI BLAST+ 2.17.0 fetches every subject these
  reports show, and its cache can then write the taxonomy database's offset error again; LOSAT wrote
  it fewer times with exit status 0. The refusal comes before any byte; blastx, which writes each
  query batch's report as it goes, holds the document of such a run.
