- With `-parse_deflines`, queries whose Seq-ids match without being equal (an accession with and
  without its version, a PDB entry with and without its chain) are shown in the XML and JSON formats
  (outfmt 5, 13 to 16) as NCBI BLAST+ 2.17.0 shows them: each result carries the id, title and
  length of the first query of its batch with a matching Seq-id, and against `-subject` the XML2
  and JSON formats write one report for the results of one Seq-id instead of one per query (line
  RG-KB, item 7c).
