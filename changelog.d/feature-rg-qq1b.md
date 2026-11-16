- psiblast and tblastn `-in_pssm` read a PSSM whose query Seq-id is not one local id (the `gi`
  of a PSSM made with `-parse_deflines`, swissprot, pdb, patent and every other kind, and several
  ids), shown as NCBI BLAST+ 2.17.0 shows them with or without `-parse_deflines`, the first id
  naming the query in the alignments, later PSSMs and warnings (line RG, Q1b). A database sequence
  with one of those ids is refused as for a `-parse_deflines` query.
