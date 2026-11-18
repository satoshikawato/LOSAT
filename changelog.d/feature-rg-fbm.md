- psiblast and tblastn `-in_pssm` decode the BitMagic data of a Seq-table in the PSSM query's annotations
  (`bit-bvector` values, a `bit-set-bvector`'s sparse index) as NCBI BLAST+ 2.17.0's object manager decodes it, where
  LOSAT refused it: the table's rows, values and messages come as in NCBI, and NCBI's `BM::Invalid serialization
  format` error ends the run with exit status 255 (or is the message of a `Sorted, max length` or `Seq-table location`
  column, after which the run goes on). Malformed data that NCBI reads past, decodes with undefined behaviour or from
  memory it did not write stays refused, with its reason (line RG item Q-BM, decided by the Owner on 2026-11-16). A row
  in a GAP block that is not one is answered where NCBI's scalar, SSE2 and SSE4.2 builds of BitMagic agree, and refused
  only where they differ, as NCBI's answer then depends on the machine it was built for.
