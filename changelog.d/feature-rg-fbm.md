- psiblast and tblastn `-in_pssm` decode the BitMagic data of a Seq-table in the PSSM query's annotations
  (`bit-bvector` values, a `bit-set-bvector`'s sparse index) as NCBI BLAST+ 2.17.0's object manager decodes it, where
  LOSAT refused it: the table's rows, values and messages come as in NCBI, and NCBI's `BM::Invalid serialization
  format` error ends the run with exit status 255 (or is the message of a `Sorted, max length` or `Seq-table location`
  column, after which the run goes on). Malformed data that NCBI reads past, decodes with undefined behaviour or from
  memory it did not write stays refused, with its reason (line RG item Q-BM, decided by the Owner on 2026-11-16). A GAP
  block that is not one is decoded, and its rows answered, where NCBI's scalar, SSE2 and SSE4.2 builds of BitMagic
  agree, and refused only where they differ, as NCBI's result then depends on the machine it was built for.
  The decoder is a Rust port of parts of BitMagic 8.0.1 (Apache-2.0, bundled in NCBI BLAST+); its notice and licence
  text are in THIRD_PARTY_NOTICES, section 10 (BitMagic), and the License sentence of the READMEs now names BitMagic.
