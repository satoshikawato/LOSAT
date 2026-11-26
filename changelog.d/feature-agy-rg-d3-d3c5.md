- `losat makeblastdb -input_type asn1_txt|asn1_bin`: ported NCBI's `CDeflineGenerator` across full
  ASN.1 Bioseq inputs (`D3c-5`, building on `#1627`), splitting `report::defline` into a single-owner
  core shared with `blastdbcmd` and search reports (`x_SetFlags`, `x_SetPrefix`, `x_SetSuffix`,
  `x_AdjustProteinTitleSuffix`, `x_SetTitleFromBioSrc`, `x_SetTitleFromProtein`, `NC`/`NM`/`NR`,
  `WGS`/`TSA`/`TLS`, `PDB`, patent, segmented, and map titles, backed by a compact per-TSE feature
  index that drops raw `Seq-annot` trees immediately). Only `AutoDef` (`AutodefOptions` on an
  untitled non-TLS, non-RefSeq-contig nucleotide sequence) remains refused (`D3c-6`).
