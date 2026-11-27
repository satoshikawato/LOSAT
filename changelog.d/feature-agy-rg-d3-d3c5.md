- `losat makeblastdb -input_type asn1_txt|asn1_bin`: ported NCBI's `CDeflineGenerator` across full
  ASN.1 Bioseq inputs (`D3c-5`, building on `#1627`), splitting `report::defline` into a single-owner
  core shared with `blastdbcmd` and search reports (`x_SetFlags`, `x_SetPrefix`, `x_SetSuffix`,
  `x_AdjustProteinTitleSuffix`, `x_SetTitleFromBioSrc`, `x_SetTitleFromProtein`, `NC`/`NM`/`NR`,
  `WGS`/`TSA`/`TLS`, `PDB`, patent, segmented, and map titles, backed by a compact per-TSE feature
  index that drops raw `Seq-annot` trees immediately). Feature locations are read as NCBI's
  object manager reads them: `IsSameBioseq` (`CSeq_id::Match`, the Seq-id handle, then the
  scope) for the strand, the product and overlaps, the scope's weak matching (an accession
  without a version, a name, another text alternative) to find a location's Bioseq, and the
  reverse matching of a Bioseq's own Seq-ids for the features it has. `AutoDef`
  (`AutodefOptions` on a nucleotide sequence without a title, not TLS and not a RefSeq contig)
  remains refused (`D3c-6`), and so is a title made from the features that changes with the
  Bioseq NCBI picks for a location by the order of handles in memory (LOSAT makes the title
  with each possible Bioseq and writes it when they agree). PDB Seq-ids and a `pir` or `prf`
  accession alone match as in NCBI, with the handles of segments and references that earlier
  titles resolve. An untitled protein (or `NM`/`NR` sequence) whose segment or reference no
  input sequence matches ends with NCBI's `Cannot resolve` exception (exit status 255), and a
  segmented set whose master has such a segment gets NCBI's two `null scope pointer` errors.
