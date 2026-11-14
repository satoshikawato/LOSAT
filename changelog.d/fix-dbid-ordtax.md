- `-db` on a database without parsed Seq-ids whose search leaves sequences out (`-taxids`, `-taxidlist`,
  `-negative_taxids`, `-negative_taxidlist`, or an alias file's `FIRST_OID`/`LAST_OID`), in every search
  program: BLAST XML and XML2 (outfmt 5, 14, 16), JSON (13, 15), the Seq-annots (8, 9, 12), the BLAST archive
  (11), SAM (17) and the organism report (18) now name a sequence `gnl|BL_ORD_ID|n` (`BL_ORD_ID:n` in SAM, the
  accession `n`) with its OID n, as NCBI BLAST+ 2.17.0 does. LOSAT wrote its number among the searched
  sequences there with exit status 0 (for example `gnl|BL_ORD_ID|0` for OID 1 under `-taxids`); outfmt 0 and
  the tabular formats were already right.
