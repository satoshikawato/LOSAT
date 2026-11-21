- blastdbcmd: `-outfmt %d` writes the definition lines of each sequence as ASN.1 text
  (`Blast-def-line-set`, or with `-target_only` the target's `Blast-def-line`) and `-outfmt %b`
  the sequence as an ASN.1 text Bioseq (ids, title, the definition lines' header, the taxonomy
  names, the sequence data), byte for byte as NCBI BLAST+ 2.17.0 writes them, with
  `-entry`, `-entry_batch`, `-ipg`, `-taxids`, `-target_only`, `-get_dups` and `-out`, on
  protein and nucleotide databases of version 4 and 5. Both were refused before. As in NCBI,
  `%b` and `%d` replace the rest of the format, `-range` and `-strand` do not change `%b`, and
  two of `%f`, `%b` and `%d` together are an error.
