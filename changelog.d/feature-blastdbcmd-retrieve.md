- `losat blastdbcmd`: the retrieval of sequences as NCBI BLAST+ 2.17.0 blastdbcmd writes it, byte for byte:
  `-entry all` (a database to FASTA), `-entry` and `-entry_batch` of `gnl|BL_ORD_ID|n` ids in a version 4
  database and of any id in a database built without `-parse_seqids` (NCBI's answers), every `-outfmt` specifier
  but the ASN.1 texts `%d` and `%b` (`%f`, `%s`, `%a`, `%i`, `%g`, `%o`, `%t`, `%l`, `%h`, `%T`, `%X`, `%N`,
  `%C`, `%S`, `%L`, `%B`, `%K`, `%P`, `%e`, `%n`, `%%`, with NCBI's errors for an unknown specifier, `%m` and
  `-mask_sequence_with` on a database without masks), `-range`, `-strand`, `-line_length`, `-ctrl_a`,
  `-target_only`, `-get_dups`, `-long_seqids`. Two outputs that NCBI gets wrong are written correctly (approved
  differences): `-long_seqids` with `-range` on the minus strand writes the reverse complement under NCBI's `:c`
  header (NCBI writes the forward letters), and `gnl|BL_ORD_ID|n` of a version 4 database in volumes is the one
  sequence `n` of the database (NCBI writes the sequence `n` of every volume). The lookup of an accession, a name
  or a GI in the tables of a database built with `-parse_seqids` (LMDB, ISAM), `-outfmt %d` / `%b`, `-range ""`
  (NCBI crashes), `-list`, the taxonomy options, `-metadata` and `-ipg` are refused.
