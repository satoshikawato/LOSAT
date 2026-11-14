- `-db` with version 5 databases made with and without `makeblastdb -parse_seqids` side by side (`-db "a b"`, or an alias
  file that lists them), in every search program: LOSAT now searches such a list as NCBI BLAST+ does, where it refused the
  list when it opened it. Every run whose reports look no parsed Seq-id up in a database without the accession table is
  written as NCBI writes it: the Seq-annots and the archive (`-outfmt 8`, `9`, `11`, `12`), hits only in the databases
  without parsed Seq-ids or in ones made with `-input_type blastdb`, an alias file whose `FIRST_OID`-`LAST_OID` range or a
  `-taxids` list keeps no parsed sequence. Where a reported parsed Seq-id would be looked up in a database without the
  table, NCBI stops with `BLAST Database error: DB contains no accession info.` and LOSAT refuses the run, now at the
  reports instead of at the opening.
