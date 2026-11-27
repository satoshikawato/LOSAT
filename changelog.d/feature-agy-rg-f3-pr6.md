- `blast_formatter -archive <file> -outfmt 11`: writes the BLAST archive (`Blast4-archive`) again,
  read in ASN.1 text, binary ASN.1 or XML, and written in ASN.1 text, binary
  (`ARCHIVE_FORMAT=bin*`), XML (`ARCHIVE_FORMAT=xml`) or nothing for other `ARCHIVE_FORMAT`
  values, with the search programs' archive writer; for single- and multi-object `-db` and
  `-subject` archives of blastn, blastp, blastx, tblastn, tblastx, psiblast, deltablast,
  rpsblast and rpstblastn, with `HitlistSize` set from `-max_target_seqs` / `-num_alignments`
  and the alignments kept whole. Refused: archives of the mapper, archives whose queries are a
  list of sequence locations (`seq-loc-list`), values that LOSAT does not write yet (a
  UTF8String, a NULL value, an OCTET STRING), and a `FilterString` number on which NCBI BLAST+
  has no defined result.
