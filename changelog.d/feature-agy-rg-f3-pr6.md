- `blast_formatter -archive <file> -outfmt 11`: re-emits BLAST archives (`Blast4-archive`) in
  ASN.1 text, binary (`ARCHIVE_FORMAT=bin*`), XML (`ARCHIVE_FORMAT=xml`), or silent mode for
  other `ARCHIVE_FORMAT` values, across single- and multi-object `-db` and `-subject` archives
  of every search program, updating `HitlistSize` from `-max_target_seqs` / `-num_alignments`
  without pruning alignments.
