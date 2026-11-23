- The library's `Hit::subject_title` of a FASTA subject (its title after the first word) has one producer
  (`blastinput::fasta_reader::hit_subject_title`) in place of five copies in BLASTN, BLASTP, BLASTX,
  TBLASTN and TBLASTX. Values and outputs do not change.
