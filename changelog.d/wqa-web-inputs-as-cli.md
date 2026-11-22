- LOSAT Web: a FASTA input with a `>?` gap line (for example `>?100` between two lines of a
  record) is read as on the command line, with the gap's letters in the record, where the input
  was refused ("gap lines are not supported in the browser").
- LOSAT Web: a FASTA input whose lines the command line joins (a line ending in CR alone in a file
  of LF line ends) is read as on the command line, where it could be refused ("NCBI BLAST+ reads two of its
  lines as one … such input is not supported in the browser").
- LOSAT Web: an input whose first line may be a sequence identifier that NCBI BLAST+ would fetch,
  and a record longer than 2147483647 letters, are refused with the command line's words ("… is
  not supported by LOSAT's BLASTP (start the query with a '>' defline, or set DATA_LOADERS=none in
  the [BLAST] section of .ncbirc)"), where the Web said "not supported in the browser".
- LOSAT Web: with "Additional options" `-conffile NAME` and a file that sets `DATA_LOADERS=none` in
  its `[BLAST]` section, a first line that may be a sequence identifier is read as data, as on the
  command line, where the Web refused the input. ABI v2 is version 8 (`losat_web2_scan_begin` and
  `losat_web2_register` take the input's program, role and data loaders; `losat_web2_validate`
  gives the data loaders of the argv).
- LOSAT Web: an input that is still being read when the program or the data loaders of `-conffile`
  change is read again as the new setting reads it, where its card could keep what the old
  setting read (a refusal, or records).
- LOSAT Web: the card of an input whose first line may be a sequence identifier says that the
  browser reads `DATA_LOADERS` from the file of `-conffile`, not from `.ncbirc`, and gives its
  remedies ("Put '>' before the first line", the sentence on Fetch) also when the engine's check
  refuses the input after `-conffile` is taken away.
- LOSAT Web: while another option of "Additional options" is refused, the inputs are read with
  `DATA_LOADERS` of the file of `-conffile`, where the card asked for `DATA_LOADERS=none` that the
  file already set.
