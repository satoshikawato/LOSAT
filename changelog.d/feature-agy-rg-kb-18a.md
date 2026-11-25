- `blastx` and `rpstblastn` read their FASTA input (queries, `-subject`, the pre-read of the XML
  formats and LOSAT Web's ABI v1 inputs) with the reader of `blastn`, `blastp`, `tblastn` and
  `tblastx`, so BLASTX no longer keeps a FASTA reader of its own, and ABI v1 reads BLASTX input as
  the command line reads a file with the same bytes. Outputs do not change. A 35 MB `-subject` file
  takes 0.33-0.87x the time (1 or 16 threads, with or without `-parse_deflines`) with 5-16 % less
  peak memory, and large query files take 1-7 % less time (20 % with `-lcase_masking` at 16
  threads). Two LOSAT-only refusals, a sequence longer than `Int4` and a `-parse_deflines` Seq-id
  that names the other molecule, now have the other programs' wording and name the program that
  refused (`RPSTBLASTN` for `rpstblastn`).
