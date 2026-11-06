- LOSAT Web: "Additional options" takes `-outfmt` of every single-file format that the CLI writes:
  5 (BLAST XML), 8 and 9 (Seqalign, text and binary ASN.1), 10 and 20 (CSV, with fields too),
  11 (BLAST archive), 12 (Seqalign JSON), 15 and 16 (single-file BLAST JSON and
  XML2), and BLASTN's 17 (SAM, with its tags). The run still writes outfmt 0, 6 and 7 and adds that
  output, byte for byte the CLI's output of the same command line (SAM's `@PG` line shows the
  run's command; the archive holds the whole records and their deflines as read, with
  `-query_loc` and `-subject_loc` too). The Outputs show it as "outfmt 5" (outfmt 9 is saved, not
  shown), Export saves it with its format's extension (`…outfmt5.xml`, `…outfmt9.asnb`,
  `…outfmt10.csv`, `…outfmt17.sam`), and "Download All" names it as BLAST+'s `-outfmt` help does
  ("BLAST XML (outfmt 5)"). A run with such an output is "Not yet supported" until the support records compare
  it. ABI v2 stays version 2.
