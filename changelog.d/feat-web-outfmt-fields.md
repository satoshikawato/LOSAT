- LOSAT Web: "Additional options" takes `-outfmt` with fields for outfmt 6 and 7, as the CLI does
  (`-outfmt '6 qseqid sseqid pident'`, `-outfmt '7 std qcovs'`, `delim=`), in every mode. The run
  still writes outfmt 0, 6 and 7, and adds the output of that `-outfmt`, byte for byte the CLI's
  output of the same command line: the Outputs show it as "outfmt 6 qseqid sseqid pident", Export
  and "Download All" save it as `…outfmt6-custom.txt`, and Run details gives its command. The
  results views keep reading outfmt 6. The formats that the browser does not write (5, 10, …) say
  so ("other formats need the CLI") instead of "-outfmt is managed by LOSAT", and `-num_threads`
  says "-num_threads is set in the form: Threads". A run with such an output is "Not yet supported"
  until the support records compare it. The engine's BLASTX for the browser writes custom fields
  as its command line does. ABI v2 gains stream 4 and stays version 2.
