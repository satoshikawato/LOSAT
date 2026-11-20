- LOSAT Web: BLASTN's `-window_masker_taxid`. The "Species-specific repeats filter for:" row gains a
  "Taxonomy ID" box beside the "Counts file" choice; "Download" fetches that taxonomy ID's
  `wmasker.obinary` from ftp.ncbi.nlm.nih.gov only when pressed, keeps it in the browser (a stopped
  download continues with "Continue"), and searches read it from there with the same output as the
  command line. The run details show the taxonomy ID, the file's date and its SHA-256.
  The engine looks for the counts file of a taxonomy ID through the same counts source as
  `-window_masker_db` (the file system for the command line, unchanged; the registered files for the
  browser, under `windowmasker/<taxid>/wmasker.obinary`); ABI v2 stays version 2.
