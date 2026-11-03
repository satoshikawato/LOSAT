- LOSAT Web: "Fetch from" starts at "Fastest": the first fetch of nucleotide or protein accessions asks NCBI, EBI and
  CNCB for one fixed public record (never the user's accession) and sends the accessions to the first that answers,
  among those that serve them; the tab keeps that order, and a site can still be chosen by hand. An Example's saved
  copy can be written into a folder of the user's choice ("Save to folder…", Chromium and Edge) or downloaded
  ("Download", other browsers), to be opened again with "Open FASTA files…" (line WF-3, docs/web/fetch_design.md §10).
