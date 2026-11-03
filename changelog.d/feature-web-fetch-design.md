- LOSAT Web: design of the fetch of accessions and Example data (docs/web/fetch_design.md): the promise becomes "research
  data never leaves the browser; a fetch the user asks for sends only the accession or the file's URL to NCBI or EBI";
  the CSP's `connect-src` names those hosts. No fetch is made yet (lines WF-1 to WF-3 follow).
