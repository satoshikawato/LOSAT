- LOSAT Web: a box of accessions that no one site serves (a GenBase `C_` accession with a UniProt one) is no longer
  stopped: under "No one site has all of these accessions" it has one Fetch per group of accessions, each naming the
  one site it sends to ("Fetch C_AA004835.1 from CNCB", "Fetch 2 accessions from NCBI", with NCBI or EBI to choose).
  A press sends only that group's accessions, only to that site; the box is written, in its order, once every group
  has come, and a failed group can be fetched again alone (line WF-5, docs/web/fetch_design.md §11).
