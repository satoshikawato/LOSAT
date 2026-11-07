- LOSAT Web: when the first line of a pasted input is an accession that NCBI BLAST+ would fetch as a
  sequence identifier and the input cannot be read, the input's card now points to Fetch ("To fetch
  accessions, put only the accessions in the box, one per line."), and, for any such first line (an
  accession, a GI number, a PDB id), offers "Put '>' before the first line" in place of "Add the line
  >pasted_query", which left the identifier as the first line of data and failed again or turned it
  into residues. A pasted list of accessions on one line, separated by spaces or commas, offers a
  "Put each accession on its own line" button, which rewrites the box so that Fetch takes them, also
  when the list reads as a sequence (accessions without versions). The buttons that rewrite the box
  ("Add the line", "Put '>'", "Remove the byte order mark") keep the keyboard focus on the input's
  card instead of dropping it to the page. The command line's message is unchanged.
