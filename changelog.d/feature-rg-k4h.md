- BLAST XML (`-outfmt 5`), XML2 (`14`, `16`) and JSON (`13`, `15`): a title, message, `-db` name or taxonomy
  name with a byte outside US-ASCII or a control character is written as NCBI BLAST+ 2.17.0 writes it, instead of
  being refused: JSON copies the bytes (`\u00XX` for a control character); XML decodes them as NCBI's stream does,
  writes a control character as `#` with NCBI's error message, and where NCBI's stream stops on a character it
  cannot write, writes the report up to it with NCBI's error and exit status (`16` goes on, `14` and outfmt 5's
  first query end the run with exit status 255, a later outfmt 5 query's iteration is left out). A title whose
  encoding NCBI cannot guess is written as NCBI's reports write it (`lcl|Query_N`, `No definition line`,
  `Unknown`).
- When such a report ends the run (`14`, outfmt 5's first query), the `-export_search_strategy` file (a FIFO or a pipe too) holds the
  requests NCBI wrote, up to the query batch of that report, and the taxonomy database's messages of that report
  come before the error, as in NCBI; with the prefetch off (`PRE_FETCH_SEQS_LIMIT`), the messages that the loads of
  the subjects of an outfmt 5 iteration left out or cut would have written come at NCBI's next load that makes the
  same lookup, of whichever subject and in whichever later query or batch (a later batch's prefetch
  too), or not at all.
  A library or Web caller gets these messages on its diagnostics, not in the error's text.
- With a version 5 list of databases with and without the accession table and a taxonomy database that writes a
  message, a later query batch whose prefetch stops (NCBI's `BLAST Database error`) after an outfmt 5 iteration that
  NCBI's stream left out at a query title outside US-ASCII is refused before any output where that batch's other
  loads look a tax id up (`blastp`, `blastn`, `tblastn`, `tblastx`): NCBI writes the messages that the left-out loads
  did not write before the error or not, as the order of its loads decides.
