- BLAST XML (`-outfmt 5`), XML2 (`14`, `16`) and JSON (`13`, `15`): a title, message, `-db` name or taxonomy
  name with a byte outside US-ASCII or a control character is written as NCBI BLAST+ 2.17.0 writes it, instead of
  being refused: JSON copies the bytes (`\u00XX` for a control character); XML decodes them as NCBI's stream does,
  writes a control character as `#` with NCBI's error message, and where NCBI's stream stops on a character it
  cannot write, writes the report up to it with NCBI's error and exit status (`16` goes on, `14` and outfmt 5's
  first query end the run with exit status 255, a later outfmt 5 query's iteration is left out). A title whose
  encoding NCBI cannot guess is written as NCBI's reports write it (`lcl|Query_N`, `No definition line`,
  `Unknown`).
