- `losat rpsbproc` holds a long string of its input once, so an input with a long title or sequence needs
  less memory, and it reads a BlastOutput XML input with fewer instructions. The query's title and residues move
  out of the archive's Bioseq instead of being copied, the report writes a long string straight to the output
  instead of copying it into its own buffer first, and the stream of an XML input no longer keeps the character
  data it has read. The XML reader takes the text that needs no conversion, and the text of a CDATA section, in
  runs instead of one byte at a time, and the characters that it converts (entities, Latin-1 or cp1252 bytes)
  are no longer an allocation each. Outputs do not change. On an archive whose query title is one line of
  64 MiB, max RSS falls from 200 MiB to 71 MiB (NCBI rpsbproc 0.5.1: 161 MiB); the same holds for a title taken
  from the query's protein feature, and a sequence of 64 Mi residues falls from 136 MiB to 71 MiB. On an XML
  input whose `Iteration_query-def` is 64 MiB, max RSS falls from 199 MiB to 71 MiB (NCBI: 160 MiB) and the run
  takes 0.08 s instead of 0.49 s. With an `Iteration_query-def` of 2 MiB, the run takes 0.03 to 0.24 of the
  instructions for text in runs (letters, tabs, UTF-8, line ends between them), 0.18 to 0.75 for CDATA sections
  (one long section, sections holding `]`, or many small ones), and 0.55 to 0.89 for text of entities, Latin-1 or
  cp1252 bytes, line ends or NULs. An XML input of ordinary size runs 4 % fewer instructions (13 % with a title of
  4,000 bytes of UTF-8), one of 5,000 hits or 5,000 queries whose strings are empty 0.1 to 0.3 % fewer, and one of
  20,000 hits half the instructions. Archive inputs of ordinary size run no more instructions than before, and an
  archive of 5,000 queries of one residue 4 to 6 % fewer.
