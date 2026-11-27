- `losat rpsbproc` holds a long string of its input once, so an input with a long title or sequence needs
  less memory, and a BlastOutput XML input is read faster. The query's title and residues move out of the
  archive's Bioseq instead of being copied, the report writes a long string straight to the output instead of
  copying it into its own buffer first, and the XML input reads the query's id and title into their own
  strings and no longer keeps the text it has read. Outputs do not change. On an archive whose query title is
  one line of 64 MiB, max RSS falls from 199 MiB to 71 MiB (NCBI rpsbproc 0.5.1: 161 MiB); the same holds for a
  title taken from the query's protein feature, and a sequence of 64 Mi residues falls from 136 MiB to 71 MiB
  (NCBI: 161 MiB). On a BlastOutput XML input whose `Iteration_query-def` is 64 MiB, max RSS falls from
  199 MiB to 71 MiB (NCBI: 160 MiB) and the run takes 0.09 s instead of 0.49 s (0.14 of the instructions); an
  XML input of ordinary size runs 3 % fewer instructions, and the archive inputs of ordinary size are
  unchanged.
