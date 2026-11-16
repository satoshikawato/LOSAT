- Speed: a BLASTN database search no longer allocates and clears a new buffer for every block of the
  database it reads; each search thread reads its blocks into one buffer. Searches of many short
  database sequences run up to about 1 % faster. Output is unchanged.
