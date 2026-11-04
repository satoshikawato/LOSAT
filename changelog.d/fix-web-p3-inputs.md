- LOSAT Web, inputs: a file or pasted text that looks like GenBank (`LOCUS`), EMBL / UniProt (`ID   `)
  or FASTQ (`@` records) gets a warning on its card, since LOSAT reads it as FASTA text as the command
  line does; an input that starts with a byte order mark (UTF-8, or UTF-16) says so after the reader's
  error, and a pasted one can drop the mark. A subrange whose start is after its stop, or equal to it,
  is refused under the fields with the meaning of the engine's message (the engine refuses both in
  every mode). Reloading or closing the tab asks first while the form holds inputs or a run has
  results, not only while a search runs, since nothing is kept for a later visit. Searches and
  outputs are unchanged.
