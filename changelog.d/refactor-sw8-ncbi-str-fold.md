- Internal: the remaining copies of NCBI's `NStr` string and number functions in LOSAT (the argument
  readers' integer, genetic code, DUST and SEG conversions, the FASTA readers' space trimming and gap
  sizes, the format guess, the database date, the accession guide and `StringToBool`) now call the
  one port, `ncbi_str`. No change in behaviour.
