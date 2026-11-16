- BLASTX (and rpstblastn, the library API and LOSAT Web): a query or `-subject` FASTA title with bytes that are not
  UTF-8 (Latin-1, cp1252, overlong or other invalid sequences) is read and written as NCBI BLAST+ 2.17.0 writes it,
  instead of refusing the run: `Query=`, `# Query:` and the tabular ids keep the bytes, the archive and the search
  strategy write `#` with NCBI's error messages, and outfmt 0 shows a `-subject` record's description and heading as
  NCBI's definition line generator makes them (entities decoded, a last byte >= 0x80 dropped, `Unknown` and the
  skipped-alignment line for a title whose encoding NCBI cannot guess), which also corrects valid UTF-8 and ASCII
  titles that it changes. The default outfmt 6 and 7 are written by the same writer as a field list.
