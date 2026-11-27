- `makeblastdb -input_type asn1_txt` and `asn1_bin` now give NCBI BLAST+ 2.17.0's messages about the
  annotations of each object (for example `Invalid 'ids' size in denseg`), its exceptions (a
  Bioseq-set id twice), and its output after `Failed to parse location of`: later messages are their
  headers alone, and a later empty location ends the build with exit status 0, keeping the volumes
  written, as NCBI's does. An object that NCBI reports with a memory address (`Empty region in 0x...`,
  an element stated as `NULL`) or crashes on is refused where NCBI meets it and leaves no files; such
  objects were built without a message before, and some descriptors stated as `NULL` that NCBI does
  not read were refused.
