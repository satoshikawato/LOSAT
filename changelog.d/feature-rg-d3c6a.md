- `makeblastdb -input_type asn1_txt` and `asn1_bin` now give NCBI BLAST+ 2.17.0's messages about the
  annotations of each object (for example `Invalid 'ids' size in denseg`), its exceptions (a
  Bioseq-set id twice), and its output after `Failed to parse location of`: later messages are their
  headers alone, and a later empty location ends the build with exit status 0, keeping the volumes
  written, as NCBI's does. An object that NCBI reports with a memory address (`Empty region in 0x...`,
  an element stated as `NULL`) or crashes on is refused where NCBI meets it and leaves no files; such
  objects were built without a message before, and some descriptors stated as `NULL` that NCBI does
  not read were refused. After `Failed to parse location of`, an element stated as `NULL` in a
  `Delta-ext` or a binary definition line header ends the build as NCBI's does, keeping the volumes
  written. An error that NCBI reports with a standard library message (for example
  `Error: vector::reserve` of an object that asks for more elements than a list can hold) now keeps
  the volumes written so far, as NCBI's does; they were removed before.
- A `makeblastdb` build that LOSAT refuses now leaves no file of the build, whatever the input
  (FASTA, ASN.1, or a database with `-input_type blastdb`) and wherever the refusal comes, the lock
  files of the version 5 lookup tables included; before, a `-gi_mask` refusal kept its GI-mask files
  and alias file, and a refusal when a version 5 build closes (a sequence with GI ids only between
  sequences with other ids, `MAX_LMDB_TXN_ENTRY=0`) kept some of the build's files.
