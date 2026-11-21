- A database file that opens but whose bytes cannot be read (an I/O error of the disk, or in LOSAT Web a picked file
  that changed after it was picked) gives NCBI BLAST+'s error, `BLAST Database error: Database memory map file error -
  Cannot map file '<path>' to memory (offset=0, length=<bytes>): Input/output error` (exit status 2; the sequence file
  in the search: `BLAST engine error: Database memory map file error`, exit status 3; a header file after the outfmt 0
  prolog, as NCBI), where LOSAT said the database "is not supported" (exit status 1). The lookup files that `-taxids`
  and `-seqidlist` open give NCBI's own errors when they cannot be opened or read: the taxonomy LMDB file
  (`.ptf`/`.ntf`), and the LMDB file (`.pdb`/`.ndb`) with too many open files, say `LMDB runtime error: mdb_env_open:
  <reason>` (exit status 2), and the `.pos`, `.pot` and `.pto` files give NCBI's memory map exception (exit status
  255). blastdbcmd writes NCBI's text for a header, sequence or lookup file that it cannot read while it retrieves
  sequences (`Error: [blastdbcmd] CMemoryFile: Cannot memory map file '<path>':Permission denied`, without the search
  programs' `Database memory map file error - `). LOSAT Web gives the command line's texts for a file it cannot open
  or read; before, it said the database "is not supported", with the WebAssembly C library's words.
