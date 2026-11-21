- Web engine (npm package `losat`, ABI version 12): a search of `-db` runs in the browser with the same
  bytes as the command line, with its filters (`-seqidlist`, `-taxids`, `-taxidlist`, `-gilist` and their
  negatives) and BLASTN's database masks. The host gives the run the database's files as `File`s
  (`EngineRunRequest.files`, named as the command line names them from its directory), and the engine
  reads each file in place, only the bytes that it needs: no database is copied into memory or into the
  browser's storage, and nothing is uploaded. All threads of a threaded search read the same files.
  LOSAT Web's form does not offer a database yet. `-remote` and `-entrez_query` stay refused, as on the
  command line, and a search must fit in the 4 GiB of WebAssembly memory.
