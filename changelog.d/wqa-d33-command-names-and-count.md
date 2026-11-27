- Web: a list file of the Database filters (`-seqidlist`, `-negative_seqidlist`, `-gilist`, `-taxidlist` and the
  others) whose name starts or ends with a space was run, and shown in Run details, under the name without
  the space (`-seqidlist ids.txt`), so the shown command failed on the command line. The file is now run
  and shown under its own name, quoted as a query file is (`-seqidlist 'ids.txt '`).
- Web: after a help or version word (`-version`, `-h`) is typed in Additional options, the Run LOSAT button
  no longer says "(2 runs)" for a moment: it names no count until LOSAT has read the new words, then the
  count of the argv it will run.
