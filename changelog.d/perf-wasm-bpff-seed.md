- Web (WebAssembly): BLASTP with `-task blastp-fast` takes less time in Node, with the same reports. At
  one thread against a 6,600-protein proteome, 1,000 or 3,000 query proteins run in about 20 % less
  time, and 1 to 100 proteins or one protein of 1,800 residues about 5 to 10 %. Its seed loop now takes
  the staged copy and the two-hit test without branches that native builds already use. Native builds
  run as before.
