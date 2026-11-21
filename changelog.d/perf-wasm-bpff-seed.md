- Web (WebAssembly): BLASTP with `-task blastp-fast` takes less time in Node, with the same reports. At
  one thread against a 6,600-protein proteome, 1,000 to 8,000 query proteins run in about 10 to 12 % less
  time and 100 proteins about 11 % less; 1 to 10 proteins or one protein of 1,000 residues about 4 to 8 %
  less (these runs last 0.04 to 0.08 s and vary). Its seed loop now takes the staged copy and the two-hit
  test without branches that native builds already use. Native builds run as before.
