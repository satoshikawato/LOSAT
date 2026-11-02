- The variables `LOSAT_X_BPFSEED`, `LOSAT_X_BPFSTAGE` and `LOSAT_X_BPFLEAN` are gone: `blastp -task blastp-fast`
  always takes the loops they chose by default (on a table that takes the branch-free cell list, natively: the
  branch-free two-hit test, and the staged hits with the list two offsets per step for a subject of a word or more;
  WebAssembly keeps its loop). Reports do not change. `LOSAT_X_BPCLUTSCANBF=0` still gives the loop with the branch.
