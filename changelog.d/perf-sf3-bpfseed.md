- The variables `LOSAT_X_BPFSEED`, `LOSAT_X_BPFSTAGE` and `LOSAT_X_BPFLEAN` are gone: `blastp -task blastp-fast`
  always takes the loops they chose by default (on a table that takes the branch-free cell list, natively: the
  branch-free two-hit test, and the staged hits with the list two offsets per step for a subject of a word or more;
  WebAssembly keeps its loop). Reports do not change. `LOSAT_X_BPCLUTSCANBF=0` still gives the loop with the branch.
- `blastp -task blastp-fast` finds the subject words with lookup hits without a branch per word on every query's compressed lookup table, not only from 500 used cells per 10,000: one query of 10 to 739 residues against UniRef50 volume 00 takes about 0.6-0.7 of the CPU cycles at one thread. `LOSAT_X_BPCLUTSCANBFMIN` is gone; `LOSAT_X_BPCLUTSCANBF=0` still gives the loop with the branch. Reports do not change.
