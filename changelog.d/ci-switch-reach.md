- CI: `LOSAT_X_REACH=<dir>` records which side each conditional performance switch took (its own path or the
  reference path while on; `x_reach!` at the decision points of 75 switches, output unchanged), and
  `LOSAT/tests/switch_reach_check.py` fails a pull request when no fast fixture reaches a side that
  `LOSAT/tests/switch_reach_expected.tsv` requires; 25 sides no fixture reaches yet are listed there as gaps.
