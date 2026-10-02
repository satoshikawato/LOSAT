# E2g round 3, part A (r3a) notes
Binary AC7 (ad5fa9c85). Harness: lib.py (copy of b/lib.py; LOSAT=AC7, REF dropped), loadall.py+runmain.py (3728 ids of b/CASES.tsv, saved as CASES_R1.tsv),
CASES.tsv (new run), CASES_class.tsv (round-3 classes), CHANGES.tsv (outcome/class changes vs round 1), classify3.py.
Extra sets: CASES_X.tsv (976 merged/devfull/-out devfull etc.), CASES_MINI.tsv (360 K-A scoring x warnings x devfull), CASES_EXTRA.tsv (163 T7ac1 from AC1 round),
CASES_EV.tsv (127 -evalue inf/nan forms), CASES_LONG.tsv (5 reruns with 3600 s timeout), exc1.tsv (exception-1 reference runs).
Defect: ka-scoring + outfmt 0 write failure + title warning -> extra FASTA-Reader line (diff_x, diff_mini).
