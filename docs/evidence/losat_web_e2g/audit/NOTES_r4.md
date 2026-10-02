# E2g round 4 (r4) notes
Binary AC8 (f4057718a). Files: CASES.tsv (all CLI compared runs with class: X 976, MINI 360, R4 2038, R4G 584 = 3958), CASES_*_class.tsv, DIFFS_R4.tsv,
WEB.tsv (966 web-path comparisons: webh harness = run_local_blastn with in-memory writers for formats 0/6/7 + diagnostics + observer, vs NCBI -outfmt 0/6/7 and stderr),
cases4.py/cases4g.py (generators), webcmp.py, webh/ (harness crate, built in e2g-r4-webh), classify4.py.
Result: no defect. 3280 identical; 100 approved exc1 (-evalue nan, merged/sep text), 12 exc2, 187 exc3 (6/7 write failure), 379 explicit rejections.
