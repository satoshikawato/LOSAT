# deltablast databases and queries (line PS, series PR, item PS13a)

NCBI BLAST+ 2.17.0's unit-test domain databases (public domain, `c++/src/algo/blast/unit_tests/api/data`): `deltatest`
(three serpin domains, with `.freq`, `.wcounts`, `.obsr`) and `deltatest_nocbs` (the same without `.freq`); `prof_rps` is
an RPS-only database (no `.wcounts`, `.obsr`); `dbs` is a protein database of 380 kB (`dbs.faa`) that the queries search;
`q_*.faa`, `serpin3.faa`, `multi.faa` are queries and `subj_*.faa` subject sets. The rows `db.dl.*` of
`../../db_search_manifest/deltablast.tsv` run `deltablast -db tests/db_search/delta/dbs -rpsdb tests/db_search/delta/deltatest`.

PS13a fix round 1 (2026-11-24) adds the queries of its NCBI 2.17.0 probes: `q_batch2.faa` and `q_batch3.faa` (more
than NCBI's 10000 letters of a query batch), `q_gap2.faa`, `q_lc.faa` (lower case), `q_none.faa` (no alignments),
`q_x1.faa`, `q_x2.faa`, `q_xs.faa` (a query of X alone, after and before a query with letters) and `q_e1.faa`,
`q_e1s.faa`, `q_e2.faa`, `q_ee.faa`, `q_ees.faa`, `q_see.faa` (queries without letters, alone and mixed).
