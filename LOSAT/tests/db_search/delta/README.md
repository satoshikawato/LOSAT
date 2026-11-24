# deltablast databases and queries (line PS, series PR, item PS13a)

NCBI BLAST+ 2.17.0's unit-test domain databases (public domain, `c++/src/algo/blast/unit_tests/api/data`): `deltatest`
(three serpin domains, with `.freq`, `.wcounts`, `.obsr`) and `deltatest_nocbs` (the same without `.freq`); `prof_rps` is
an RPS-only database (no `.wcounts`, `.obsr`); `dbs` is a protein database of 380 kB (`dbs.faa`) that the queries search;
`q_*.faa`, `serpin3.faa`, `multi.faa` are queries and `subj_*.faa` subject sets. The rows `db.dl.*` of
`../../db_search_manifest/deltablast.tsv` run `deltablast -db tests/db_search/delta/dbs -rpsdb tests/db_search/delta/deltatest`.
