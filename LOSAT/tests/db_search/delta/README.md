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

PS13a fix round 2 (2026-11-24) adds `q_latin.faa` (the query `t5/hd_q1.faa` with the byte 0xE9 in its title, which the
PSSM file's ASN.1 text writes with a post); its rows search the databases `../t5/phd_vstr` and `../t5/phd_utf8` (a header
that NCBI cannot read, and a UTF8String title) with the queries `../t5/hd_q*.faa`.

PS13b (2026-11-27) adds the queries of `-query_loc` and `-show_domain_hits`: `q_tail.faa` (a serpin then a kinase, 1094
letters), `q_head.faa` (a kinase then a serpin), `q_dd.faa` (serpin letters 1-200, a kinase, a serpin; 1294 letters),
`q_tt.faa` (two records like `q_tail.faa`), `q_tp.faa` (`q_tail.faa` and the kinase of `q_pk.faa`), `q_three.faa`
(`q_serpin.faa`, `q_pk.faa`, `q_tail.faa`) and `q_empty.faa` (no bytes), letters of `q_serpin.faa` and `q_pk.faa`.
The `o_*.faa` files are the oracles of the rows `db.dl.ql.or.*`: each record is a record of the named query cut to the
range of its name (`o_serp_100_364.faa` = `q_serpin.faa` letters 100-364, `o_m_100_500.faa` = each record of
`multi.faa` cut to 100-500, as NCBI's `-query_loc` cuts them), under the record's definition line, so NCBI's run of
the range's letters as a query of their own (PD-LOSAT-NCBI-DEFECTS exception 35).
