# Inputs of the makeblastdb version 5 cases

`b2.fa`: three lines (`>x`, `>y`, `MK`) that `src/writedb/tests/makeblastdb_v5_cases.rs` case `M_ignored` gives to
`makeblastdb` with `a.fa`: an empty record `x`, then `y`; origin not recorded (added in commit 5f19ec192, 2026-10-29).

`oid_*.fa`: the inputs of the cases `O_*` (`-oid_masks`, line RG item D3b, 2026-11-21), written by hand for the probes of
`remaining-gaps-20261112/D/d3/b/probe/`: predicted RefSeq accessions (`XP_`, `XM_`, `XR_`) beside others, in one and in
several definition lines (joined by Ctrl-A), with a GI (`oid_pm_gi.fa` has a record of a GI alone, which makeblastdb 2.17.0
reads only for version 4), over several volumes (`oid_pmv.fa`), and accessions with the protein flag `P` that make NCBI's
accession guide warn (`oid_pflag2.fa`). `oid_sc.fa` (case `O_warn_all`, from
the independent review of the item): a record whose first id is no model and whose second makes the guide warn.
