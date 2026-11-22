# Inputs of the makeblastdb version 5 cases

`b2.fa`: three lines (`>x`, `>y`, `MK`) that `src/writedb/tests/makeblastdb_v5_cases.rs` case `M_ignored` gives to
`makeblastdb` with `a.fa`: an empty record `x`, then `y`; origin not recorded (added in commit 5f19ec192, 2026-10-29).

`oid_*.fa`: the inputs of the cases `O_*` (`-oid_masks`, line RG item D3b, 2026-11-21), written by hand for the probes of
`remaining-gaps-20261112/D/d3/b/probe/`: predicted RefSeq accessions (`XP_`, `XM_`, `XR_`) beside others, in one and in
several definition lines (joined by Ctrl-A), with a GI (`oid_pm_gi.fa` has a record of a GI alone, which makeblastdb 2.17.0
reads only for version 4), over several volumes (`oid_pmv.fa`), and accessions with the protein flag `P` that make NCBI's
accession guide warn (`oid_pflag2.fa`). `oid_sc.fa` (case `O_warn_all`, from
the independent review of the item): a record whose first id is no model and whose second makes the guide warn.

`raw_*`: the inputs of the cases `B_*` (`-input_type blastdb`, line RG item D3c, 2026-11-21), written by hand for the probes
of `remaining-gaps-20261112/D/d3/c/probe/`: `raw_n.fa` (five nucleotides: IUPAC ambiguity codes, lower case, a run of N,
RefSeq, GenBank with a GI, a local id), `raw_n2.fa` (the same sequences under other accessions), `raw_m.fa` (one record of
three definition lines joined by Ctrl-A, and two more), `raw_p.fa` (three proteins with X, B and Z), the taxid maps
`raw_tax.map`, `raw_ptax.map`, `raw_part.map` (one accession) and `raw_none.map` (none of the input's), and the mask files
`raw_n.dust.asnb`, `raw_n.dust2.asnb`, `raw_n2.dust2.asnb` (dustmasker 2.17.0, `-outfmt maskinfo_asn1_bin`, window 64
and window 32 as their options record) and `raw_p.seg.asnb` (segmasker 2.17.0, the default options). `raw_pig/`: `e4.asn`, a Bioseq-set written by
hand whose definition lines (an `ASN1_BlastDefLine` user object) carry membership bits and a PIG in `other-info`, and the
database `db.*` that makeblastdb 2.17.0 built from it (`-input_type asn1_txt -dbtype prot -title E4`; LOSAT does not read
ASN.1 input yet).
