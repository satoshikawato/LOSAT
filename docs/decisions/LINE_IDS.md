# Line and stage IDs

Code comments, tests, changelog entries and documents in this repository cite the maintainer's work
lines and their stages by short IDs: "line DBID, N1b", "DI-8f", "OF8a", "plan v4 line R", "SF S3".
The working records of those lines (plans, state, decision logs) are kept outside the repository.
This file gives each cited ID one row: what it is, its state, its main pull requests, and the document
in this repository that holds its fullest record ("external" when there is none).

As of dev 193657f3a (2026-11-06): 154 IDs cited in the tree, 99 of them in `CHANGELOG.md`, `changelog.d/`
and `docs/`. State: 142 done, 12 in progress, 0 stopped.

- `python3 scripts/line_ids.py list --all` lists the IDs the tree cites, with the number of citing
  lines and the first one; `python3 scripts/line_ids.py check` names the IDs that `CHANGELOG.md`,
  `changelog.d/` and `docs/` cite without a row here.
- State: done (its work merged or the stage closed), in progress (work remains), stopped (abandoned,
  folded into other work or superseded; the meaning says which).
- Families with a registry of their own are not repeated here: DW-n (decision log of
  [`docs/losat_web_gui_plan.md`](../losat_web_gui_plan.md)), TD-n (the same plan, section 0.5), REQ-nn
  ([`docs/web/requirements_trace.tsv`](../web/requirements_trace.tsv)), PD-* ([`docs/product_decisions/`](../product_decisions/))
  and LIB-n ([`docs/losat_library_surface_plan.md`](../losat_library_surface_plan.md)).
- The Owner's decisions that the repository cites are in [`OWNER_DECISIONS.md`](OWNER_DECISIONS.md).

## Engine and command line

### BDC: `losat blastdbcmd`

| ID | Meaning | State | Main pull requests | Record in this repository |
| --- | --- | --- | --- | --- |
| `BDC` | The line that added `losat blastdbcmd`, NCBI BLAST+ 2.17.0 blastdbcmd as a LOSAT subcommand: -info, sequence retrieval, id lookups, taxonomy, masks, several -db names, batch files. | done | #616, #647, #665, #676, #681, #695, #704, #728 | [docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md](../../docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md) |
| `BDC2a` | BDC stage 2a: `losat blastdbcmd` retrieves sequences (-entry, -entry_batch, every -outfmt, -entry all) from databases without parsed ids, with version 4 ISAM lookup. | done | #647 | [docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md](../../docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md) |
| `BDC2b` | BDC stage 2b: `losat blastdbcmd` looks up ids in version 5 databases (LMDB accession index, GI ISAM). | done | #665 | [docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md](../../docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md) |
| `BDC2c` | BDC stage 2c: `losat blastdbcmd` looks up ids in version 4 databases built with -parse_seqids (string ISAM), and records the evidence of items 1 and 2 of the deviations decision. | done | #676 | [docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md](../../docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md) |
| `BDC3a` | BDC stage 3a: `losat blastdbcmd` taxonomy and discovery options: -taxids, -taxidlist, -tax_info, -metadata, -list, -ipg. | done | #681 | [docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md](../../docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md) |
| `BDC3b` | BDC stage 3b: `losat blastdbcmd` reads the masks of makeblastdb -mask_data (-mask_sequence_with, %m). | done | #695 | [docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md](../../docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md) |
| `BDC4` | BDC stage 4: `losat blastdbcmd -db "a b"` reads several database names as one virtual database. | done | #704 | [docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md](../../docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md) |
| `BDC5` | BDC stage 5: `losat blastdbcmd` reads the lines of a -entry_batch file as NCBI does on each platform (line ends). | done | #728 | [LOSAT/src/blastdbcmd/args.rs](../../LOSAT/src/blastdbcmd/args.rs) |

### CLI: NCBI command lines with `losat` in front

| ID | Meaning | State | Main pull requests | Record in this repository |
| --- | --- | --- | --- | --- |
| `CLI` | Line CLI: make NCBI BLAST+ command lines work with only `losat ` prefixed, and keep a ledger of 210 downstream-tool command lines frozen against NCBI 2.17.0. | done | #613, #614, #623, #628, #694, #691, #696 | [README.md](../../README.md), “Scope & Differences from NCBI BLAST+” |
| `CLI C2` | Line CLI stage C2: NCBI custom tabular fields, delim=, and -outfmt 10 and 20 for every program, through one shared tabular writer. | done | #614, #623, #628 | [README.md](../../README.md), “Output formats” |
| `CLI C2b` | Line CLI stage C2b: custom tabular fields, delim= and -outfmt 10 and 20 for TBLASTN and TBLASTX. | done | #623 | [README.md](../../README.md), “Output formats” |
| `CLI C2d` | Line CLI stage C2d: integrate the TBLASTN/TBLASTX tabular PR with the BLASTP/BLASTX one (taxonomy-name wiring, rejection rows, matrices) on the later branch. | done | #628 | [README.md](../../README.md), “Output formats” |

### DBID: BLAST databases with parsed Seq-ids

| ID | Meaning | State | Main pull requests | Record in this repository |
| --- | --- | --- | --- | --- |
| `DBID` | Work line that made all six search programs read BLAST databases built with makeblastdb -parse_seqids (versions 4 and 5) and write NCBI BLAST+ 2.17.0's subject ids and titles. | in progress | #536, #550, #572, #677, #683, #701 | [docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md](../../docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md), “Databases with parsed Seq-ids (version 1.3, line DBID)” |
| `DBID H3` | Made outfmt 0 refuse an empty stored title whose tax id the taxonomy database names, in databases without parsed Seq-ids, because NCBI generates that title. | done | #794 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md), “Databases with parsed Seq-ids (line DBID, ID2a)” |
| `DBID ID2a` | Let blastp and blastx search databases with parsed Seq-ids (versions 4 and 5, several definition lines, GIs, taxids) and write NCBI's subject ids and titles in every format. | done | #550 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md), “Databases with parsed Seq-ids (line DBID, ID2a)” |
| `DBID ID2b` | Let blastn (megablast, blastn, dc-megablast), tblastn and tblastx search databases with parsed Seq-ids and write NCBI's subject ids and titles in every format. | done | #572 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md), “Databases with parsed Seq-ids (line DBID, ID2a)” |
| `DBID ID4f` | Made blastn's outfmt 0 build its alignment rows for a window of queries at a time, as NCBI does, cutting memory on databases with parsed Seq-ids. | done | #598 | changelog entry `fix-dbid-pairwise-rows` |
| `DBID ID4t` | tblastx -db memory fixes: hit lists kept inside the batch search (#604), outfmt 0 rows built per query window (#607), reports written as each batch is taken (fix 2b, #619). | done | #604, #607, #619, #673 | changelog entry `feature-dbid-tbxstream` |
| `DBID L1` | Made a tabular report whose only alignment field is nident read a subject, and so write taxonomy messages, only on rows with identities, as NCBI does with the prefetch off. | done | #830 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md), “`nident` alone and HSPs without identities (line DBID, L1)” |
| `DBID L1b` | Made tblastn and tblastx without the report stream refuse a run for NCBI's taxonomy cache only where NCBI fetches a subject, not for fields that read none or for batches without hits. | done | #864 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md), “Rows `db.nz.cq.*` (line DBID, L1b)” |
| `DBID N1` | Made -db searches refuse a reported subject whose Seq-align id NCBI's lookup finds in another sequence first, instead of printing the wrong sequence's names. | done | #683 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md), “Databases with parsed Seq-ids (line DBID, ID2a)” |
| `DBID N1b` | Made the same record with other tax ids in databases read together refuse only where reports show the copies apart; later also for titles made from tax ids. | done | #765, #785 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md), “Databases with parsed Seq-ids (line DBID, ID2a)” |
| `DBID N2a` | Made -db searches refuse a stored gnl\|BL_ORD_ID\|n whose n is not the sequence's own number, since NCBI looks it up by that id and reads another sequence. | done | #730 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md), “Databases with parsed Seq-ids (line DBID, ID2a)” |
| `DBID N2b` | Left every kind of Seq-align id (bbs, nat, gnl\|PIG, lcl integer) to NCBI's database lookup instead of refusing by kind, and refused ids another sequence also holds. | done | #739 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md), “Databases with parsed Seq-ids (line DBID, ID2a)” |
| `DBID N2c` | Made -db searches refuse the weak Seq-id matches of NCBI's scope (accession without version, name-only id) that make NCBI write another sequence's names. | done | #744 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md), “Databases with parsed Seq-ids (line DBID, ID2a)” |
| `DBID N2d` | Made blastx refuse PRE_FETCH_SEQS_LIMIT values NCBI cannot convert and BL2SEQ_LEGACY with -subject, which NCBI reads but LOSAT ignored with exit 0. | done | #754 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md), “The environment around BLASTX's search and report (line DBID, N2d)” |
| `DBID N2e` | Made the other programs refuse PRE_FETCH_SEQS_LIMIT and BL2SEQ_LEGACY only where NCBI's app reads them; outfmt 11 and -db now write NCBI's bytes. | done | #759 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md), “The same environment in the other programs (line DBID, N2e)” |
| `DBID N3` | Made outfmt 0, 5 and 13 to 16 show copies of one record that follow each other in a hit list as one hit, as NCBI does, instead of refusing. | done | #781 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md), “Databases with parsed Seq-ids (line DBID, ID2a)” |
| `DBID P1` | Made psiblast and blastp write NCBI's taxonomy database prefetch messages, including when PRE_FETCH_SEQS_LIMIT turns the prefetch off. | done | #802 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md), “The taxonomy database's prefetch messages (line DBID, item T, review 1)” |
| `DBID P2` | Wrote the taxonomy messages in blastn, tblastn, blastx and tblastx where PRE_FETCH_SEQS_LIMIT turns the prefetch off, after each query's warnings as NCBI's formatters fetch. | done | #807 | [docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md](../../docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md), “Databases with parsed Seq-ids (version 1.3, line DBID)” |
| `DBID P3` | Wrote the taxonomy messages in outfmt 8, 9 and 12 at PRE_FETCH_SEQS_LIMIT=2147483647, where NCBI prefetches before testing the format, and refused subjects it cannot read as NCBI does. | done | #817 | [docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md](../../docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md), “Databases with parsed Seq-ids (version 1.3, line DBID)” |

### DI: the drop-in survey (items DI-n)

| ID | Meaning | State | Main pull requests | Record in this repository |
| --- | --- | --- | --- | --- |
| `DI-1` | Drop-in item DI-1: `losat blastdbcmd -info`, which ABRicate and rBLAST read (it is stage BDC1). | done | #616 | [LOSAT/tests/blastdbcmd/README.md](../../LOSAT/tests/blastdbcmd/README.md) |
| `DI-2` | Drop-in item DI-2: blastx -task blastx-fast and -query on standard input (ABRicate protein databases, Galaxy), then memory work on blastx-fast. | done | #637, #648, #684, #706, #727, #735, #792, #803, #821 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md) |
| `DI-3` | Drop-in item DI-3: blastp -comp_based_stats 0 (AMRFinderPlus, Galaxy), NCBI ordinary traceback. | done | #624 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md) |
| `DI-4` | Drop-in item DI-4: tblastn -task tblastn-fast and -word_size 5 to 7 with the compressed lookup, plus -window_size 15 -threshold 100 (AMRFinderPlus). | done | #621, #632 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md) |
| `DI-5` | Drop-in item DI-5: blastn -ungapped (mlst, ngmaster, Galaxy). | done | #629 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md) |
| `DI-6` | Drop-in item DI-6: blastn -xdrop_ungap, -xdrop_gap, -xdrop_gap_final, -window_size, -off_diagonal_range and -no_greedy (pyani ANIb). | done | #638 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md) |
| `DI-6b` | Drop-in item DI-6b: blastn -min_raw_gapped_score as NCBI 2.17.0 replaces the E-value cutoffs. | done | #655 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md) |
| `DI-7` | Drop-in item DI-7: -show_gis, -searchsp, -parse_deflines on the five programs, blastn -soft_masking, tblastx -lcase_masking (MAKER, Galaxy). | done | #634, #635, #656, #642, #680, #698, #700 | [docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md](../../docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md) |
| `DI-8` | Drop-in item DI-8: restrict a -db search to a subset with -taxids, -taxidlist, -seqidlist, their negative forms, -no_taxid_expansion, and several names in -db (BlobToolKit, MOB-suite, nf-core, Galaxy). | done | #630, #640, #654, #664, #668, #685 | [docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md](../../docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md) |
| `DI-8a` | DI-8 stage a: -taxids, -taxidlist, -negative_taxids, -negative_taxidlist and -no_taxid_expansion with -db on version 5 databases (LMDB and taxid-list readers). | done | #630 | [docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md](../../docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md) |
| `DI-8b` | DI-8 stage b: -seqidlist and -negative_seqidlist with -db (acc2oid lookup, text and binary lists, string ISAM for version 4). | done | #640 | [docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md](../../docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md) |
| `DI-8c` | DI-8 stage c: a taxid or seqid list that keeps no sequence of the database is searched as NCBI does (nothing searched, warnings, same report). | done | #654 | [docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md](../../docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md) |
| `DI-8d` | DI-8 stage d: a positive taxid or seqid list writes the use_this_seqid extension on each Seq-align (outfmt 8, 9, 11, 12). | done | #664 | [docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md](../../docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md) |
| `DI-8f` | DI-8 stage f: -db with several names ("a b") is searched as one virtual database, as NCBI CSeqDB does (Galaxy multiple-database choice). | done | #685 | [docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md](../../docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md) |
| `DI-9` | Drop-in item DI-9: blastdbcmd -entry, -entry_batch and -outfmt retrieval (nf-core, Galaxy), delivered as stages BDC2a, BDC2b and BDC2c. | done | #647, #665, #676 | [LOSAT/tests/blastdbcmd/README.md](../../LOSAT/tests/blastdbcmd/README.md) |
| `DI-12` | Drop-in item DI-12: blastp and blastx search protein databases and queries holding the gap character `-` as NCBI does (stored as a residue; dropped with a warning in a query). | done | #627 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md) |
| `DI-16` | Drop-in item DI-16: speed up TBLASTX all-versus-all (ViPTreeGen) with every matrix on the unit path and as few groups as the pool needs (LOSAT_X_TBXFEWGROUPS). | done | #652 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md) |

### GZ: compressed FASTA input

| ID | Meaning | State | Main pull requests | Record in this repository |
| --- | --- | --- | --- | --- |
| `GZ` | Line GZ: FASTA inputs named .gz (then .bz2 and .zst) are read through decompression by the CLI (-query, -subject, makeblastdb -in) and by LOSAT Web (plan row DW-45). | done | #782, #780, #791, #788, #797, #795 | [docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md](../../docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md) |

### NCMEM: maximum RSS above NCBI's

| ID | Meaning | State | Main pull requests | Record in this repository |
| --- | --- | --- | --- | --- |
| `NCMEM1` | NCMEM stage 1: blastn streams the custom tabular reports batch by batch and drops the pre-search check of reported sequences, cutting max RSS. | done | #773 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md) |
| `NCMEM2` | NCMEM stage 2: megablast reads its queries batch by batch, and the second QSTREAM reading re-reads the parts of the first in parallel, lowering max RSS. | done | #783 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md) |

### NP: NCBI features LOSAT lacked

| ID | Meaning | State | Main pull requests | Record in this repository |
| --- | --- | --- | --- | --- |
| `NP` | Line NP: give LOSAT the NCBI BLAST+ 2.17.0 features it lacked (WindowMasker, -qcov_hsp_perc, -strand, -culling_limit, -ungapped, soft masking, scoring matrices) with byte-identical output. | in progress | #329, #338, #345, #361, #364, #501, #527, #629, #656 | [docs/web/ncbi_ui_mapping.md](../../docs/web/ncbi_ui_mapping.md) |
| `NP M2` | Stage M2: blastn -window_masker_db, masking queries with a WindowMasker counts file in all four blastn tasks, byte-identical to NCBI BLAST+ 2.17.0. | done | #501, #512 | changelog entry `feature-np-m2` |
| `NP M3` | Stage M3: blastn -window_masker_taxid, finding the counts file through WINDOW_MASKER_PATH, .ncbirc or the current directory. | in progress | none | [docs/evidence/ncbi_unit_cases/LEDGER.tsv](../../docs/evidence/ncbi_unit_cases/LEDGER.tsv) |
| `NP O1b` | Stage O1b: -num_descriptions and -num_alignments (outfmt 0-4 hit counts) in all five programs, matching NCBI BLAST+ 2.17.0. | done | #345, #529 | [LOSAT/tests/outfmt0_manifest.tsv](../../LOSAT/tests/outfmt0_manifest.tsv) |
| `NP O3b` | Stage O3b: -culling_limit for blastp and tblastn, using NCBI hspfilter_culling on the culling core shared with blastn and tblastx. | done | #527 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md), “-culling_limit rows (line NP, O3b)” |
| `NP O4a` | Stage O4a: blastn -ungapped, finishing the search in the preliminary stage as NCBI does; O4a-1 covers outfmt 0, 5, 7, 13-16. | in progress | #629 | [LOSAT/tests/db_search/README.md](../../LOSAT/tests/db_search/README.md), “Ungapped nucleotide databases (line NP, O4a-1)” |
| `NP SC1` | Stage SC1: tblastn accepts all nine NCBI matrices, every gap-cost pair of their tables, recommended word thresholds and composition-based statistics modes 0 and 2. | done | #631 | [LOSAT/tests/outfmt0_manifest.tsv](../../LOSAT/tests/outfmt0_manifest.tsv) |
| `NP SC2a` | Stage SC2a: blastp accepts BLOSUM62 with every gap-cost pair of NCBI's table, with -comp_based_stats 2 and 0, -task blastp and blastp-fast, -subject and -db. | done | #824 | [LOSAT/tests/outfmt0_manifest.tsv](../../LOSAT/tests/outfmt0_manifest.tsv) |
| `NP SC4a` | Stage SC4a: blastp -comp_based_stats 0 runs NCBI ordinary traceback (no composition adjustment), also for -db. | done | #624 | [LOSAT/tests/outfmt0_manifest.tsv](../../LOSAT/tests/outfmt0_manifest.tsv) |
| `NP SC8a` | Stage SC8a: blastp -word_size 6 and 7 with NCBI compressed lookup table (BLOSUM62). | done | #822 | [LOSAT/tests/outfmt0_manifest.tsv](../../LOSAT/tests/outfmt0_manifest.tsv) |
| `NP SC11` | Stage SC11: tblastn compressed lookup table (tblastn-fast, word size 5-7) for every matrix except IDENTITY. | done | #823, #871 | [LOSAT/tests/outfmt0_manifest.tsv](../../LOSAT/tests/outfmt0_manifest.tsv) |
| `NP SM-1` | Stage SM-1: blastp -soft_masking true/false (SEG masks only the lookup table when true), matching NCBI BLAST+ 2.17.0. | done | #492 | [LOSAT/tests/outfmt0_manifest.tsv](../../LOSAT/tests/outfmt0_manifest.tsv) |
| `NP SM-2` | Stage SM-2: tblastx -soft_masking true/false, matching NCBI BLAST+ 2.17.0 including sum statistics. | done | #499 | [LOSAT/tests/outfmt0_manifest.tsv](../../LOSAT/tests/outfmt0_manifest.tsv) |
| `NP SM-3` | Stage SM-3: blastn -soft_masking, with false writing DUST, lower-case and WindowMasker masks as N into the searched query. | done | #656 | changelog entry `feat-blastn-soft-masking` |
| `NP WW0` | Stage WW0: bring threaded Wasm shared memory of 65,535 pages (earlier PRs #136, #137) from the Web branches into dev. | done | #488 | [docs/evidence/wasm_memory_max_20261012/README.md](../../docs/evidence/wasm_memory_max_20261012/README.md) |

### OF: output formats (`-outfmt`)

| ID | Meaning | State | Main pull requests | Record in this repository |
| --- | --- | --- | --- | --- |
| `OF` | Work line that made LOSAT write every NCBI BLAST+ 2.17.0 -outfmt (1 to 18 and 20, delim=, custom fields) byte for byte in all five programs, with -subject and -db. | in progress | #319, #337, #341, #346, #352, #363, #373, #385, #389, #404, #416, #463, #471, #524, #540 | [README.md](../../README.md), “Output formats” |
| `OF3` | Added outfmt 8 and 12 (Seq-annot as ASN.1 text and JSON) and then 9 (binary ASN.1) in all five programs, with the Seq-align model and NCBI's real-number digits. | done | #337, #341 | [README.md](../../README.md), “Output formats” |
| `OF3b` | Added -outfmt 8 (Seq-annot, ASN.1 text) and 12 (Seq-annot, JSON) in all five programs. | done | #337 | [README.md](../../README.md), “Output formats” |
| `OF3c` | Added -outfmt 9 (Seq-annot, binary ASN.1) in all five programs. | done | #341 | [README.md](../../README.md), “Output formats” |
| `OF4` | Added -outfmt 5 (BLAST XML) with per-query statistics, filter string and masked sequences to all five programs. | done | #346, #352 | [README.md](../../README.md), “Output formats” |
| `OF4a` | Added -outfmt 5 (BLAST XML) to blastn and blastp, with -subject and -db. | done | #346 | [README.md](../../README.md), “Output formats” |
| `OF4b` | Added -outfmt 5 (BLAST XML) to blastx, tblastn and tblastx, with every refusal before the first output byte. | done | #352 | [README.md](../../README.md), “Output formats” |
| `OF5` | Added -outfmt 15 and 16 (single-file BLAST JSON and XML2) and 13 and 14 (multi-file, with -out) to all five programs. | done | #363, #373, #385 | [README.md](../../README.md), “Output formats” |
| `OF5a` | Added -outfmt 15 and 16 (BLAST JSON and XML2) to blastp and blastn, with taxonomy database names. | done | #363 | [README.md](../../README.md), “Output formats” |
| `OF5b` | Added -outfmt 15 and 16 (BLAST JSON and XML2) to blastx, tblastn and tblastx. | done | #373 | [README.md](../../README.md), “Output formats” |
| `OF5c` | Added -outfmt 13 and 14 (multi-file BLAST JSON and XML2: per-query side files plus a master file, -out only) to all five programs. | done | #385 | [README.md](../../README.md), “Output formats” |
| `OF6` | Added -outfmt 17 (SAM) to blastn, with NCBI's CIGAR, clipping, sort order, SQ and SR options and the header written last. | done | #389 | [README.md](../../README.md), “Output formats” |
| `OF8` | Added -outfmt 11 (Blast4 archive; ASN.1 text, binary or XML via ARCHIVE_FORMAT) to all five programs, with BATCH_SIZE batching. | done | #404, #416, #463, #471, #524, #540 | [docs/product_decisions/PD-LOSAT-ARCHIVE-DEFLINE.md](../../docs/product_decisions/PD-LOSAT-ARCHIVE-DEFLINE.md) |
| `OF8a` | Added -outfmt 11 (BLAST archive) to blastn (OF8a1, #404) and blastp (OF8a2, #416). | done | #404, #416 | [README.md](../../README.md), “Output formats” |
| `OF8b` | Added -outfmt 11 (BLAST archive) to blastx with BATCH_SIZE (OF8b1), tblastn (OF8b2) and tblastx (OF8b3). | done | #463, #471, #524 | [README.md](../../README.md), “Output formats” |
| `OF8c` | Added ARCHIVE_FORMAT=xml, writing NCBI's XML form of the -outfmt 11 archive, in all five programs. | done | #540 | [README.md](../../README.md), “Output formats” |

### PS: psiblast

| ID | Meaning | State | Main pull requests | Record in this repository |
| --- | --- | --- | --- | --- |
| `PS` | Line PS: add psiblast and related tools (iterations, PSSM files, tblastn -in_pssm, rpsblast, rpstblastn, deltablast, makeprofiledb) to LOSAT. | in progress | #503, #874, #892, #906, #917 | [README.md](../../README.md), “7. PSI-BLAST (psiblast)” |
| `PS1` | PS stage 1: `losat psiblast` runs the first round (-num_iterations 1) of NCBI BLAST+ 2.17.0 psiblast with -subject or -db, for outfmt 0, 6 and 7, on the blastp engine. | done | #503 | changelog entry `feature-ps-ps1` |
| `PS2` | PS stage 2: the PSSM engine (NCBI core/blast_psi.c and blast_psi_priv.c) and psiblast's -out_ascii_pssm after the first round (-save_pssm_after_last_round, -save_each_pssm), byte for byte as NCBI BLAST+ 2.17.0. | done | #874 | changelog entry `feature-ps-ps2` |
| `PS3` | PS stage 3: psiblast's -out_pssm, the PSSM as NCBI's ASN.1 text `PssmWithParameters` after the first round, byte for byte as NCBI BLAST+ 2.17.0, and the reader of that text (for -in_pssm, PS4). | done | #892 | changelog entry `feature-ps-ps3` |
| `PS4a` | PS stage 4a: psiblast -in_pssm with -comp_based_stats 0, the search with a PSSM query (lookup table, ungapped and gapped extensions, traceback from the PSSM's rows), byte for byte as NCBI BLAST+ 2.17.0. | done | #906 | changelog entry `feature-ps-ps4` |
| `PS4b` | PS stage 4b: psiblast -in_pssm with the position-based composition-based statistics (-comp_based_stats 1, and 2 and 3 run as 1 with NCBI's warning; the default), the zero-ratio warning and PSSMs of frequency ratios only, byte for byte as NCBI BLAST+ 2.17.0. | done | #917 | changelog entry `feature-ps-ps4b` |
| `PS5a` | PS stage 5a: psiblast's rounds of FASTA queries (-num_iterations N, convergence, the two-part description table, the PSSM files of every round) and -comp_based_stats 0, with -subject or -db, for outfmt 0, 6, 7, 10 and 20, byte for byte as NCBI BLAST+ 2.17.0. | in progress | none | changelog entry `feature-ps-ps5` |
| `PS11` | PS stage 11 (series PR): `losat rpsblast`, NCBI BLAST+ 2.17.0 rpsblast against profile databases (makeprofiledb -dbtype rps); split into PS11a and PS11b. | in progress | none | changelog entry `feature-ps-ps11` |
| `PS11a` | PS stage 11a: `losat rpsblast` without composition-based statistics (-comp_based_stats 0, F, f, D, d; any mode on a database without frequency ratios, as NCBI runs it) on one-volume profile databases (version 4 and 5, 26 and 28 columns, any whole scale), outfmt 0, 6, 7 and 10, byte for byte as NCBI BLAST+ 2.17.0 at every thread count. | in progress | none | changelog entry `feature-ps-ps11` |
| `PS11b` | PS stage 11b: `losat rpsblast` with composition-based statistics 1, 2 and 3 (the default 1; `.freq` ratios, the PSSM redo), databases of several volumes (merged per query as NCBI does), `-dbsize` and the 32-bit score bound, byte for byte as NCBI BLAST+ 2.17.0 at every thread count. | in progress | none | changelog entry `feature-ps-ps11b` |

### RL: the v0.2.0 release

| ID | Meaning | State | Main pull requests | Record in this repository |
| --- | --- | --- | --- | --- |
| `RL1` | Release line session RL1 of v0.2.0: lineage and certification set, release-candidate contract and assembler, readiness workflow, first switch-graduation wave, documentation fixes. | done | #505, #508, #510, #506, #507, #511 | [docs/release/v0.2.0_certification_plan.md](../../docs/release/v0.2.0_certification_plan.md) |

### SF: structural work

| ID | Meaning | State | Main pull requests | Record in this repository |
| --- | --- | --- | --- | --- |
| `SF` | Structural line SF (Oct-Nov 2026): flatten fixed costs (S1), unify the parallel drivers (S2), fold conditional speed switches (S3). One test comment instead means the earlier FASTA-reader port stage. | in progress | #653, #663, #682, #699, #709, #720, #724, #737, #746, #756, #772, #786, #804, #748 (closed) | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md) |
| `SF S1` | SF stage S1: flatten per-unit fixed costs (allocation, zeroing, re-reading, serial head and tail) so the ratio to NCBI stays within x1.15 as the number of work units grows tenfold. | in progress | #653, #724, #737, #746, #772, #786, #804 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md) |
| `SF S3` | SF stage S3: fold conditional speed switches family by family into one form, aiming at 68 or fewer conditional switches and no timing-midpoint thresholds. | in progress | #663, #682, #709, #720, #756, #748 (closed) | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |

### Plan v4 (speed and memory lines)

| ID | Meaning | State | Main pull requests | Record in this repository |
| --- | --- | --- | --- | --- |
| `v4-A` | Plan v4 line A: lane-wise seed stage and ungapped extension, tblastx first, then blastp, blastx and tblastn (prefetch, byte-table walks, extension batches, query slices). | done | #300, #301, #307, #320, #325, #474, #548, #612 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-A10` | Item A10: branch-free two-hit hit loop (too-far/overlap test as a value select), switches TBXHITBF, TBNHITBF, BPHITBF, BXHITBF; later sessions also did blastn memory and fill items. | done | #318, #322, #344, #353 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-A14` | Item A14: blastx cell test for short reads against protein databases (deferred calls, cross-cell batching, masked batches). | done | #394, #433, #444, #447, #528 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-A15` | Item A15: dc-megablast chains as CSR tables and the three-step fold core, plus follow-up seed-stage list scans. | done | #366, #441, #594, #606 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-B` | Plan v4 line B: linking and thread scaling (tblastx linking merges, report stages on the pool, blastn -db result hand-over, X-drop traceback memory). | done | #299, #303, #312, #314, #328, #349, #357, #390, #402, #412 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-B10` | Item B10: run a long blastx query chunks one after another through the subject flow, then tblastn -db flow, redo and memory items (B11). | done | #393, #583, #421, #549, #586, #611 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-BPF` | Plan v4 line BPF: make -task blastp-fast faster under today bytes-identical rules, measured against NCBI fast and DIAMOND (seed loop, table build, query-block slots). | done | #472, #532, #543, #551, #597, #579, #589 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-C` | Plan v4 line C: blastp composition-based redo (Newton and log cost, traceback DP, SEG trim memo), plus a greedy-row fusion for megablast. | done | #308, #315, #326, #477, #552 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-D` | Plan v4 line D: fixed cost, build and database costs for blastn (lookup-table keeping, QBATCHPAR slots, CSR heads, report formatting on the pool), later memory only. | done | #306, #311, #317, #330, #438, #475, #514, #539, #580, #617 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-D12` | Item D12: bring blastp -db anonymous memory (UniRef50 scale) below NCBI, then the same for tblastn, blastx and blastn -db (streamed reads, pools, block sizes). | done | #343, #351, #359, #443, #576, #601, #697, #743 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-E` | Plan v4 line E (measurement only): unmeasured use shapes, namely short-read blastx against large protein databases and two-socket NUMA sensitivity. | done | none | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-F` | Plan v4 line F (measurement only): the exact floor, i.e. work per stage times cost per unit, for each program, to rank remaining speed-up items. | done | none | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-F P3` | Line F probe P3: tblastn lookup scan with a 16-word AVX2 index and a one-byte presence table (V2), later A16. | done | #362, #370 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-G` | Plan v4 line G (gate-gap audit): find input shapes and options where an existing default speed-up silently did not run, and widen it (-subject packing, block checks, few-subject searches). | done | #325, #334, #347, #377, #391, #397, #410, #437, #479, #525 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-MEM` | Plan v4 memory line: audit max RSS in every case against NCBI (required: at or below NCBI at n <= 16) and fix the cells above it. | in progress | #340, #544, #545, #608, #625, #633, #697, #743, #753, #792, #803, #821 | [docs/release/v0.2.0_certification_plan.md](../../docs/release/v0.2.0_certification_plan.md) |
| `v4-MF` | Plan v4 line MF (measurement only): memory floor, the bytes each case needs while keeping output and speed, per component, from a tagging allocator. | done | none | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-R` | Plan v4 line R: root-cause counters without code changes (X-drop DP kernel, blastp cell test, CPU growth at 16 threads, megablast fill, lookup CSR). | done | none | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-R6` | Report R6: instruction budget of the tblastx x_hit/scan_unit loop (counts only), showing a branch-free form would run about x0.8. | done | #318 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-R R1` | Report R1: per-row time of the X-drop DP kernels align8 and align16, limited by the row-to-row dependency chain rather than ports. | done | #308 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-R R10` | Report R10: the one-thread megablast lookup fill, 0.190 s on one case, about 87 instructions per word, compared with NCBI. | done | #408 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default)” |
| `v4-W` | Plan v4 line W: Wasm and browser speed (tier-up avoidance, one worker pool per run, warm-up, startup cost, simd128 kernels, memory of thread workers). | done | #316, #323, #348, #358, #679, #688, #722 | external |
| `v4-W W7` | Item W7: close the warm Wasm versus native gap function by function (byte-table two-hit walks and the hit loop in two parts on wasm32, report rows). | done | #323, #348, #358 | [docs/agents/diagnostics.md](../../docs/agents/diagnostics.md), “Performance switches (on by default) (wasm32 note under TBXUNGAP)” |

## LOSAT Web

### SX: BLASTX in LOSAT Web

| ID | Meaning | State | Main pull requests | Record in this repository |
| --- | --- | --- | --- | --- |
| `SX` | Line SX: make BLASTX usable in LOSAT Web (engine entry, ABI adapter, app and browser test cells) and give BLASTX -query_loc and -subject_loc in the CLI and the Web. | done | #731, #733, #734, #745, #760, #749, #736 | [docs/losat_web_gui_plan.md](../../docs/losat_web_gui_plan.md) |
| `SXa` | Stage SXa of line SX: BLASTX gets run_local, the engine entry LOSAT Web uses, with the query batch loop shared with the CLI path, plus native checks (V-NAT). | done | #731 | [docs/web/abi_v2.md](../../docs/web/abi_v2.md), “1. Purpose” |
| `SXb` | Stage SXb of line SX: the ABI v2 adapter runs BLASTX (run, validate, count, describe) and the V-ABI suite covers it. | done | #733 | [docs/web/abi_v2.md](../../docs/web/abi_v2.md), “1. Purpose” |
| `SXc` | Stage SXc of line SX: the Web app offers BLASTX with formats 0, 6 and 7, its option lists and V-BR browser cells; badge words became Supported. | done | #734 | [docs/losat_web_gui_plan.md](../../docs/losat_web_gui_plan.md) |
| `SXd` | Stage SXd of line SX: port NCBI's -query_loc and -subject_loc to BLASTX in three steps: SXd1 CLI formats 0/6/7/10, SXd2 other formats, SXd3 the Web. | done | #745, #760, #749 | [docs/losat_web_gui_plan.md](../../docs/losat_web_gui_plan.md) |
| `SXd1` | Step SXd1 of line SX: BLASTX reads -query_loc and -subject_loc on the command line for outfmt 0, 6, 7 and 10, matched to NCBI 2.17.0 frozen fixtures. | done | #745 | [docs/losat_web_gui_plan.md](../../docs/losat_web_gui_plan.md) |
| `SXd2` | Step SXd2 of line SX: BLASTX ranges also in outfmt 5, 8, 9 and 11 to 16, with the split mask and lower-case handling under -subject_loc. | done | #760 | [docs/evidence/losat_web_e2d/README.md](../../docs/evidence/losat_web_e2d/README.md) |
| `SXd3` | Step SXd3 of line SX: LOSAT Web shows the range fields for BLASTX and runs BLASTX ranges as the CLI does. | done | #749 | [docs/losat_web_gui_plan.md](../../docs/losat_web_gui_plan.md) |

### WA: the LOSAT Web audit

| ID | Meaning | State | Main pull requests | Record in this repository |
| --- | --- | --- | --- | --- |
| `WA` | Line WA: a browser audit of LOSAT Web at one dev commit (input, form, run, results), a second audit round, and pull requests that fix the findings. | in progress | #799, #800, #805, #806, #808, #809, #810, #811, #812, #816, #818, #820, #831, #867 | [docs/web/ncbi_ui_mapping.md](../../docs/web/ncbi_ui_mapping.md), “4. 文言の方針（stop-describing-ui、2026-11-02）” |
| `WA2-BR-01` | Second-audit finding (P2): when stored results were removed while open, Outputs kept the previous format's text and Export did nothing; it now says why the results cannot be read. | done | #818 | [docs/web/ncbi_ui_mapping.md](../../docs/web/ncbi_ui_mapping.md), “4. 文言の方針（stop-describing-ui、2026-11-02）” |
| `WA2-OP-03` | Second-audit finding: -outfmt with column fields in Additional options was refused though the CLI takes it; Web now writes outfmt 6 or 7 with fields as an extra output. | done | #831 | [docs/web/ncbi_ui_mapping.md](../../docs/web/ncbi_ui_mapping.md), “4. 文言の方針（stop-describing-ui、2026-11-02）” |
| `WA2-SC-05` | Second-audit finding (P2): after a tab hidden during a run returned, the notice said results were stored but the run said they were removed; the notice now checks the store. | done | #820 | [docs/web/ncbi_ui_mapping.md](../../docs/web/ncbi_ui_mapping.md), “4. 文言の方針（stop-describing-ui、2026-11-02）” |
| `WA W1` | Stage W1 of line WA: remove the word verification from Web screens and READMEs, mark results only Supported or Not yet supported, and fix the default-run marks (WA-RS-03). | done | #805, #806, #808 | [docs/evidence/losat_web_wa_w1/README.md](../../docs/evidence/losat_web_wa_w1/README.md) |
| `WA-FM-01` | Audit finding: the Web validate refused BLASTN match/mismatch scores without a Karlin-Altschul table row, though CLI and NCBI only warn; validate now leaves the tables to the run. | done | #799 | [docs/web/abi_v2.md](../../docs/web/abi_v2.md), “4. Exports” |
| `WA-FM-02` | Audit finding: search form selects listed values that LOSAT's engine always refuses; they now offer only the NCBI values the engine reads, with Other for the rest. | done | #812 | [docs/web/ncbi_ui_mapping.md](../../docs/web/ncbi_ui_mapping.md), “4. 文言の方針（stop-describing-ui、2026-11-02）” |
| `WA-FM-03` | Audit finding: dc-megablast Template length and Template type could be set one without the other, which BLAST+ refuses; changing one now writes both. | done | #812 | [docs/web/ncbi_ui_mapping.md](../../docs/web/ncbi_ui_mapping.md), “1. 検索画面” |
| `WA-FM-04` | Audit finding: -h or -help in Additional options printed the CLI's whole help; it is now refused with a short message that help is not available in the browser. | done | #812 | [docs/web/ncbi_ui_mapping.md](../../docs/web/ncbi_ui_mapping.md), “4. 文言の方針（stop-describing-ui、2026-11-02）” |
| `WA-FM-05` | Audit finding: options meaningless in a browser (-db, -mt_mode, -subject, -num_threads) gave clap's messages; they are now refused with short messages naming the limit. | done | #812 | [docs/web/ncbi_ui_mapping.md](../../docs/web/ncbi_ui_mapping.md), “4. 文言の方針（stop-describing-ui、2026-11-02）” |
| `WA-FM-07` | Audit finding: a negative value such as -3 after -penalty looked like an option name in the shown command; the command text now reads the grammar and shows a value. | done | #812 | [docs/web/ncbi_ui_mapping.md](../../docs/web/ncbi_ui_mapping.md), “1. 検索画面” |
| `WA-IN-01` | Audit finding (P2): switching through a mode that reads a source as the other kind (BLASTN to BLASTP and back) silently dropped excluded records and subranges; now kept per kind. | done | #800 | external |
| `WA-IN-02` | Audit finding (P2): a run with an empty query completed with the misleading note No queries match the query filters; the notice now says the query is empty. | done | #800 | external |
| `WA-IN-03` | Audit finding (P2): after OPFS storage was cleared outside the page, every Run failed with a raw DOMException text; the directory is now recreated once, else a plain message shows. | done | #800 | external |
| `WA-IN-04` | Audit finding: a GenBank, EMBL or FASTQ file passed as FASTA text without comment; the source card now warns that it is read as FASTA, as the CLI does. | done | #809 | [docs/web/ncbi_ui_mapping.md](../../docs/web/ncbi_ui_mapping.md), “4. 文言の方針（stop-describing-ui、2026-11-02）” |
| `WA-IN-05` | Audit finding: inputs with a UTF-8 byte order mark or in UTF-16 were refused without explanation; the source card now names the cause. | done | #809 | [docs/web/ncbi_ui_mapping.md](../../docs/web/ncbi_ui_mapping.md), “4. 文言の方針（stop-describing-ui、2026-11-02）” |
| `WA-IN-06` | Audit finding: reloading the page discards all inputs and resets the mode (by design); the page now asks before a reload loses inputs. | done | #809 | [docs/web/ncbi_ui_mapping.md](../../docs/web/ncbi_ui_mapping.md), “4. 文言の方針（stop-describing-ui、2026-11-02）” |
| `WA-IN-07` | Audit finding: a reversed (300 to 100) or one-character subrange showed no message beside the field; it now does, and a later PR offers to swap a reversed range. | done | #809, #816 | [docs/web/ncbi_ui_mapping.md](../../docs/web/ncbi_ui_mapping.md), “4. 文言の方針（stop-describing-ui、2026-11-02）” |
| `WA-RN-01` | Audit finding: Cancel pressed just after a run entered Preparing had no effect until the search ended. | done | #810 | external |
| `WA-RN-02` | Audit finding: a failed run's message ended with two periods because the results code always appended one. | done | #800 | external |
| `WA-RN-04` | Audit finding: when the browser refused a wake lock, the page re-requested it at every run state change and the failure message flickered. | done | #810 | external |
| `WA-RS-01` | Audit finding: the filter E value at most 0 hid every HSP displayed as 0.0; the E value and bit score filters now compare the values the table shows. | done | #810 | [docs/web/results_columns.md](../../docs/web/results_columns.md) |
| `WA-RS-02` | Audit finding: the Dot Plot popup covered other HSP lines, and clicking there jumped to Hits; the popup no longer takes the clicks of lines under it. | done | #810 | [docs/web/ncbi_ui_mapping.md](../../docs/web/ncbi_ui_mapping.md), “Dot Plot” |
| `WA-RS-03` | Audit finding: default BLASTN and BLASTX runs (outfmt 7) and TBLASTN at 8 threads were marked Not supported; test cells now cover those defaults at 1 to 32 threads. | done | #806 | [docs/evidence/losat_web_wa_w1/README.md](../../docs/evidence/losat_web_wa_w1/README.md) |
| `WA-RS-04` | Audit finding: the Dot Plot tab was blank for a run with no hits; it now shows the no-hit state. | done | #810 | [docs/web/results_columns.md](../../docs/web/results_columns.md) |
| `WA-RS-05` | Audit finding: Download All only switched to the Outputs tab, and file names omitted the task so the four BLASTN tasks shared a name; both were fixed. | done | #810 | [docs/web/ncbi_ui_mapping.md](../../docs/web/ncbi_ui_mapping.md), “2. 結果画面” |

### WF: fetching accessions and Example data

| ID | Meaning | State | Main pull requests | Record in this repository |
| --- | --- | --- | --- | --- |
| `WF` | Line WF: LOSAT Web fetches accessions and per-mode Example data sets from public sites (NCBI, EBI, CNCB) only when the user asks, and saves them in the browser. | done | #761, #764, #766, #768, #777, #771 (closed) | [docs/web/fetch_design.md](../../docs/web/fetch_design.md) |
| `WF-0` | Step WF-0 of line WF: the design memo, the rewritten no-network promise and AGENTS rule 4, and the CSP connect-src for the fetch sites; no fetch yet. | done | #761 | [docs/web/fetch_design.md](../../docs/web/fetch_design.md), “6. 進め方（PR の切り方）” |
| `WF-1` | Step WF-1 of line WF: accession numbers in the query or subject box are fetched from NCBI, EBI or CNCB when the user presses Fetch. | done | #764 | [docs/web/fetch_design.md](../../docs/web/fetch_design.md), “8. WF-1 の実装（2026-11-03、Owner-delegated）” |
| `WF-2` | Step WF-2 of line WF: an Example line per mode, with data sets downloaded only on request and cached in the browser's OPFS. | done | #766 | [docs/web/fetch_design.md](../../docs/web/fetch_design.md), “9. WF-2 の実装（2026-11-03、Owner-delegated）” |
| `WF-3` | Step WF-3 of line WF: Fetch from Fastest races NCBI, EBI and CNCB, and an Example's saved copy can be written to a folder the user picks. | done | #768 | [docs/web/fetch_design.md](../../docs/web/fetch_design.md), “10. WF-3 の実装（2026-11-03、Owner-delegated）” |
| `WF-5` | Step WF-5 of line WF: a box of accessions that no one site serves gets one Fetch button per group, each sending its group to the one site it names. | done | #777 | [docs/web/fetch_design.md](../../docs/web/fetch_design.md), “11. WF-5 の実装：どこも全部は扱えない欄の組ごとのボタン（2026-11-03）” |

### WG, WO: threaded Wasm in browsers

| ID | Meaning | State | Main pull requests | Record in this repository |
| --- | --- | --- | --- | --- |
| `WG` | Line WG: guard the threaded Wasm artifacts for browsers at build time against bulk-memory traps after another thread grew shared memory, marking a guarded module with a custom section. | done | #693 | [docs/evidence/wasm_browser_guard_20261101/README.md](../../docs/evidence/wasm_browser_guard_20261101/README.md) |
| `WO` | Line WO: investigate threaded Wasm out-of-bounds traps (V8 compares atomics and ARM64 plain stores with a stale memory size) and fix them by widening the guard. | done | #703, #710, #713, #716 | [docs/evidence/wasm_oob_20261101/README.md](../../docs/evidence/wasm_oob_20261101/README.md) |

### WR: the LOSAT Web refresh

| ID | Meaning | State | Main pull requests | Record in this repository |
| --- | --- | --- | --- | --- |
| `WR` | Line WR: bring LOSAT Web in line with the current LOSAT and NCBI's classic pages (results page without a queue, mode wording, new look, engine-owned defaults, Additional options). | in progress | #717, #718, #719, #725, #726, #732, #738, #742, #747, #755, #757, #763 | [docs/losat_web_gui_plan.md](../../docs/losat_web_gui_plan.md) |
