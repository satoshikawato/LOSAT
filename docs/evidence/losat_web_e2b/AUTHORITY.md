# TBLASTX `-outfmt 0` and `-outfmt 7` authority (Session S08, stage E2b)

This record fixes what TBLASTX `-outfmt 0` (pairwise report) and `-outfmt 7` (tabular with comments) reproduce, and where LOSAT does it. Section A traces NCBI's application path, section B what the report receives from the search, section C the pairwise report, section D the tabular formats, section E the display translation and the genetic codes, section F the decisions of this session, and section G the inventory of NCBI's path and what the port did with each row. The gate record is [README.md](README.md).

| Item | Pinned value |
| --- | --- |
| NCBI source (sole authority) | `/mnt/c/Users/genom/GitHub/ncbi-blast/`, commit `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`. Paths below are relative to `c++/`. Line numbers count lines after removing CR (`tr -d '\r' \| cat -n`). |
| NCBI BLAST+ (comparison oracle only) | `tblastx`, `tblastn`, `blastn` and `makeblastdb` 2.17.0+ (package `blast 2.17.0`, build `Aug 11 2025 09:46:06`) in `/home/kawato/micromamba/bin/`. No `.ncbirc`; `BL2SEQ_LEGACY`, `CTOOLKIT_COMPATIBLE`, `OLD_FSC`, `BATCH_SIZE` and `CHUNK_SIZE` unset unless a case sets one (`run_oracle.py` and the fixture scripts refuse otherwise). |
| LOSAT | Branch `feature/losat-web-gui`; the session's engine commits are `fd896a52c` (port) and `0d533ba76` (regression fixtures and CI). Line numbers below are those of `0d533ba76`. |
| Fixtures | `tblastx.*` rows of [`LOSAT/tests/outfmt0_manifest.tsv`](../../../LOSAT/tests/outfmt0_manifest.tsv) (26 rows, frozen NCBI output in `LOSAT/tests/fixtures/outfmt0/`, inputs made by [`make_inputs.py`](make_inputs.py)); 40 cases of [`LOSAT/tests/tblastx_regression_fixtures.py`](../../../LOSAT/tests/tblastx_regression_fixtures.py) (frozen in `LOSAT/tests/fixtures/tblastx_regression/`). |
| Inventory | [`INVENTORY.tsv`](INVENTORY.tsv): 422 rows over 7 ranges (A application and orchestration, B alignments, C description table, D outfmt 7/6, E what the report receives from the search, F residues and genetic codes, G the hit list), each with the state before the port and the result of the port (§G). |

---

## A. NCBI application path (`tblastx -query Q -subject S`, local)

| Step | NCBI file:line | Snippet | LOSAT |
|---|---|---|---|
| Subject adapter before the query | `src/app/blast/tblastx_app.cpp:117-118` | `InitializeSubject(db_args, opts_hndl, m_CmdLineArgs->ExecuteRemotely(),` | `algorithm/tblastx/blast_engine/run_impl.rs:644` (`run`: subject read, then `check_subjects_not_empty` at 740) |
| Empty query file | `src/app/blast/tblastx_app.cpp:132-135` | `if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())){` / `ERR_POST(Warning << "Query is Empty!");` | `run_impl.rs:775` (`search_cli`; `blastn::input::is_blank`) |
| Query batches | `src/app/blast/tblastx_app.cpp:136-137` | `CBlastInput input(&fasta, m_CmdLineArgs->GetQueryBatchSize());` | `run_impl.rs:1233` (`tblastx_query_batch_size`), 1271 (`run_in_pool`) |
| `BATCH_SIZE` | `src/algo/blast/blastinput/blast_input_aux.cpp:86-88` | `char* batch_sz_str = getenv("BATCH_SIZE");` / `retval = NStr::StringToInt(batch_sz_str);` | `blastinput/query_batch.rs` (`get_query_batch_size`); a value that is not an integer is rejected (§F) |
| Default size | `src/algo/blast/blastinput/blast_input_aux.cpp:130-134` | `case eTblastx:` … `retval = 10002;` | `run_impl.rs:1233` |
| Grouping | `src/algo/blast/blastinput/blast_input.cpp:138-171` | `while (size_read < GetBatchSize()) {` | `blastinput/query_batch.rs` (`next_query_batch_end`) |
| Prolog before the first batch | `src/app/blast/tblastx_app.cpp:173` | `formatter.PrintProlog();` | `run_impl.rs:994` (`write_pairwise_prologs`); a `-out` file keeps the report order (`prolog_before_search`, 974) |
| A batch of size 0 | `src/algo/blast/api/objmgr_query_data.cpp:375-380` | `NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");` | `run_impl.rs:1271` (after the prolog; exit 3) |
| One search per batch | `src/app/blast/tblastx_app.cpp:176-195` | `for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {` / `CLocalBlast lcl_blast(queries, opts_hndl, db_adapter);` | `run_impl.rs:1373` (`search_query_batch`) |
| A batch without a valid context is not searched | `src/algo/blast/api/local_blast.cpp:166-177` | `int status = m_PrelimSearch->CheckInternalData();` | `search_query_batch` (`searched: false`) |
| Its warning (translated programs) | `src/algo/blast/core/blast_stat.c:2815-2823` | `if (valid_context == FALSE) {` | `run_impl.rs:1036` (`write_tblastx_outputs`, per query of the batch, through `report::query_warnings`) |
| Reports per query | `src/app/blast/tblastx_app.cpp:203-205` | `formatter.PrintOneResultSet(**result, query_batch);` | `report.rs` and `report/pairwise.rs` (§C, §D) |
| Epilog | `src/app/blast/tblastx_app.cpp:209` | `formatter.PrintEpilog(opt);` | `report/pairwise.rs:2664` (`write_tblastx_epilog`), `algorithm/tblastx/report.rs:544` (outfmt 7 footer) |
| Write failure (outfmt 0) | `src/app/blast/blast_app_util.hpp:252-255` | `catch (const std::ios::failure&) {` | `run_impl.rs:644` (`crate::cli::ReportStream`; "BLAST failed to write output", exit 6) |
| Descriptions and alignments | `src/algo/blast/blastinput/blast_args.cpp:2910-2928`; `src/objtools/align_format/format_flags.cpp:219,221` | `m_NumDescriptions = args[kArgMaxTargetSequences].AsInteger();`; `kDfltArgNumDescriptions = 500;`, `kDfltArgNumAlignments = 250;` | `algorithm/tblastx/args.rs` (`max_target_seqs: Option<usize>`), `run_impl.rs:1036` |
| Fewer than 5 matches | `src/algo/blast/blastinput/blast_args.cpp:2975-2977` | `ERR_POST(Warning << "Examining 5 or more matches is recommended");` | `run_impl.rs:775` (`report::query_warnings::few_matches_warning`) |

The order of NCBI's checks and errors, as LOSAT's `run` follows it: options (`-outfmt`, `-max_target_seqs`), the subject file (an empty subject: `BLAST engine error: Empty CBlastQueryVector`, exit 3), the query file, `-out`, "Query is Empty!", then `BATCH_SIZE` and the batches. Warnings go through `QueryWarnings`, which flushes the report stream first (NCBI posts them on `cerr`, tied to `cout`), so a merged stream (`2>&1`) has NCBI's order (fixtures `*.merged`).

## B. What the report receives from the search

| Item | NCBI file:line | Snippet | LOSAT |
|---|---|---|---|
| Hit list per query | `src/algo/blast/core/blast_hits.c:3243-3300` | `Blast_HitListUpdate(BlastHitList* hit_list,` | `algorithm/tblastx/report.rs:765` (`final_hit_order`, with `blastn::hsp::HitList`) |
| Subject order | `src/algo/blast/core/blast_hits.c:3071-3107` | `s_EvalueCompareHSPLists(const void* v1, const void* v2)` | `common.rs` (`evalue_compare_hsps`); `report.rs:765` |
| The list of a subject is in score order until the hit list overflows, then in e-value order | `src/algo/blast/core/blast_hits.c:3269-3285` | `Blast_HSPListSortByEvalue(hsp_list);` | `report.rs:765` |
| Preliminary search on ncbi2na with random bases for ambiguity letters | `src/algo/blast/core/blast_engine.c:772-775`; `src/algo/blast/api/blast_objmgr_tools.cpp:427-431` | `BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,`; `CRandom random(base_length);` | `run_impl.rs:546` (`preliminary_subject_bases`); re-evaluation on the ncbi4na translation |
| Sum statistics `num` | `src/algo/blast/core/link_hsps.c:1775-1776`, `1044-1061` | `hsp_list->hsp_array[index]->num = 1;` | `algorithm/tblastx/sum_stats_linking/linking.rs:557`, 2121, 847 (`replay_ncbi_output_list`) |
| Query footer: first valid context | `src/algo/blast/api/blast_results.cpp:82-100` | `// find the first valid context corresponding to this query` | `run_impl.rs:1373` (`TblastxQueryStats`) |
| The init hit list of a subject chunk, sorted before the HSP list is made (ties in different query frames by the absolute query offset) | `src/algo/blast/core/aa_ungapped.c:234-235`; `src/algo/blast/core/blast_extend.c:306-310` | `Blast_InitHitListSortByScore(init_hitlist);`; `qsort(init_hitlist->init_hsp_array, init_hitlist->total,` | `algorithm/tblastx/blast_gapalign.rs` (`sort_init_hsps_by_score_ncbi`), called before `get_ungapped_hsp_list` |
| The HSP list sorted by score after the first `BLAST_LinkHsps`, before the re-evaluation | `src/algo/blast/core/link_hsps.c:1802-1803` | `Blast_HSPListSortByScore(hsp_list);` | `run_impl.rs` (after the preliminary linking) |
| Linking groups: query and strand (six contexts per query) | `src/algo/blast/core/blast_query_info.c:76-77` | `retval->last_context = retval->num_queries * kNumContexts - 1;` | `sum_stats_linking/linking.rs` (`translated_context_group`) |
| Sum statistics of large gaps, in NCBI's order of evaluation | `src/algo/blast/core/blast_stat.c:4560-4561` | `xsum -= num*log(lcl_subject_length*lcl_query_length)` | `stats/sum_statistics.rs` (`ncbi_large_gap_sum_e`) |
| SEG keeps only the head of the left recursion's segments | `src/algo/blast/core/blast_seg.c:2096-2100` | `leftsegs->next = *segs;` | `utils/seg.rs` (`seg_seq`; BLASTX keeps the former behavior until SX) |
| SEG masks back to nucleotides | `src/algo/blast/core/blast_filter.c:892,925-931` | `Int2 BlastMaskLocProteinToDNA(BlastMaskLoc* mask_loc,`; `to = dna_length - CODON_LENGTH*seq_range->left + frame;` | `algorithm/tblastx/report.rs:283` (`query_dna_masks`, with `blastx::query_setup`) |

Sum statistics give every HSP of a chain the chain's count (`hsp_link.num` of the head); NCBI prints it as `Expect(n)` above 1 and in the description table's `N` column (the first HSP of the subject in Seq-align order).

## C. The pairwise report (`-outfmt 0`)

| Part | NCBI file:line | Snippet | LOSAT |
|---|---|---|---|
| Prolog: version, reference, database | `src/algo/blast/format/blast_format.cpp:386-441` | `CBlastFormatUtil::BlastPrintVersionInfo(m_Program, m_IsHTML, m_Outfile);` | `report/pairwise.rs:2368` (`write_tblastx_pairwise_prolog`) |
| Query preamble, no hits | `src/algo/blast/format/blast_format.cpp:1478-1515` | `<< "***** " << CBlastFormatUtil::kNoHitsFound << " *****" << "\n"` | `report/pairwise.rs:2717` (`write_tblastx_pairwise_report`) |
| Description table with `N` | `src/algo/blast/format/blast_format.cpp:497-527`; `src/objtools/align_format/showdefline.cpp:1055-1062,1100-1124,1133-1158` | `if (m_ShowLinkedSetSize) flags \|= CShowBlastDefline::eShowSumN;` | `report/pairwise.rs:1783` (`write_blastn_description_table`, `show_sum_n` true): the highest bit score (strict `>`) and its e-value, the last row's total width (E2a §G.3) |
| Alignments limited to the subjects | `src/algo/blast/format/blast_format.cpp:1545-1561` | `CBlastFormatUtil::PruneSeqalign(*aln_set, copy_aln_set, m_NumAlignments);` | `write_tblastx_pairwise_report` (500 descriptions and 250 alignments without `-max_target_seqs`, N and N with it) |
| Protein rows, middle line of characters, genetic codes | `src/algo/blast/format/blast_format.cpp:1570-1585`; `src/objtools/align_format/showalign.cpp:1949-1953` | `display.SetMasterGeneticCode(m_QueryGenCode);` / `display.SetSlaveGeneticCode(m_DbGenCode);` | `algorithm/tblastx/report.rs:137` (`displayed_rows`), `report/pairwise.rs:2514` (`write_tblastx_alignment`) |
| Score, Expect(n), Identities, Positives, Gaps, Frame | `src/objtools/align_format/showalign.cpp:310-332,3593-3603` | ` Frame = ` | `report/pairwise.rs:2422` (`write_tblastx_hsp_info`) |
| Identities and positives of the shown letters | `src/objtools/align_format/showalign.cpp:2120-2158` | — | `algorithm/tblastx/report.rs:182` (`row_counts`), 208 (`set_displayed_identities`, also for outfmt 6/7) |
| Lowercase query residues of the SEG masks | `src/objtools/align_format/showalign.cpp:2485-2574` | — | `algorithm/tblastx/report.rs:374` (`lowercase_query_row`) |
| Query footer (ungapped block only; none and search space 0 for an invalid query) | `src/algo/blast/format/blast_format.cpp:444-478` | `m_Outfile << "Effective search space used: " << summary.GetSearchSpace() << "\n";` | `report/pairwise.rs:2612` (`write_tblastx_query_footer`); an unsearched query: the `-1.00` blocks and `Gapped` |
| Epilog: matrix, threshold, window, no gap penalties | `src/algo/blast/format/blast_format.cpp:2249-2288` | `m_Outfile << "Neighboring words threshold: " << options.GetWordThreshold() << "\n";` | `report/pairwise.rs:2664` (`write_tblastx_epilog`) |
| `CTOOLKIT_COMPATIBLE` | `src/objtools/align_format/showdefline.cpp:80` | — | shared with BLASTN (E2g T11); TBLASTX in `ctoolkit_compare.py` |

## D. The tabular formats (`-outfmt 7`, `-outfmt 6`)

| Part | NCBI file:line | Snippet | LOSAT |
|---|---|---|---|
| Header per query | `src/algo/blast/format/blast_format.cpp:794-809`; `src/objtools/align_format/tabular.cpp:1266-1320` | `string strProgVersion = NStr::ToUpper(m_Program) + " " + blast::CBlastVersion().Print();` | `algorithm/tblastx/report.rs:507` (`write_outfmt7_query_header`) |
| `# Fields:` only with hits; `# N hits found` only for a searched query | `src/objtools/align_format/tabular.cpp:1266-1284` | `if (align_set) {` | `report.rs:507`, `run_impl.rs:1036` |
| Footer | `src/algo/blast/format/blast_format.cpp:2233-2238`; `src/objtools/align_format/tabular.cpp:1322-1325` | `m_Ostream << "# BLAST processed " << num_queries << " queries\n";` | `report.rs:544` (`write_tabular`) |
| Rows: pident and mismatch from the shown letters (tblastx is not NoFetch) | `src/algo/blast/format/blast_format.cpp:791-792` | — | `report.rs:208` |
| Rows limited to the hit list | `src/algo/blast/format/blast_format.cpp:811-833` | `CBlastFormatUtil::PruneSeqalign(*aln_set, copy_aln_set, m_HitlistSize);` | `report.rs:765` |

## E. The display translation and the genetic codes

NCBI translates the shown rows again from the nucleotides with `CTrans_table` (`src/objects/seqfeat/Genetic_code_table.cpp:139-232`): an ambiguous codon gives the amino acid when every resolution agrees, B for {D, N}, Z for {E, Q}, J for {I, L}, else X; the query row uses `-query_gencode`, the subject row `-db_gencode`; minus rows are the reverse complement (not flipped). The search translation stays different (X for ambiguity). LOSAT: `algorithm/tblastx/report.rs:56` (`display_base`: X/x to N, U to T), 104 (`display_residue`, with `blastx::report::display_codon`).

**Approved exception (AGENTS.md).** NCBI's local `-subject` search translates the subject with code 1 whatever `-db_gencode` says; LOSAT searches it with `-db_gencode`. The rows `tblastx.code4.0` and `tblastx.code4.7` classify the difference: `run_oracle.py` also runs NCBI with the subject as a BLAST database (`makeblastdb -dbtype nucl` without `-parse_seqids`, then `-db`), whose search applies `-db_gencode`, and `check_losat.py` compares LOSAT with that output outside the database lines (`normalize_database_lines`: the `Database:` blocks, the posted date, the `>` heading space, `# Database:`). Both rows are equal there; NCBI's local output differs in the HSPs only.

## F. Decisions of this session

Each was taken as the recommended option (the maintainer's standing instruction) unless marked.

1. **TBLASTN uses the description table of E2a §G.3 and NCBI's nucleotide titles** (`ncbi_nucleotide_title`). Every frozen TBLASTN output is unchanged; two TBLASTN repros of inventory range C (p417 and p716, the first HSP not the best) now match NCBI. BLASTP stays on its table until S08+ (protein titles were not inventoried).
2. **The batches, the hit list, the random ncbi2na resolution and the identities of the shown letters were ported in S08** (DW-12), although they also change outfmt 6 for some inputs: outfmt 6 now matches NCBI in those cases (fixtures `env.batch*`, `hitlist.*`, `ambig.*`). Gate A and the capture of the S02 cases are unchanged.
3. **Explicit rejections** (`... not supported by LOSAT's TBLASTX`):
   - `-window_size 0` (the one-hit word finder, not ported);
   - outfmt 0/7 query or subject deflines that are empty, non-ASCII, contain a byte below a space, or have an empty id (NCBI's reader handles them differently; input reading is S08+);
   - the outfmt 0 title of a subject with hits that NCBI's `NStr::HtmlDecode` changes (`src/objmgr/util/create_defline.cpp:4066`; `report/defline.rs` `ncbi_nucleotide_title_is_decoded`, NCBI's scan and entity table, `src/corelib/ncbistr.cpp:4223-4590`; checked after the search for the shown subjects, `src/algo/blast/format/blast_format.cpp:1540`);
   - a `BATCH_SIZE`, `CHUNK_SIZE` or `OVERLAP_CHUNK_SIZE` that NCBI cannot convert (its `CStringException` names the oracle build's source path), and `BL2SEQ_LEGACY`;
   - `-culling_limit` above 0 (LOSAT's culling keeps other HSPs than `hspfilter_culling.c`);
   - outfmt specifications other than 0, 6 and 7 without fields.
4. **Pending the maintainer: outfmt 0 subject titles that NCBI's `x_CleanAndCompress` reads past** (`src/objmgr/util/create_defline.cpp:219-312`). NCBI tblastx and tblastn crash (SIGSEGV) on them as blastn does ([`punct_defline.py`](punct_defline.py)); exception 2 of `PD-LOSAT-NCBI-DEFECTS` covers BLASTN only. Until the maintainer decides, TBLASTX and TBLASTN reject such subjects with hits in outfmt 0 (`report/defline.rs`, `ncbi_nucleotide_title_reads_past_end`), as BLASTN did before its exception; TBLASTX's rejection follows its outfmt 0 prolog (NCBI had flushed its first lines before it crashed). TBLASTN printed such titles uncleaned before S08.
5. **The ABI v1 stays frozen** (TD-1): its TBLASTX entry keeps only outfmt 6 and its error text (`tblastx_v1_outfmt`) and no hit list size.
6. **Exception classification by the database oracle** (§E) instead of a hand-written diff.
7. **The inventory results and four investigations found differences that LOSAT did not report, most of them older than S08 and visible in outfmt 6 too; all were ported or rejected in S08** (DW-12; [`investigations/`](investigations/)): TBLASTX reads its inputs as BLASTN does (`blastn/input.rs`, with the program's name: deflines, residues, empty records, file names, the order of `-subject`, `-query` and `-out`, `U` as `T`, CFastaReader's title warnings); linking per query and strand (a 3- or 4-nt query shifted the groups of the later queries of its batch, [`shortq.md`](investigations/shortq.md)); the two sorts of tied HSPs above ([`linktie.md`](investigations/linktie.md), [`nisland.md`](investigations/nisland.md)); `BLAST_LargeGapSumE`'s order of evaluation; SEG's left recursion ([`segmask.md`](investigations/segmask.md), shared by BLASTP, TBLASTN and TBLASTX); `-culling_limit` rejected. BLASTX keeps its behavior for the shared SEG and sum statistics until SX (DW-10).
8. **The independent audit's findings on 968fa98c8 were ported or rejected** (`7fbbfad96`; README, "独立監査"; inventory range H): the outfmt 0 title checks of the shown subjects only and NCBI's exact `HtmlDecode` (also BLASTN's check); NCBI's `-evalue` check (`src/algo/blast/core/blast_options.c:1518-1523`); `-seg` tokens split on single spaces with an `int` window and values not above 0 kept at NCBI's defaults (`src/algo/blast/core/blast_filter.c:1147-1154`; BLASTP, TBLASTN, TBLASTX); `CHUNK_SIZE` that is not divisible by 3 (`src/algo/blast/api/local_blast.cpp:98-103`); `BATCH_SIZE` read before the queries (`src/app/blast/tblastx_app.cpp:136-137`); a DEL byte kept in the titles; `BLAST_Cutoffs`' floor of 1 (`src/algo/blast/core/blast_stat.c:4097, 4126-4129`). A standard output closed at the start is not observable (Rust's runtime opens `/dev/null` there before `main`); LOSAT writes the report there and succeeds, a difference put to the maintainer (README). NCBI's texts and exit 1 of the other option checks after parsing (`-threshold 0`, malformed `-seg`) stay with S08+ (TD-13), where LOSAT's parser stops with exit 2.

## G. Inventory of NCBI's path and what the port did

Before the port, seven read-only agents listed every NCBI function and branch on the path that can change the outfmt 0/7 bytes, stderr or the exit status (422 rows, ranges A–G; conventions in `inventory/COMMON.md`; the result agents added 2 rows and the independent audit 14, range H; 438 in all). Each row has the state of LOSAT's TBLASTX before the port (`status`) and the result of the port (`s08_result`, filled from the final code by a second set of read-only agents, with the fixture or oracle run that shows it). The counts are in the gate record ([README.md](README.md), "棚卸し").
