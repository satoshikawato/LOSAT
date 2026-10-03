# BLASTN `-task dc-megablast` and `-task blastn-short` authority (Session SD, stage E2i)

This record fixes what LOSAT's BLASTN reproduces for NCBI's discontiguous megablast and blastn-short tasks, and where. Section A gives the option values of the tasks, B the arguments and option checks, C the lookup table, D the subject scan, E the ungapped extension and the diagonal structures, F the gapped stage, G the report, H the decisions of the session, I the inventory. Every quotation is on the cited lines of the pinned source (checked by [`check_authority.py`](check_authority.py)).

| Item | Pinned value |
| --- | --- |
| NCBI source (sole authority) | `/mnt/c/Users/genom/GitHub/ncbi-blast/`, commit `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`. Paths below are relative to `c++/`. Line numbers count lines after removing CR. |
| NCBI BLAST+ (comparison oracle only) | `blastn` 2.17.0+ (package `blast 2.17.0`, build `Aug 11 2025 09:46:06`) in `/home/kawato/micromamba/bin/`. No `.ncbirc`; `BATCH_SIZE`, `CHUNK_SIZE`, `OVERLAP_CHUNK_SIZE`, `BL2SEQ_LEGACY`, `CTOOLKIT_COMPATIBLE` unset unless a case sets them. |
| LOSAT | Branch `feature/losat-web-gui`; the port is `90c5f0181`, the fixtures `af03129ee`. LOSAT line numbers below are those of `90c5f0181` (paths relative to `LOSAT/src/`). |
| Fixtures | `dc.*` and `short.*` rows of [`LOSAT/tests/outfmt0_manifest.tsv`](../../../LOSAT/tests/outfmt0_manifest.tsv) (21, outfmt 0 and 7) and of [`LOSAT/tests/blastn_regression_fixtures.py`](../../../LOSAT/tests/blastn_regression_fixtures.py) (54); inputs made by [`make_inputs.py`](make_inputs.py). |
| Inventory | [`INVENTORY.tsv`](INVENTORY.tsv) (built by [`build_inventory.py`](build_inventory.py) from `inventory/`). |

---

## A. The tasks' option values

blastn-short is NCBI's blastn handle with five overrides; dc-megablast is the megablast handle with the discontiguous overrides. `CDiscNucleotideOptionsHandle`'s base constructor calls the base class's defaults, then its own constructor calls `SetDefaults` again, which reaches the overrides; the hit-saving defaults are not overridden, so dc-megablast keeps the megablast min diag separation (6).

| Value | NCBI file:line | Snippet | LOSAT |
|---|---|---|---|
| blastn-short | `src/algo/blast/api/blast_options_handle.cpp:354-361` | `opts->SetMatchReward(1);` / `opts->SetMismatchPenalty(-3);` / `opts->SetEvalueThreshold(1000);` / `opts->SetWordSize(7);` / `opts->ClearFilterOptions();` | `algorithm/blastn/coordination.rs:333` (`task_defaults`) |
| dc-megablast handle | `src/algo/blast/api/blast_options_handle.cpp:377-380` | `retval = CBlastOptionsFactory::Create(eDiscMegablast, locality);` | `task_defaults` |
| dc template and word size | `src/algo/blast/api/disc_nucl_options.cpp:55-64` | `SetTemplateType(0);` / `SetTemplateLength(18);` / `SetWordSize(BLAST_WORDSIZE_NUCL);` | `task_defaults`, `determine_template` (`coordination.rs:483`) |
| dc window and ungapped X-drop | `src/algo/blast/api/disc_nucl_options.cpp:66-74` | `SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_NUCL);` / `SetWindowSize(BLAST_WINDOW_SIZE_DISC);` | `TaskConfig::window_size` (40); the ungapped X-drop of every task but megablast (`blast_engine/run.rs` `x_dropoff_init`) |
| dc gapped extension | `src/algo/blast/api/disc_nucl_options.cpp:76-84` | `SetGapXDropoff(BLAST_GAP_X_DROPOFF_NUCL);` / `SetGapExtnAlgorithm(eDynProgScoreOnly);` / `SetGapTracebackAlgorithm(eDynProgTbck);` | `configure_task` (`use_dp`, `x_drop_gapped` 30, `x_drop_final` 100) |
| dc scores | `src/algo/blast/api/disc_nucl_options.cpp:86-90` | `CBlastNucleotideOptionsHandle::SetScoringOptionsDefaults();` | reward 2, penalty -3, gaps 5/2 |
| dc min diag separation | `src/algo/blast/api/blast_nucl_options.cpp:250-259` | `SetMinDiagSeparation(6);` | `MIN_DIAG_SEPARATION_MEGABLAST` for dc-megablast; 50 for blastn and blastn-short (`blast_nucl_options.cpp:239`) |
| e-value of blastn-short | `src/algo/blast/blastinput/blast_args.cpp:142-146,253-255` | `des += " (1000 for blastn-short)";` / `arg_desc.AddOptionalKey(kArgEvalue, "evalue", des, CArgDescriptions::eDouble);` / `opt.SetEvalueThreshold(args[kArgEvalue].AsDouble());` | `algorithm/blastn/args.rs:82` (`evalue: Option<f64>`), `determine_evalue` (`coordination.rs:412`) |
| no DUST for blastn-short | `src/algo/blast/api/blast_options_cxx.cpp:1022-1031` | `SetDustFiltering(false);` / `SetMaskAtHash(false);` | `args.rs:302` (`resolve_dust`), `task_dust_by_default` |
| mask at hash restored | `src/algo/blast/blastinput/blast_args.cpp:387-391` | `if (args[kArgLookupTableMaskingOnly]) {` / `opt.SetMaskAtHash(args[kArgLookupTableMaskingOnly].AsBoolean());` | unchanged (the default `-soft_masking true` restores what `ClearFilterOptions` cleared; LOSAT always masks at hash and rejects `-soft_masking`) |
| query chunk size | `src/algo/blast/api/local_blast.cpp:64-72` | `case eBlastn:` / `retval = 1000000;` / `case eDiscMegablast:` / `retval = 5000000;` | `task_uses_megablast_chunks` (`coordination.rs:452`), `blast_engine/run.rs:5523` |

## B. Arguments and option checks

| Step | NCBI file:line | Snippet | LOSAT |
|---|---|---|---|
| `-template_type` | `src/algo/blast/blastinput/blast_args.cpp:708-717` | `arg_desc.AddOptionalKey(kArgDMBTemplateType, "type",` / `kTemplType_CodingAndOptimal));` / `CArgDescriptions::eRequires,` | `args.rs:55`, `blastinput/value_parsers.rs:357` (`blastn_template_type`); clap `requires` |
| `-template_length` | `src/algo/blast/blastinput/blast_args.cpp:719-730` | `allowed_values.insert(21);` / `new CArgAllowIntegerSet(allowed_values));` | `args.rs:57`, `value_parsers.rs:382` (`blastn_template_length`) |
| applied to any task | `src/algo/blast/blastinput/blast_args.cpp:743-763` | `options.SetMBTemplateType(static_cast<unsigned char>(temp_type));` / `options.SetMBTemplateLength(tlen);` | `determine_template` (`coordination.rs:483`) |
| check order | `src/algo/blast/core/blast_options.c:1759-1776` | `if ((status = LookupTableOptionsValidate(program_number,` / `if ((status = BlastHitSavingOptionsValidate(program_number, hit_options,` | `algorithm/blastn/scoring.rs` `check_scoring_options`: penalty, gap extension, word size ≤ 100, template, e-value, greedy |
| template word size | `src/algo/blast/core/blast_options.c:1252-1261` | `if (template_length == 0)` / `if (word_size != 11 && word_size != 12) {` / `"size must be either 11 or 12");` | `scoring.rs:144` |
| template lookup table | `src/algo/blast/core/blast_options.c:1399-1413` | `options->mb_template_length > 0) {` / `} else if (options->lut_type != eMBLookupTable) {` / `"Invalid lookup table type for discontiguous Mega BLAST");` | `scoring.rs` (megablast and dc-megablast have the MB table type; blastn and blastn-short the standard one) |
| task list | `src/algo/blast/api/blast_options_handle.cpp:344-351` | `!NStr::CompareNocase(task, "blastn-short") ||` / `!NStr::CompareNocase(task, "rmblastn") ||` | `value_parsers.rs:320` (`blastn_task`: rmblastn rejected explicitly) |

Parser errors (an option without its partner, values outside the constraints) are LOSAT's parser text and exit 2 under the approved `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`; the checks after parsing are NCBI's texts (`BLAST query/options error: ...`, exit 1).

## C. The discontiguous lookup table

| Step | NCBI file:line | Snippet | LOSAT |
|---|---|---|---|
| table choice | `src/algo/blast/core/blast_nalookup.c:62-67` | `if (lookup_options->mb_template_length > 0) {` / `*lut_width = lookup_options->word_size;` / `return eMBLookupTable;` | `coordination.rs:682` (`choose_na_lookup_table`), `blast_engine/run.rs:6595` |
| table size | `src/algo/blast/core/blast_nalookup.c:1256-1260` | `mb_lt->lut_word_length = lut_width;` / `mb_lt->hashsize = 1ULL << (BITS_PER_NUC * mb_lt->lut_word_length);` | `lookup.rs:2004` (`build_disc_mb_lookup`) |
| discontiguous branch | `src/algo/blast/core/blast_nalookup.c:1327-1331` | `mb_lt->scan_step = 1;` / `status = s_FillDiscMBTable(query, location, mb_lt, lookup_options);` | `coordination.rs:1162` |
| template of the options | `src/algo/blast/core/blast_nalookup.c:603-643` | `s_GetDiscTemplateType(Int4 weight, Uint1 length,` / `return eDiscTemplateContiguous; /* All unsupported cases default to 0 */` | `algorithm/blastn/disc_lookup.rs:168` (`get_disc_template_type`) |
| second template | `src/algo/blast/core/blast_nalookup.c:711-716` | `if (kTwoTemplates) {` / `int temp_int = template_type + 1;` | `disc_lookup.rs:116` (`DiscTemplateType::next`) |
| words, 1-based starts | `src/algo/blast/core/blast_nalookup.c:749-790` | `from = loc->ssr->left - (template_length - 2);` / `accum = 0;` / `ecode1 = ComputeDiscontiguousIndex(accum, template_type);` / `mb_lt->next_pos[index] = mb_lt->hashtable[ecode1];` | `build_disc_mb_lookup` (locations: the unmasked ranges of each context, newest first in the chains) |
| second table, shared PV | `src/algo/blast/core/blast_nalookup.c:793-806` | `ecode2 = ComputeDiscontiguousIndex(accum, second_template_type);` / `PV_SET(pv_array, ecode2, pv_array_bts);` / `mb_lt->hashtable2[ecode2] = index;` | `build_disc_mb_lookup` |
| longest chain | `src/algo/blast/core/blast_nalookup.c:809-822` | `mb_lt->longest_chain = longest_chain + 1;` / `mb_lt->longest_chain += longest_chain + 1;` | `build_disc_mb_lookup` (sizes the offset array only) |
| word index | `include/algo/blast/core/blast_nalookup.h:540-588` | `static NCBI_INLINE Int4 ComputeDiscontiguousIndex(Uint8 accum,` / `index = DiscontigIndex_11_18_Coding(accum);` | `disc_lookup.rs:558` and the 12 `discontig_index_*`; unit test `templates_select_the_documented_positions` reads every mask as the template strings of `blast_nalookup.h:195-210` |

## D. The subject scan

| Step | NCBI file:line | Snippet | LOSAT |
|---|---|---|---|
| word lengths | `src/algo/blast/core/na_ungapped.c:1624-1634` | `if (lookup->discontiguous) {` / `word_length = lookup->template_length;` / `lut_word_length = lookup->template_length;` | `lookup.rs:940` (`TwoStageLookup::lut_word_length`, `word_length`) |
| scan range, masked subject | `src/algo/blast/core/na_ungapped.c:1651-1667` | `((BlastMBLookupTable *) lookup_wrap->lut)->discontiguous) {` / `scan_range[2] = subject->seq_ranges[0].right - lut_word_length;` | `blast_engine/run.rs:2829` (`scan_subject_disc_words_with_ranges`) |
| scan choice | `src/algo/blast/core/blast_nascan.c:2617-2627` | `if (mb_lt->two_templates)` / `mb_lt->scansub_callback = (void *)s_MB_DiscWordScanSubject_11_18_1;` | `disc_lookup.rs:600` (`choose_disc_scan_subject`) |
| the scans | `src/algo/blast/core/blast_nascan.c:2202-2273,2289-2367` | `index = ComputeDiscontiguousIndex(accum >> 6, template_type);` / `index2 = ComputeDiscontiguousIndex(accum, second_template_type);` | `disc_lookup.rs:727,874,1044,1238`; unit test `scans_give_the_index_of_every_word_from_any_start` |
| two tables per word | `src/algo/blast/core/blast_nascan.c:1457-1469` | `#define MB_ACCESS_HITS2()` / `total_hits += s_BlastMBLookupRetrieve2(mb_lt,` | `blast_engine/run.rs:9070` (`on_word`: the first template's chain, then the second's); `lookup.rs:893` (`for_each_hit2`) |

NCBI looks the words up inside the scan and stops when its offset array is full, to resume at the same offset; LOSAT streams the words in the same order to the same per-hit extension, which gives the same hits (the extension of a hit depends only on the hits before it).

## E. Ungapped extension and the diagonal structures

| Step | NCBI file:line | Snippet | LOSAT |
|---|---|---|---|
| extension of a discontiguous table | `src/algo/blast/core/na_ungapped.c:1798-1799` | `if (lut->lut_word_length == lut->word_length || lut->discontiguous)` / `lut->extend_callback = (void *)s_BlastNaExtendDirect;` | the two-stage branch with `word_length == lut_word_length` (no word extension) |
| word length of the hit | `src/algo/blast/core/na_ungapped.c:979-983` | `word_length = (lut->discontiguous) ? lut->template_length : lut->word_length;` | as above |
| no mask check | `src/algo/blast/core/na_ungapped.c:531` | `if (word_length == lut_word_length) return 1;` | `type_of_word` |
| two hits | `src/algo/blast/core/na_ungapped.c:674-718` | `if (two_hits && (hit_saved || s_end_pos > last_hit + window_size )) {` / `hit_ready = 0;` | `blast_engine/run.rs:7881` (`two_hits = window_size > 0`); Delta = MIN(0, 40 − template length) = 0 |
| extension routine | `src/algo/blast/core/na_ungapped.c:740-750` | `(word_params->matrix_only_scoring || word_length < 11))` / `s_NuclUngappedExtend(query, subject, matrix, q_off, s_end, s_off,` | the approximate extension (template length ≥ 16) with `s_end` = start + template length |
| diagonal state | `src/algo/blast/core/na_ungapped.c:767-772` | `hit_level_array[real_diag].flag = hit_ready;` / `diag_table->hit_len_array[real_diag] = (hit_ready) ? 0 : s_end_pos - s_off_pos;` | as above |
| hash insert window | `src/algo/blast/core/na_ungapped.c:939-941` | `hit_ready, s_off_pos, window_size + Delta + 1);` | `blast_engine/run.rs:893` (`diag_hash_insert_window`) |
| table size and offset | `src/algo/blast/core/blast_extend.c:54-64` | `while (diag_array_length < (qlen+window_size))` / `diag_table->offset = window_size;` | `blast_engine/run.rs:6692` (`window_size = config.window_size`), `SubjectScratch::new` (`run.rs:3766`) |
| hit lengths | `src/algo/blast/core/blast_extend.c:147-150` | `if (word_params->options->window_size) {` / `diag_table->hit_len_array = (Uint1 *)` | the diagonal setup of the search loop |
| per subject | `src/algo/blast/core/blast_extend.c:168-173` | `ewp->diag_table->offset += subject_length + ewp->diag_table->window;` | `blast_engine/run.rs:652` (`advance_diag_table_offset`) |

## F. The gapped stage

dc-megablast and blastn-short extend with dynamic programming (score-only, then traceback), as `-task blastn`; the only new combination is dc-megablast's min diag separation 6 with DP. Every reader of the separation takes the value of the task (`blast_engine/run.rs` passes `config.min_diag_separation` to the interval tree's containment test at the preliminary, traceback and final purges); none assumes the greedy extension.

## G. The report

| Element | NCBI file:line | Snippet | LOSAT |
|---|---|---|---|
| megablast reference only for the task megablast | `src/app/blast/blastn_app.cpp:244` | `m_CmdLineArgs->GetTask() == "megablast",` | `blast_engine/run.rs` (`megablast = args.task == "megablast"`): dc-megablast and blastn-short print the Gapped BLAST reference |
| gap extension of 0 | `src/algo/blast/format/blast_format.cpp:2272` | `if ((m_Program == "megablast" || m_Program == "blastn") && options.GetGapExtensionCost() == 0)` | `report/pairwise.rs` `write_blastn_final_footer` (`zero_gap_extension_formula`) |
| window line | `src/algo/blast/format/blast_format.cpp:2285-2288` | `if (options.GetWindowSize()) {` / `m_Outfile << "Window for multiple hits: " <<` | `report/pairwise.rs:2096` |

The program line (`BLASTN 2.17.0+`), the outfmt 7 header and the Karlin-Altschul blocks do not depend on the task beyond the scores (`Matrix: blastn matrix 1 -3` for blastn-short). The template options change nothing in the report.

## H. Decisions of the session

1. A template with `-task megablast` (NCBI applies `-template_type`/`-template_length` to any task): ported, as NCBI runs it (the discontiguous table with megablast's greedy extension and one hit; fixtures `dc.megablast_template`, `dc.megablast_w11_coding_18`). With blastn and blastn-short NCBI's checks reject it, and LOSAT gives the same messages.
2. Kept as explicit rejections: `-task rmblastn` (matrix scoring and masklevel), and the options of `cli.rs` `is_unported_blastn_arg` that the tasks' defaults set internally (`-window_size`, `-off_diagonal_range`, `-xdrop_ungap`, `-xdrop_gap`, `-xdrop_gap_final`, `-no_greedy`, `-ungapped`, `-soft_masking`, `-use_index`, `-index_name`, `-min_raw_gapped_score`).
3. ABI v1 (`web_api.rs`) stays at megablast and blastn (plan TD-1 freezes v1); ABI v2 describes and runs both tasks.
4. No NCBI defect was met (`PD-LOSAT-NCBI-DEFECTS` not used).

## I. Inventory

[`INVENTORY.tsv`](INVENTORY.tsv): the rows of six ranges (A application and arguments, B the C++ options layer and the option checks, C core setup and parameters, D lookup/scan/ungapped extension, E gapped extension/traceback/hit saving, F formatting) written before the port by read-only agents (`inventory/<R>.tsv`, conventions `inventory/COMMON.md`; `status` is LOSAT before the port, `a92fa902f`), and the result of each row after the port (`inventory/result_<R>.tsv`, conventions `inventory/RESULT_COMMON.md`, on `90c5f0181`).
