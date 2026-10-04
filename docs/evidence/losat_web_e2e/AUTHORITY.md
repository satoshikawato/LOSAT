# E2e（S08+）：BLASTP・TBLASTN・TBLASTX の既定以外のオプション — NCBI の経路

NCBI のソースは固定 commit `598d8ae6`（`/mnt/c/Users/genom/GitHub/ncbi-blast/c++`。作業では LF に揃えた複製 `~/.cache/losat-web-gui-target/s08p/ncbi/c++` を読んだ）、振る舞いの確認は NCBI BLAST+ 2.17.0（`/home/kawato/micromamba/bin`）で行った。LOSAT の移植箇所には、それぞれの直上に NCBI のファイル・行と断片を書いてある（`s08p/verify_refs.py` で断片が NCBI の行にあることを確かめた）。経路の関数の棚卸しは `INVENTORY.tsv`（13 の範囲、1063 行。範囲ごとの記録は `inventory/`）。この文書は、どの NCBI の経路が LOSAT のどこに当たるかと、S08+ で直したものと拒否するものの対応表である。

範囲の記号：PA・NA・XA（BLASTP・TBLASTN・TBLASTX の app と CLI）、PV・NV（オプションの検査）、PS・NS（得点と統計）、LK（lookup と word finder）、IP・IN（蛋白・核酸の入力）、RP・RN（報告）、XC（TBLASTX の culling）。

## A. app の流れ（3 つの program に共通）

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `blast_args.cpp:3624-3627,2745-2748,2801-2851` | `-outfmt` は引数の組の処理より前に `ParseFormattingString` で読む。前後の C の空白を除き、最初の空白までを `NStr::StringToInt` で読む。読めなければ `BLAST query/options error: '<値>' is not a valid output format`（終了コード 1）、0〜21 の外は `Error: Formatting choice is out of range`（終了コード 255）。`delim=` の字句、tabular・CSV・SAM 以外では custom の指定を捨てる | `blastinput/app.rs` の `parse_formatting_string`（`FormatChoice{number, spec, delimiter}`、`normalized()`）、`report_format`（LOSAT が書かない形式は「not supported by LOSAT's <PROGRAM>」） |
| `blast_args.cpp:2553-2562` | `CBlastDatabaseArgs` が subject を読む。`-subject` が無ければ `Either a BLAST database or subject sequence(s) must be specified`、レコードが無ければ `BLAST engine error: Empty CBlastQueryVector`（終了コード 3） | `app.rs` の `missing_subject_error`・`empty_subjects_error`。BLASTP は `blastn/input.rs` の部品で蛋白を読み、TBLASTN と TBLASTX は `read_nucleotide_subjects`（S08+ で TBLASTX から共有に出した） |
| `blast_args.cpp:3456-3481` | `CStdCmdLineArgs` が query、`-out` の順に開く。`-` は標準入力・標準出力、`-query` の既定値は `-` | 3 つの program の `run`。`ncbi_input_path`・`ncbi_output_path`（256 バイト未満） |
| `blast_args.cpp:3631-3639` | 引数の組の処理（program ごとの順、§D）の後に `Validate`（§E）。失敗は `CInputException` | `BlastpArgs::check_options`、`TblastnArgs::check_options`、TBLASTX の `check_ncbi_options` |
| `blastp_app.cpp:211-214`、`tblastn_app.cpp:213-216`、`tblastx_app.cpp:132-135`、`blast_app_util.cpp:856-860` | 空の query は `Warning: [<program>] Query is Empty!` で終了コード 0。位置の無い stream（パイプ）は空とみなさない | 3 つの `search_cli`。パイプの空の query は明示的な拒否（D10、BLASTN と同じ） |
| `blast_app_util.cpp:888-899` | outfmt 13・14 を標準出力に出すと `Please provide a file name for outfmt 13.` | `app.rs` の `xinclude_check` |
| `blast_app_util.hpp:172-263` | 誤りの種類ごとの文言と終了コード（query/options 1、engine 3、書き込みの失敗 6、`std::exception` 255） | `cli.rs` の `NativeError`、`main.rs` の `exit_on_native_error`（3 つの program で呼ぶ） |
| `ncbiapp.cpp:1031-1044` ほか | app の層が環境と registry を読む | `main.rs` が 3 つの program で `check_ncbi_application_settings` を呼ぶ（S08+ の前は BLASTN だけ） |

## B. 引数の読み方

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `ncbiargs.cpp:118-129` | 整数は 10 進数、読めなければ `0x` の 16 進数 | `value_parsers.rs` の `ncbi_integer` を 3 つの program の全整数に使う（`blastn_count`、`nonnegative_ncbi_integer`、`protein_word_size` ほか）。制約のある引数は `ncbi_integer_at_least` |
| `ncbiargs.cpp:464`、`ncbistr.cpp:1312-1339` | 実数は `NStr::StringToDouble`：最初の文字が数字・小数点・符号でなければ誤り、それ以外は `strtod`（`+inf`・`+nan`・`1e400` を読む） | `ncbi_string_to_double`・`ncbi_double`（`blastp_real`・`tblastn_real`・`tblastx_real`）。16 進数の実数と `1e` の形は明示的な拒否（D6、BLASTN と同じ） |
| `blast_args.cpp:578-583` | `-threshold` は 0 以上の実数（NaN は制約で落ちる） | `ncbi_threshold`（`blastp_threshold`・`tblastn_threshold_value`・`tblastx_threshold`） |
| `ncbiargs.cpp:489-497` | 真偽値は `NStr::StringToBool`（`true/t/yes/y/on/1` と否定、大文字小文字を区別しない） | `ncbi_boolean`（TBLASTN の `-soft_masking`・`-sum_stats`） |
| `blast_args.cpp:387-408` | `-seg` は filtering の handler が読む。`no`・`yes` 以外は 1 つの空白で 3 つに分け、数でなければ `Invalid input for filtering parameters`、3 つでなければ `Invalid number of arguments to filtering option`（終了コード 1） | `app.rs` の `parse_seg_option`（3 つの program が共有。TBLASTN・TBLASTX の既定値は `12 2.2 2.5`） |
| `blast_args.cpp:834-892` | `-comp_based_stats` は最初の文字だけで決める（`0Ff` 0、`1` 1、`2DdTt` 2、`3` 3、それ以外は 0）。blastp だけ 2 文字目の `u`・`U` で unified P。`-ungapped` と 0 以外の組は `Composition-adjusted searched are not supported with an ungapped search…` | `app.rs` の `parse_comp_based_stats` |
| `blast_args.cpp:3158-3163` | `-task` は大文字小文字を区別する集合（blastp・blastp-fast・blastp-short、tblastn・tblastn-fast） | clap の `value_parser`（構文の誤りは承認済みの例外 1） |
| `tabular.cpp:70-99,1393-1408`、`format_flags.cpp:41-195` | custom の field は 1 つの空白で分け、名前に完全に一致するものだけを足す（`std` は既定の 12、`-<名前>` は消す）。知らない名前は黙って捨て、何も残らなければ既定 | BLASTP の `parse_blastp_tabular_fields`・`blastp_tabular_field`。NCBI が書き LOSAT が書かない field（`qgi`・`sallseqid`・`staxid`・`qcovs` など 20）と `delim=` は明示的な拒否。TBLASTN・TBLASTX は custom の指定を拒否する（S08 から） |

## C. task の既定値

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `blast_options_handle.cpp:381-402`、`blast_prot_options.cpp:56-82` | blastp：word size 3、threshold 11、window 40、BLOSUM62 11/1、e-value 10、cbs 2。blastp-fast：word size 5、compressed lookup、threshold 20。blastp-short：word size 2、PAM30 9/1、e-value 20000 | `BlastpTaskOptions::create`。blastp-short は word size 2 のため拒否（§K） |
| `blast_options_handle.cpp:436-447`、`tblastn_options.cpp:54-85` | tblastn：threshold 13、sum statistics、cbs 2、genetic code 1。tblastn-fast：word size 5、compressed lookup、threshold 20 | `TblastnArgs::check_options` |
| `blast_options_local_priv.hpp:612-619` | `CBlastOptionsLocal::SetWordSize`：compressed lookup で word size ≤ 4 なら通常の lookup に、通常の lookup で > 4 なら compressed に戻す（tblastn-fast の `-word_size 3`） | BLASTP の `set_word_size`、TBLASTN の `check_options` |
| `blast_args.cpp:288-301` | 蛋白の query で word size > 4 なら compressed lookup と threshold 19.3（> 5 で 21、> 6 で 20.25） | 同上 |
| `blast_args.cpp:586-623`、`blast_options.c:1174-1209` | `-threshold` を省くと、threshold の整数部が program の既定（blastp 11、tblastn 13）のときだけ行列の推奨値（翻訳した subject は +2）に置き換える | `protein_options.rs` の `suggested_threshold`（`SuggestionProgram`） |
| `blast_args.cpp:485-497`、`blast_options.c:1123-1172` | `-window_size` を省くと行列の推奨値 | `suggested_window_size` |
| `blast_args.cpp:251-286` | `-matrix` を与えると、`-gapopen`・`-gapextend` を省いた側は行列の表の `BLAST_MATRIX_BEST` の行の値 | `protein_gap_existence_extend_params`（表の 2 行目から探す） |

## D. 引数の組の処理の順

- BLASTP（`blastp_args.cpp:44-120`）：Task、BlastDatabase、StdCmdLine、GenericSearch（e-value、gap、word size、threshold の変更）、Filtering（`-seg`）、MatrixName、WordThreshold、WindowSize、HspFiltering、QueryOptions、Formatting（他の program の形式、`Examining 5 or more matches is recommended`）、MT、Remote、CompositionBasedStats、Debug。
- TBLASTN（`tblastn_args.cpp:44-125`）：Task、BlastDatabase、StdCmdLine、GenericSearch、GeneticCode（db）、Gapped、LargestIntron、Filtering、MatrixName、WordThreshold、HspFiltering、WindowSize、QueryOptions、Formatting、MT、Remote、CompositionBasedStats、PsiBlast、Debug。
- TBLASTX（`tblastx_args.cpp:44-118`）：BlastDatabase、StdCmdLine、GenericSearch、LargestIntron、Filtering、MatrixName、WordThreshold、HspFiltering、WindowSize、QueryOptions、GeneticCode（query、db）、Formatting、MT、Remote、Debug。

LOSAT は誤りを出す handler だけを同じ順に置いた（filtering、formatting、composition、Validate）。

## E. `BLAST_ValidateOptions`（`blast_options.c:1749-1810`）

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `blast_options.c:913-943` | gapped の検索で行列と gap の組を表から探す（blastp・tblastn は IDENTITY も可）。行列が無ければ `BLAST_PrintMatrixMessage`（9 の行列の一覧、`snprintf` の 1024 バイト）、組が無ければ `BLAST_PrintAllowedValues`（2048 バイト） | `protein_options.rs` の `karlin_blk_gapped_load_from_tables`・`print_matrix_message`・`print_allowed_values`、表は `protein_tables.rs`（`gen_protein_tables.py` で `blast_stat.c` から生成、`--check` で一致を確かめる） |
| `blast_options.c:1303-1395` | threshold ≤ 0 は `Non-zero threshold required`。blastp・tblastn・blastx は word size > 7 で `Word-size must be less than 8…`、> 5 で compressed でなければ、compressed で 5〜7 でなければ誤り。tblastx は > 4 で `Word-size must be less than 6 for protein comparison` | `validate_protein_options`（BLASTP と TBLASTN が共有）、TBLASTX の `check_ncbi_options` |
| `blast_options.c:1518-1523` | e-value ≤ 0（cutoff の option は無い）は `expect value or cutoff score must be greater than zero` | 同上 |
| `blast_options.c:1783-1795` | IDENTITY で word size > 5 は誤り | 同上 |

## F. lookup の threshold と `(Int4)` の変換

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `blast_aalookup.c:245` | 通常の lookup の threshold は `(Int4)opt->threshold`。x86-64 の `cvttsd2si` は範囲の外と inf で `INT_MIN`（全ての word が隣接語になる） | BLASTP・TBLASTX の lookup で `core/blast_util.rs` の `ncbi_int4_from_double`（D4） |
| `blast_aalookup.c:1292` | compressed lookup は `(Int4)(kMatrixScale * opt->threshold)`（100 倍） | BLASTP の `try_build_blosum62_compressed_lookup`（BLASTX は従来どおり） |
| `blast_aalookup.c:759-768` | overflow の cell の bank は 1024 個。release build に ASSERT は無く、使い切ると `overflow_banks` の外に書いて heap を壊す | 使い切る表は `overflow_banks_exhausted` にし、BLASTP は明示的な拒否（NCBI の欠陥、D3）。BLASTX は従来の panic のまま |
| `blast_parameters.c:457-466`、`blast_kappa.c:2372-2384` | 予備と最後の gapped の X-drop は `(Int4)` の `MAX(double, double)` | TBLASTN の `local_extension_final_xdrop`・`local_kappa_redo_params`・予備の X-drop（`-xdrop_gap 1e9` の panic と `-xdrop_gap_final 1e9` の差を直した。`-comp_based_stats 0` の差は下の fence の行） |
| `blast_hits.c:1224-1228`、`blast_traceback.c:425-433,1651-1686` | 部分翻訳（HSP の前後 99 塩基）の両端に fence（先頭からの翻訳も左端に）。traceback が fence に触れると subject を全翻訳でやり直し、`stat_length` は最後に翻訳した長さ | TBLASTN の `TargetTranslation::get`（`base` は符号付き、先頭からの窓の fence は位置 −1）。第 1 回の監査 TN-1（`-comp_based_stats 0` の小さい X-drop で e-value が違った） |
| `blast_extend.c:52-57` | 対角の表の長さは query の長さ＋window を超える 2 の冪。和が 2^30 を超えると `Int4` で終わらないか回り込む | BLASTP・TBLASTX は明示的な拒否（D8・D11、`app.rs` の `check_diag_table_window`） |

## G. hit list の大きさ

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `blast_hits.c:43-70` | 予備の hit list の大きさは `int` で計算して折り返す（cbs ありで 500 以下は 1050、それ以外は 2h+50。cbs 無しは MIN(MAX(2h,10), h+50)） | TBLASTN は BLASTN の `get_prelim_hitlist_size` を使う。正でない大きさ（NCBI は落ちる）は明示的な拒否 |
| `blast_args.cpp:2894-2978` | `-max_target_seqs` は hit list と outfmt 0 の説明・alignment の数を決める | 3 つの program |

## H. 報告

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `blast_format.cpp:2266-2285` | epilog の `Matrix:` は入力どおりの名前、threshold は double を C++ の stream の既定（`%g`、`inf`） | BLASTP・TBLASTN の `matrix_name`、`BlastpPairwiseReport.word_threshold: f64`・TBLASTX の `TblastxPairwiseReport.word_threshold: f64`、`cpp_default_double`（inf・nan を足した） |
| `blast_seqalign.cpp:672-674,1484-1490` | 得点 0 の HSP は Seq-align にならず、どの報告にも出ない。HSP が全て得点 0 の subject は報告に出ない | BLASTP（最終の hit list から除く）、TBLASTN（`drop_zero_score_hsps`） |
| `blast_format.cpp:1547-1581`、`blast_aux.cpp:936-951`、`showalign.cpp:2518-2521,4305-4318` | outfmt 0 は SEG で mask した query の列を小文字にする（query 全体を覆う mask は出さない。中線は大文字で比べる） | BLASTP の `masked_query_rows`、`write_alignment_with_sequences` の中線 |
| `blast_kappa.c:331-342`、`showalign.cpp:3599-3604` | `Method:` は HSP ごとの行列の調整の規則 | BLASTP（S08+ の最初のコミット `68469680b`） |
| `blast_format.cpp`（`PrintProlog`） | prolog は最初の query の batch を読む前に書く | BLASTP・TBLASTN の `write_*_pairwise_prolog`（`BATCH_SIZE=0` の経路） |

## I. 環境

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `blast_app_util.cpp:204-210,732-737` | `BL2SEQ_LEGACY`（subject ごとの検索と bl2seq の報告）、`PRE_FETCH_SEQS_LIMIT`（`int` でなければ例外） | `app.rs` の `check_unsupported_environment`（TBLASTX から共有に出し、BLASTP・TBLASTN でも呼ぶ。明示的な拒否） |
| `blast_stat.c:922-923` | `OLD_FSC` があると Gumbel の block（有限の長さの補正）を作らない | `check_old_fsc`（BLASTP・TBLASTN。明示的な拒否） |
| `blast_input_aux.cpp:85-135` | `BATCH_SIZE`（blastp 10000、tblastn 20000、tblastx 10002）。0 は最初の batch が空で、prolog の後に `Empty CBlastQueryVector`（終了コード 3）。`int` でなければ例外 | `query_batch_size`。0 は NCBI と同じ（BLASTP・TBLASTN・TBLASTX）、`int` でない値は明示的な拒否 |
| `split_query_cxx.cpp:55-61`、`local_blast.cpp:98-103` | `CHUNK_SIZE`・`OVERLAP_CHUNK_SIZE`（翻訳する query は 3 の倍数） | `check_query_split_environment(program, translated_query)` |
| `blast_hits.c:43-70`、`blast_hspstream.c:816-864` | `ADAPTIVE_CBS` | BLASTP は予備の大きさを移植済み、TBLASTN は明示的な拒否 |

## J. 入力

| NCBI | 振る舞い | LOSAT |
|---|---|---|
| `fasta.cpp:966-979` | 蛋白の残基でない文字は警告して落とす | `blastn/input.rs` の `check_protein_input_of`・`check_protein_residues_of`（明示的な拒否。BLASTP の query・subject、TBLASTN の query） |
| `fasta.cpp:375-384` ほか | 核酸の subject：最初の定義行の前の空行と注釈、IUPAC でない残基、中身の無いレコード、`U`、題の警告 | `read_nucleotide_subjects`（TBLASTN は S08+ で TBLASTX と同じ読み方になった。最初の作業 2） |

## K. LOSAT が対応する値と拒否する値（S12 の検索画面の入力）

明示的な拒否の文言は全て「… is not supported by LOSAT's <PROGRAM>」。

**BLASTP**：対応する値は、`-task blastp`・`blastp-fast`、BLOSUM62 と gap 11/1、word size 3 と compressed の 5（`-task blastp-fast` か `-word_size 5`）、threshold は任意の正の実数（`+inf` を含む。compressed で overflow の bank を使い切る値は拒否）、`-window_size` は 0 以上の整数（0 は one-hit。query の batch の長さとの和が 2^30 を超える値は拒否、D11）、`-comp_based_stats` は mode 2（`2`・`D`・`d`・`T`・`t` と、それで始まる任意の文字列）、`-seg`（`no`・`yes`・3 つの値）、`-evalue`（正の有限の実数。`+inf`・`1e999` は拒否、D12）、`-max_target_seqs`（予備の hit list の大きさが NCBI の `int` で正になる値。2147483624 以上は NCBI と同じく 2〜48 に回り込む）・`-max_hsps`（1 以上）、`-outfmt` の 0・6・7（NCBI の書き方。6・7 は LOSAT が書く field の custom の指定）。拒否する値は、BLOSUM62 以外の行列と 11/1 以外の gap（NCBI の表にある組）、word size 2・4・6・7（`-task blastp-short` を含む）、`-comp_based_stats` の 0・1・3 と unified P（`u`）、`-ungapped`、`-use_sw_tback`、§B の field と `delim=`、NCBI の blastp の option のうち LOSAT に無いもの（`cli.rs` の `is_unported_blastp_arg`：`-db` の一族、`-culling_limit`、`-best_hit_*`、`-subject_besthit`、`-lcase_masking`、`-soft_masking`、`-xdrop_*`、`-searchsp`、`-dbsize`、`-qcov_hsp_perc`、`-query_loc`、`-subject_loc`、`-num_descriptions`・`-num_alignments` ほか、`-h`、`-version`）、NCBI C++ Toolkit の option（`-help-full`・`-xmlhelp`・`-logfile`・`-conffile`・`-version-full*`。4 つの program）。

**TBLASTN**：対応する値は、`-task tblastn`、BLOSUM62 の 11/1・word size 3・threshold 13・window 40（cbs 0 と 2）と、BLOSUM45 の 14/2・word size 2・threshold 16・window 60（cbs 0）、`-db_gencode`（承認済みの例外。32 を含む）、`-seg`、`-soft_masking`・`-sum_stats`（NCBI の真偽値）、`-lcase_masking`、`-xdrop_gap`・`-xdrop_gap_final`（任意の実数。1e7 bit 以上の巨大な値は LOSAT が遅く大量の記憶を使う、S08+b の TN-4）、`-evalue`（正の有限の実数。無限大は拒否、D12）、`-max_target_seqs`（予備の大きさが正になる値）、`-outfmt` の 0・6・7（custom の指定は拒否）。拒否する値は、上の 2 つ以外の行列・gap・word size・threshold・window の組（`-task tblastn-fast`、compressed の lookup、`-window_size 0`、`-threshold 13.5` を含む）、`-comp_based_stats` の 1・3、`-ungapped`、`-max_intron_length` の 0 以外、NCBI の tblastn の option のうち LOSAT に無いもの（`is_unported_tblastn_arg`：`-db`・`-in_pssm`・`-remote`・`-subject_loc`・`-use_sw_tback`・`-xdrop_ungap`・`-max_hsps`・`-culling_limit` ほか）。

**TBLASTX**：対応する値は、word size 3、threshold は任意の正の実数（`+inf` を含む）、`-window_size` の 1 以上（query の batch の翻訳の長さとの和が 2^30 を超える値は拒否、D8・D11）、`-seg`、`-query_gencode`・`-db_gencode`（承認済みの例外）、`-culling_limit`（0 以上。S08+ で移植）、`-evalue`、`-max_target_seqs`、`-outfmt` の 0・6・7（custom の指定は拒否）。拒否する値は、word size 2・4、`-window_size 0`、NCBI の tblastx の option のうち LOSAT に無いもの（`is_unported_tblastx_arg`：`-matrix`・`-max_hsps`・`-sum_stats`・`-lcase_masking`・`-soft_masking`・`-strand`・`-max_intron_length`・`-xdrop_ungap`・`-best_hit_*`・`-subject_besthit`・`-num_descriptions`・`-num_alignments` ほか）。

## L. TBLASTX の culling（範囲 XC）

`hspfilter_culling.c` を移植した（`algorithm/tblastx/hsp_culling.rs`、`08396413b`）：予備の段の writer は merit を N+3（N > 1。`Int4` で折り返す）にして全 subject の HSP を context ごとの 1 つの木に入れ、hit list の切り詰めは traceback の段で行い、pipe の段は merit N で e-value の順（`report.rs` の `culled_hit_order`）。NCBI の fixture 98 件のうち 81 件が一致、残りは LOSAT に無い option（`-max_hsps`・`-best_hit_*`・`-subject_besthit`）と引数の誤り（例外 1）。

## M. 承認済みの例外と判断

- D1：TBLASTX の `-culling_limit` を移植（NCBI は決定的）。D3：NCBI が落ちる組（IDENTITY の word size 5・6、`-use_sw_tback -ungapped -cbs 0`、compressed の lookup の bank の使い切り、予備の hit list の大きさが正でない値）は明示的な拒否。D4：`double` から `Int4` への未定義の変換は x86-64 の `cvttsd2si`（`INT_MIN`）を再現。D5：unified P は NCBI の結果が実行ごとに変わるので拒否。D6：実数は NCBI の `CArg_Double` の書き方、16 進数の実数は拒否。D7：未知の行列と `-ungapped` は NCBI の engine error（現在は `-ungapped` 自体を拒否）。D8：query の長さ＋window が (2^30, 2^31−1] で NCBI が止まらない組は拒否（第 1 回の監査の後に実装、`check_diag_table_window`）。D9：`-remote`・`-db` の一族・検索の戦略・`-html`・`-parse_deflines`・隠れた toolkit の option・`-version` は拒否、`-h`・`-help` は LOSAT の help（例外 1）。D10：パイプの空の query は拒否。D11（第 1 回の監査 TX-1）：query の長さ＋window が 2^31 − 1 を超えると NCBI の `Int4` が回り込み、対角の表が 1 つの cell になり、後の offset も回り込んで hit が無い（決まった結果）。実用が無く再現の費用が高いので、DW-16 の扱いで明示的な拒否（D8 と合わせて和が 2^30 を超える値。BLASTP・TBLASTX）。D12（第 1 回の監査 BP-4・TN-6）：無限大の `-evalue`（`+inf`、`1e999`）で NCBI の blastp と tblastn は入力によって SIGSEGV（300 の subject、compressed の lookup）。LOSAT はいつ落ちるかを事前に知れないので、明示的な拒否（D3 の扱い。有限の `1e308` は実行し NCBI と一致）。推奨の別案は、近い有限の値の NCBI の出力と一致を確かめた LOSAT の結果を承認済みの例外にすること（DW-15）。D11・D12 は保守者の確認を待つ（まとめて諮る）。TBLASTX は無限大の `-evalue` で NCBI も一致するので受け付ける。
- TBLASTN の `-db_gencode` の承認済みの例外（`AGENTS.md`）：sweep の oracle（NCBI の `-db` の検索）は、NCBI の `-db` と `-subject` の結果が既定の遺伝暗号でも e-value・cbs によって違う（予備の cutoff と連結の違い）ため、厳密な oracle でない。TLOSAN Stage E/G の比較専用の C++ API の oracle（subject の検索の経路で、選んだ遺伝暗号を sequence source に与える。`docs/evidence/tlosan_stage_e/tblastn_stage_e_local_oracle.cpp`、[`gates/build_api_oracle.sh`](gates/build_api_oracle.sh) で Stage G の実行ファイルとバイト一致に作り直した）で確かめた：遺伝暗号 1 で CLI と一致（較正）、sweep の 27 の綴りの outfmt 0/6/7 の 81 件で LOSAT と一致（[`gencode_api/check.tsv`](gencode_api/check.tsv)、[`gencode_api_check.py`](gencode_api_check.py)）。

## N. 残り（次のセッション S08+b）

第 1 回の独立監査（[`audit/ROUND1.md`](audit/ROUND1.md)）で見つかり、このセッションで直さなかったもの。どれも NCBI とバイト一致にするか明示的に拒否する。

- TN-2：TBLASTN の `-comp_based_stats 0`、hard mask の `-lcase_masking -seg no`、少ない `-max_target_seqs` で、和の統計の e-value が 1〜4 % 違う（`-sum_stats false` では一致）。NCBI の linking の `Blast_HSPListGetEvalues` の呼び出し（subject の長さ/3）と Spouge の `db_length` を、LOSAT の linking の値と比べる（TLOSAN Stage D の trace の shim `ncbi_d_call_trace.c` で NCBI の値が取れる）。
- TN-4：巨大な `-xdrop_gap`・`-xdrop_gap_final`（1e7〜5e8 bit）で LOSAT の DP の帯の `resize` が全 cell を書き、数十秒・数 GB（NCBI は `malloc` で 1 秒未満）。BLASTP・TBLASTN の共有の経路（`blastp/gapalign.rs` の `gap_dp_reserve_initial`・`gap_dp_reserve_band` と traceback の状態）を、書く前に読まない cell を初期化しない形にするか、閾値で明示的に拒否する。V-PERF が要る。
- TN-5：BLOSUM45 の組で `-evalue` 5000 以上の同じ得点の HSP の frame・座標。
- RP-4（既定の option、前から）：60000 残基の query の長い一致の composition の行列の調整で bit score が違う（入力は `audit/round1/rp4/`）。
- `-out -version`（BP-8 の残り）：NCBI は `-version` を version の option として読む。
- 第 2 回の独立監査（上を直した後、4 観点）。
- 保守者の確認：D11、D12（推奨の別案の承認済みの例外を含む）。
