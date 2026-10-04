# Session SX — BLASTX の統合（条件付き）

## INSTRUCTION PROMPT

LOSAT の BLASTX（LOSATX）を LOSAT Web の構成に組み入れる。先に [セッション README](README.md) の共通規則を読み、特に規則 4 に従う。完了条件の正本は、総合計画書 §7 の SX の行である。背景は計画の DW-10（BLASTX は LOSATX 計画の認証の後に扱う）と DW-11（BLASTX も領域の指定の対象にする）である。

入口の条件：LOSATX 計画（`docs/losatx_blastx_v0.2.0_plan.md`）の v0.2.0 の認証（Stage G）が合格し、その記録が `main` に入っていること。入っていなければ、何も変えずに止め、README の表の状態を「条件待ち」のままにする。保守者が、BLASTX の範囲の拡大（`-query_loc` / `-subject_loc`）を LOSATX 計画の範囲の記録に書いていることも確かめる。

1. `git merge origin/main` で取り込み、衝突の解消を独立したコミットにする。取り込んだ後、S02〜S08 のゲートのうち、取り込みで変わったコードに関わるものを再実行する。
2. 基準：取り込んだ後のこの worktree で、LOSATX 計画の比較ゲートの出力の SHA-256 と、BLASTX の 1 つの fixture の性能（1 回の暖機と 3 回の計測、ネイティブと command-WASI）を取る。
3. BLASTX を `run_local` に通す。BLASTX は 10,002 文字の query バッチごとに検索と整形を繰り返す（`LOSAT/src/algorithm/blastx/web.rs`、`native.rs`、`query_setup.rs`）。この順序（prolog → バッチごとの検索と整形 → epilog）を保つ。formatter の中で書いている query ごとの警告（`LOSAT/src/algorithm/blastx/report.rs`、`LOSAT/src/report/pairwise.rs`）を `diagnostics` に移し、形式の数だけ重ならないようにする。outfmt 0 の既定のアラインメント数（`LOSAT/src/algorithm/blastx/args.rs` の `num_alignments`）は形式ごとに解決する。CLI の `run` と v1 の BLASTX の経路を、`run_local` の薄い層にする。
4. `hits` と観測者をつなぎ、`docs/web/verification_cells.tsv` の BLASTX の行を埋める。ABI v2 の `register` と `scan`（NCBI 型の解析器）を BLASTX で使えるようにし、V-ABI と V-BR に BLASTX の升目を加える。
5. `-query_loc` / `-subject_loc` を BLASTX に移植する（S11 と同じ手順。NCBI の範囲指定の意味は `c++/src/algo/blast/blastinput/blast_args.cpp:1946-1997`、`:2373`、`c++/src/algo/blast/blastinput/blast_fasta_input.cpp:433-459`）。`LOSAT/src/cli.rs` の BLASTX の範囲外の一覧から 2 つを外す。比較用の fixture を固定し、NCBI とバイト一致させる。
6. 既定以外のオプション（TD-13）：S08+ と同じ手順で、BLASTX の検索のオプションの既定以外の値を NCBI と比べ、NCBI と同じにするか、明示的に拒否する（`docs/losat_web_gui_sessions/session_s08p_e2e_protein_options.md`）。直した値と拒否する値を S12 の指示書（または、S12 が終わっていれば検索画面）に反映する。
7. 試験：LOSATX 計画の比較ゲートが第 2 項の基準と一致すること（移植した範囲指定を除く）。範囲指定の fixture で NCBI とバイト一致。BLASTX の全升目（0/6/7 × スレッド 1/2/4）の V-ABI と V-BR。v1 の reactor の検査。性能の中央値が基準の +5% 以内。`cargo fmt --check`・`clippy`・`cargo test --all-features`。
8. 独立監査を受ける。

記録は `docs/evidence/losat_web_sx/README.md`。1 セッションで終わらなければ、まとめ直し（SXa）と範囲指定（SXb）に分ける。

## S08 からの引き継ぎ（2026-10-02）

S08（E2b、[ゲート記録](../evidence/losat_web_e2b/README.md)）は、BLASTX と共有のコードで次を見つけ、DW-10 に従って BLASTX の振る舞いを変えずに残した。取り込んだ後に、LOSATX の比較ゲートで確かめて合わせる。

- `BLAST_LargeGapSumE` の計算の順：NCBI は `xsum -= num*log(...) - BLAST_LnFactorial(num)`（`c++/src/algo/blast/core/blast_stat.c:4560-4561`）で、差を引く。`LOSAT/src/stats/sum_statistics.rs` の `large_gap_sum_e`（BLASTX の `algorithm/blastx/even_gap.rs` が使う）は積を引いてから階乗を足すので、最後の桁が違い、HSP の順が変わることがある（TBLASTX では 700 の乱数の入力のうち 4 件）。TBLASTX は S08 で `ncbi_large_gap_sum_e` に移した。
- 説明の一覧：BLASTX は今も `write_subject_summary_table_with_sum_n`（最初の HSP と素の幅）を使う。NCBI の規則は `x_InitDeflineTable`（E2a §G.3、`write_blastn_description_table`）。
- `score_compare_match` の同点：TBLASTX は S08 で `Blast_InitHitListSortByScore` を移植した（`algorithm/tblastx/blast_gapalign.rs`）。BLASTX の `algorithm/blastx/seed.rs` にある注釈の写しは確かめていない。
- 入力の読み方：`algorithm/blastn/input.rs` の部品は program の名前を受け取るようになった（BLASTN・TBLASTX が使う）。
- SEG：NCBI の `s_SegSeq` は左の部分の再帰で見つけた区間のうち先頭の 1 つだけを残す（`c++/src/algo/blast/core/blast_seg.c:2086-2101` の `leftsegs->next = *segs`）。S08 で `LOSAT/src/utils/seg.rs` を NCBI と同じにし、BLASTX の 2 か所（`algorithm/blastx/query_setup.rs`、`algorithm/blastx/kappa.rs`）だけが `keeping_all_left_segments()` で前の振る舞いを保つ。取り込んだ後にこの 2 つを外し、LOSATX の比較ゲートで確かめる（S08 の調べでは、低複雑度の領域を持つ蛋白の 388 件の組で BLASTX の 344 件が NCBI と違い、直すと 0 件）。
- `-seg` の値（S08 の独立監査の後、`7fbbfad96`）：BLASTP・TBLASTN・TBLASTX は NCBI と同じに 1 つの空白で分け、0 以下の窓・locut・hicut を NCBI の既定のままにする（`c++/src/algo/blast/core/blast_filter.c:1147-1154`、`blastinput/value_parsers.rs` の `SegSpec::params`）。BLASTX は自分の `seg` の解析器（`algorithm/blastx/args.rs` の `BlastxSeg`）のままなので、取り込むときに同じにする。閉じた標準出力（`cli.rs` の `report_standard_output`）も BLASTX の報告はまだ使っていない。

## S08+ からの引き継ぎ（2026-10-04）

S08+（E2e、[ゲート記録](../evidence/losat_web_e2e/README.md)、[権威の記録](../evidence/losat_web_e2e/AUTHORITY.md)）は BLASTP・TBLASTN・TBLASTX を NCBI の app の層（`LOSAT/src/blastinput/app.rs`：`-outfmt` の読み方、NativeError の誤り、`-seg`・`-comp_based_stats` の読み方、環境の検査）と引数の読み方（`blastinput/value_parsers.rs` の `ncbi_*`）に移した。BLASTX は DW-10 に従い変えていない。第 6 項で BLASTX に同じことをするときの材料：

- `app.rs` の `parse_formatting_string`・`report_format`・`formatting_handler_check`・`xinclude_check`、`parse_seg_option`、`parse_comp_based_stats`（BLASTX も最初の文字の規則）、`check_unsupported_environment`・`check_old_fsc`・`check_query_split_environment(program, translated_query)`・`query_batch_size`、`stats/protein_options.rs`（行列の表、BEST の gap、推奨の threshold と window、`validate_protein_options`、NCBI の行列と gap の誤りの文言）、`cli.rs` の `is_ncbi_toolkit_arg`（`-help-full`・`-xmlhelp`・`-logfile`・`-conffile`・`-version-full*`。BLASTX は未適用）。
- 共有のコードで BLASTX の振る舞いを保った箇所：`tblastx/lookup/compressed.rs` の `build_blosum62_compressed_lookup` は overflow の bank を使い切ると今も panic する（BLASTP は `try_build_blosum62_compressed_lookup` で明示的に拒否。NCBI は heap を壊す）。threshold は今も Rust の飽和する変換（BLASTP・TBLASTX は x86-64 の `(Int4)` を再現、`core/blast_util.rs` の `ncbi_int4_from_double`）。`BlastpPairwiseReport.num_descriptions`・`num_alignments` は BLASTX では `usize::MAX`（BLASTX は自分の writer を使う）。`write_blastn_description_table` に `protein` の引数ができた（蛋白の題は `report/defline.rs` の `ncbi_protein_title`、`check_shown_protein_subject_title`）。BLASTX の説明の一覧と蛋白の subject の題（BLASTX の subject は蛋白）は、これに移せる。
- S08+ の第 1 回の独立監査（`docs/evidence/losat_web_e2e/audit/ROUND1.md`）で直した共有の部品の規則：query の警告の並べ替え（`query_warnings.rs` の `query_warning`、NCBI の `RemoveDuplicates`）、蛋白の配列の行で NCBI が警告して落とす byte を読むときに拒否する `check_protein_sequence_lines_of` と、50 文字の英字の後の末尾の空白の定義行を拒否する `check_protein_deflines_of`（`blastn/input.rs`。BLASTX の subject は蛋白）、query の batch ごとの題の警告（`prepend_batch_title_warnings`）、`-seg` の locut・hicut の `strtod` の規則（`app.rs` の `parse_seg_option`、`value_parsers.rs` の `strtod_reads_whole`）、`-dryrun` の明示的な拒否（`cli.rs` の `is_ncbi_toolkit_arg`。BLASTX には未適用）、対角の表の window の上限（`app.rs` の `check_diag_table_window`）。BLASTX の `report.rs:115` の `blast_seqalign.cpp:1484-1490` は 1485-1490 が正しい（BLASTP・TBLASTN は直した）。
- BLASTP で直した NCBI の規則で、BLASTX にも当てはまりうるもの：得点 0 の HSP は報告に出ない（`blast_seqalign.cpp:672-674`）、SEG で mask した query の列の小文字（BLASTX は既にする）、epilog の行列の名前は入力どおり・threshold は `%g`、`O` は X として検索し警告（`blast_setup_cxx.cpp:894-932`。BLASTX の subject は蛋白）、残基の無い subject の「Subject sequence contains no data」の警告、蛋白の題の 50 文字の警告（`fasta.cpp:1650-1673`）、query の警告を各 query の報告の前に書く（`QueryWarnings`）。

## S08+a・S08+b からの引き継ぎ（2026-10-04）

S08+a（[記録](../evidence/losat_web_e2e/s08pa/NOTES.md)）と S08+b は BLASTP・TBLASTN だけを変え、BLASTX は DW-10 に従い変えていない。BLASTX と共有する部品と、BLASTX にも当てはまりうる NCBI の規則：

- gapped DP の確保（TN-4、`20856865e`）：`LOSAT/src/algorithm/blastp/gapalign.rs` の `gap_dp_ensure_cells` は、`dp_mem_alloc` を NCBI の値のまま、実際に触れる cell を `MIN(dp_mem_alloc, N + 1)` にした（NCBI の `malloc` は触れない page を確保しない。`blast_gapalign.c:797-808,923-931`）。BLASTX も同じ経路を使い、出力は変わらない（巨大な `-xdrop_gap` の時間と記憶だけが変わる）。LOSATX を取り込むときに、BLASTX の巨大な X-drop の時間と記憶を NCBI と比べる。
- 窓の右端の traceback（TN-1 の残り、`68154c73e`）：TBLASTN は `blast_gapped_alignment_with_traceback_in_window` で、右の伸長の `N` を翻訳の窓の長さから計算し、DP の最後の列で窓の後ろの番兵を読む（`align_ex_protein` の `read_end_sentinel`。`blast_gapalign.c:429-432,563-578`）。BLASTP・BLASTX の従来の入口は `read_end_sentinel` が偽（最後の列は 0）。BLASTX の subject は蛋白なので窓の部分翻訳は無いが、traceback が subject の終わりに触れる組（`-comp_based_stats 0` と巨大な最終の X-drop、大きい `-evalue`）で NCBI と比べる。
- 同順位の並べ替え：NCBI の `qsort` は固定した NCBI BLAST+ の glibc（2.39）では安定な merge sort。S08+a（TN-5）が TBLASTN の、S08+b が BLASTP の、比べ方が全てを区別しない並べ替えを安定な `sort_by` にした。BLASTP・TBLASTN・BLASTX が共有する composition の窓の並べ替え `LOSAT/src/core/composition_adjustment/redo_alignment.rs`（1385・1403・1503・1522・1827 行付近）と BLASTX 自身の並べ替えは今も `sort_unstable_by`。NCBI の比べ方（`redo_alignment.c` の窓・HSP の比較）が全てを区別しないものは、同順位が出力に届くかを確かめ、届きうるなら安定にする。
- query の分割（RP-4）：NCBI は blastx の query も分ける（`local_blast.cpp:83-86` の chunk の大きさ 10002 塩基、`split_query_aux_priv.cpp:51-69` の重なり 297。`SplitQuery_ShouldSplit` は blastx で query が 1 つの batch だけを分ける、`split_query_aux_priv.cpp:72-97`）。LOSATX がこれを移植しているかを確かめる。BLASTP・TBLASTN の移植は `LOSAT/src/algorithm/common/protein_query_split.rs`・`blastp/query_split.rs`・`tblastn/stage_d_results.rs` の `merge_query_chunk`（`BlastHSPStreamMerge`・`Blast_HitListMerge`・`Blast_HSPListsMerge`）。
- `cli.rs` の `ncbi_preparsed_toolkit_word`（`6fa06cb6a`）：NCBI の argv の前処理（`ncbiapp.cpp:926-1001`）と同じく、option の値の位置の `-version`・`-dryrun` などを拒否する。BLASTX は未適用。
- compressed の lookup の走査（S08+b、`0e59e2f64`、監査 R2A-1）：word より短い subject は NCBI では一度だけ prime され、その word は subject の後ろの NULLB（compressed の文字が無い）を含むので何も見つけない（`aa_ungapped.c:496-500`）。LOSAT の `tblastx/lookup/compressed.rs` の `scan_subject` は word が範囲の外だと abort していたので、範囲の外を NULLB として読むようにした。BLASTX も同じ `scan_subject` を使う（`algorithm/blastx/seed.rs`）ので、1〜2 残基の subject と blastx-fast（compressed の lookup）で abort していた組は、NCBI と同じヒット無しになる。LOSATX の比較ゲートで確かめる。
- 末尾の `--`（S08+b、`13493774f`、監査 R2c-2）：NCBI の引数の解析は `--` を位置の引数の始まりとして読み（`ncbiargs.cpp:2866-2872`）、最後の `--` は何もしない。`cli.rs` は BLASTN・BLASTP・TBLASTN・TBLASTX だけで最後の `--` を落とす。BLASTX は未適用。
- BLASTP の one-hit の拡張の負の長さ（S08+b、`0e59e2f64`、監査 R2A-2）：`s_BlastAaExtendOneHit` の Int4 の長さは負になりうる（`aa_ungapped.c:1054,1083`）。BLASTX の one-hit（`-window_size 0`）の移植が同じ長さを `usize` にしていないかを確かめる。

## 終了・引き継ぎ

README の規則 8 に従う。SX を終えた後は、中断していた順番のセッションに戻る。
