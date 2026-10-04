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
- BLASTP で直した NCBI の規則で、BLASTX にも当てはまりうるもの：得点 0 の HSP は報告に出ない（`blast_seqalign.cpp:672-674`）、SEG で mask した query の列の小文字（BLASTX は既にする）、epilog の行列の名前は入力どおり・threshold は `%g`、`O` は X として検索し警告（`blast_setup_cxx.cpp:894-932`。BLASTX の subject は蛋白）、残基の無い subject の「Subject sequence contains no data」の警告、蛋白の題の 50 文字の警告（`fasta.cpp:1650-1673`）、query の警告を各 query の報告の前に書く（`QueryWarnings`）。

## 終了・引き継ぎ

README の規則 8 に従う。SX を終えた後は、中断していた順番のセッションに戻る。
