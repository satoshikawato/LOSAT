# LOSAT Web E2i（Session SD）ゲート記録

- 段階：E2i BLASTN の `-task dc-megablast` と `-task blastn-short`（[総合計画書](../../losat_web_gui_plan.md) §7 の SD、指示書 [SD](../../losat_web_gui_sessions/session_sd_e2i_blastn_dc_megablast.md)、計画 DW-18・DW-12・TD-1）
- ブランチ：`feature/losat-web-gui`。変更前はセッションの開始の `a92fa902f`（エンジンは S08b の最後の `24fcfe41b` と同じ。native の SHA-256 `2a46c2e4…`、S08b のゲートの成果物）。変更後は `eeeea4fb2`（最後のエンジンのコミット。移植は `90c5f0181` と `c45e57d85`、`eeeea4fb2` は監査の指摘の使われない定数と注釈で、実行ファイルは `c45e57d85` と同じ。native `487ac387…`）。セッションは 2026-10-03〜04（途中で WSL の `/mnt/c` の I/O の誤りで止まり、WSL の再起動の後に再開した）
- 状態：**完了（2026-10-04）**。計画 §7 の SD の行の完了条件をすべて満たした：固定した fixture のバイト一致（template の全組、2 つの入力）、両 task の sweep、BLASTN の既存のゲートの非退行、V-PERF、両 task の V-ABI、独立監査（第 1 回で 4 観点とも supported）。保守者の判断の要る項目は無い

## 完了条件（計画 §7 の SD の行）

| 条件 | 状態 | 根拠 |
|---|---|---|
| 固定した fixture で NCBI とバイト一致（outfmt 0/6/7、スレッド 1/2/4、dc-megablast は template の種類・長さ・word size の全組、blastn-short は短い query と既定値の上書き） | 満たした | outfmt 0/7 の manifest の 21 行（1・2・4 スレッドで差 0）、回帰 fixture 76 件（template の 3 種類 × 3 つの長さ × word 11・12 を 2 つの入力で。新しい入力では 18 の組の NCBI の出力がすべて違う）すべて差 0。変更前の実行ファイルはすべてで違い、壊した実行ファイルは新しい入力の該当する組すべてで違う（下の「fixture」） |
| NCBI が受け付ける両 task の組合せの sweep が、同じ拒否、バイト一致、明示的な拒否のどれか | 満たした | 得点 3 × 880、word size 256、入力 771、batch 3 × 120、slice 3 × 120 で、一致か同じ誤り（予期しない 0）。独立監査 (a)〜(c) の約 19000 組も、一致、同じ誤り、承認済みの例外、明示的な拒否だけ（下の「ゲート」と「独立監査」） |
| BLASTN の既存のゲートと fixture に退行なし | 満たした | capture 236 件（Gate A の BLASTN の行を含む）が S02 の基準とも変更前とも差 0、CI の速い検査の全件 236 件で失敗 0、E2c・E2f・E2g の sweep と fixture に差 0、v1 の WASI の行列で形式の失敗 0 |
| V-PERF の非退行（megablast と blastn）、dc-megablast の NCBI との時間の比 | 満たした | 12 組すべてが閾値以内（`blastn` native の 1 回目の ×1.095 は `--repeat 5` で ×0.981）、出力はすべて同じ。dc-megablast は NCBI の ×0.565〜×0.953（4 スレッドは ×0.701）、blastn-short は ×0.382〜×0.781（下の「V-PERF」） |
| 両 task の升目の V-ABI | 満たした | V-ABI full 600 実行で失敗した部分 0、dc-megablast と blastn-short の 18 の検索 × serial n1・threads n1/2/4 の 72 実行が native の CLI と一致し、凍結の NCBI のハッシュ 80/80。独立監査 (d) の 248 実行。`docs/web/verification_cells.tsv` の SD の 4 行を `checked` |
| 独立監査（4 観点とも supported） | 満たした | 第 1 回で 4 観点とも supported。低い指摘は `eeeea4fb2`・`c0849e604` で直し、ほかは S08+ に記録した（下の「独立監査」） |

## NCBI の経路の記録（指示書の 1.）

[`AUTHORITY.md`](AUTHORITY.md)（英語）。§A 両 task の option の値（`CBlastOptionsFactory::CreateTask`、`CDiscNucleotideOptionsHandle`、`CBlastNucleotideOptionsHandle` の blastn-short、CLI の引数が task の既定値に重なる順序）、§B 引数と option の検査（`CDiscontiguousMegablastArgs`、`CArgAllowIntegerSet` の 10 進の変換、`s_DiscWordOptionsValidate` と検査の順）、§C discontiguous の lookup の表（`s_FillDiscMBTable`、12 の template、`ComputeDiscontiguousIndex`、2 つめの表、`longest_chain`）、§D subject の走査（`s_MB_DiscWordScanSubject_1`・`_TwoTemplates_1`・`_11_18_1`・`_11_21_1` と選び方、mask された subject）、§E ungapped の拡張と対角線の構造（`s_BlastNaExtendDirect`、two-hit の窓 40、`hit_len_array`、`s_NuclUngappedExtend`）、§F gapped の段階、§G 報告（epilog の「Window for multiple hits」、gap extension 0 の式、megablast の文献）、§H このセッションで決めたこと、§I 棚卸し。ソースの引用はすべて固定 commit（598d8ae6）のファイルと行で、引用の文字列が行と一致することを [`check_authority.py`](check_authority.py) で確かめた（45 件、誤り 0）。

## 棚卸し（DW-12、指示書の 1.）

[`INVENTORY.tsv`](INVENTORY.tsv)（285 行、[`build_inventory.py`](build_inventory.py) が作る）。NCBI の dc-megablast と blastn-short の経路を 6 つの範囲（A アプリと引数、B C++ の option の層と検査、C core の準備と parameter、D lookup・走査・ungapped の拡張、E gapped の拡張・traceback・hit の保存、F 報告）に分け、読み取り専用の agent（sonnet）が移植の前に行を作り（`inventory/{A..F}.tsv`、規則は `inventory/COMMON.md`。列 `status_before` は SD の前の LOSAT、`a92fa902f`）、別の agent が移植の後の `90c5f0181` で行ごとの結果を NCBI との比較で確かめた（`inventory/result_{A..F}.tsv`、規則は `inventory/RESULT_COMMON.md`。足した行 `X1`… は棚卸しが落とした NCBI の関数）。megablast・blastn と共有で E2g が faithful とした行は、dc-megablast と blastn-short で通る分岐が同じかを確かめた。GAP の 2 行（範囲 A の 27・X1：`-template_length 0x12` を受け付けた）は `c45e57d85` で直した（`RESOLUTION`、列 `e2i_final`）。独立監査 (a) の後に 4 行（B-AU1〜AU4、すべて n/a）を足した（`AUDIT_ROWS`）。

| | faithful | ported | rejected | exception | n/a | 移植の前の残り |
|---|---|---|---|---|---|---|
| 移植の前（`status_before`） | 85 | — | 46 | — | 25 | unported 53、divergent 45。result の agent が足した 27 行と監査の後の 4 行は前の状態なし |
| 最後（`e2i_final`） | 102 | 99 | 49 | 3 | 32 | 0 |

`exception` の 3 行は承認済みの `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`（引数の依存と制約の誤りは終了コード 2 と clap の文言、`-num_threads` の警告）。`rejected` の 49 行は §H の 2. の明示的な拒否（task が内部で決める option と、`-task rmblastn`）。

## 移植（指示書の 2.）

| コミット | 内容 |
|---|---|
| `90c5f0181` | NCBI の経路を一括で移した（DW-12）。`disc_lookup.rs`（新規）：12 の template（`EDiscTemplateType`）、`s_GetDiscTemplateType`（coding_and_optimal の 2 つめは +1）、12 の `ComputeDiscontiguousIndex`、4 つの走査（`s_MB_DiscWordScanSubject_1`・`_TwoTemplates_1`・`_11_18_1`・`_11_21_1`）を NCBI の余りの塩基の switch と同じ phase の状態機械で。`lookup.rs`：`s_FillDiscMBTable`（1 から数える word の始まり、曖昧な文字での初期化、新しい順の鎖、共有の PV、`hashtable2`・`next_pos2`、`longest_chain`）と `BlastMBLookupTableNew` の discontiguous の分岐（scan_step 1、hashsize 4^word）。`coordination.rs`：task の既定値（dc-megablast は word 11・template coding/18・窓 40・ungapped の X-drop 20・DP・gap の X-drop 30/100・得点 2/−3・gap 5/2、blastn-short は word 7・reward 1・penalty −3・e-value 1000・filter なし・`-soft_masking` の既定の mask at hash）、`-template_*` はどの task にも効く。`run.rs`：template の長さを word length と lut word length にする `BlastNaWordFinder`、`s_BlastNaExtendDirect`、two-hit の窓（Delta 0、`hit_len_array` は u8、hash の挿入の窓 41）、対角線の配列の長さと offset、mask された subject でも discontiguous の走査。`scoring.rs`：`s_DiscWordOptionsValidate` を word size ≤ 100 の検査の後に（「word size must be either 11 or 12」、「Invalid lookup table type for discontiguous Mega BLAST」）。`args.rs`・`value_parsers.rs`：`-task` は 4 つ、`-template_type`・`-template_length`（互いに依存、`cli.rs` の拒否の一覧から外した）、`rmblastn` は明示的な拒否。`pairwise.rs`：epilog の「Window for multiple hits」、gap extension 0 の式は megablast と blastn の program だけ。`web_api.rs`：ABI v1 は megablast と blastn のまま（TD-1） |
| `af03129ee` | NCBI BLAST+ 2.17.0 を凍結した fixture（下の「fixture」）と入力（[`make_inputs.py`](make_inputs.py)） |
| `c45e57d85` | 棚卸しの結果の GAP の修正：`-template_length` を NCBI の `CArgAllowIntegerSet` と同じ 10 進の変換で確かめる（`0x12` は引数の誤り、`018`・`+18` は受け付ける）。`-task` の help に 4 つの task |
| `2bb2609fd` | アダプタ：`validate` の試験に両 task、`v_abi.js` の `describe` の検査（template の option、`-task` の help）、V-ABI の quick に `dc.default`・`short.default` |
| `b99f42493` | 記録：`AUTHORITY.md`、棚卸し、sweep と gate の script |
| `010d26ff1` | 独立監査の指示（`audit/`） |
| `eeeea4fb2` | 監査 (b) F1・F2、(c) F7：使われない `TWO_HIT_WINDOW` とその単体試験を外し、`scan_range` と `-penalty -32768` の注釈を直した（実行ファイルは変わらない） |
| `c0849e604` | 監査 (a)：棚卸しに 4 行と C-15 の証拠 |
| `f26a7d017` | template の区別の fixture `dc.div_*` 22 件と入力（下の「fixture」） |
| `4d3c9cc26` | 壊した実行ファイルの確かめ（`mutants/`）、E2c の `check_inputs.py` の期待値、`sd_postgate.sh` |

## fixture（指示書の 3.）

| 集まり | 件数 | 内容 |
|---|---|---|
| `LOSAT/tests/outfmt0_manifest.tsv` の `dc.*`・`short.*` | 21 | NCBI の outfmt 0 と 7 を凍結（`run_oracle.py`）。dc-megablast：既定、小さな query、optimal 16、coding_and_optimal 21、word 12 の coding 18（outfmt 7）、`-lcase_masking`、`-dust no`、得点の上書き、`-max_target_seqs`・`-max_hsps`、megablast に template、ヒット無し。blastn-short：既定（outfmt 0 と 7）、word 11、得点、`-evalue 10`、`-dust yes`、`-lcase_masking`、`-max_target_seqs`、ヒット無し |
| `LOSAT/tests/blastn_regression_fixtures.py` の `dc.*`・`short.*` | 54 | outfmt 6 を中心に stdout・stderr・終了コードを凍結。dc-megablast：word 11・12 × coding・optimal・coding_and_optimal × 16・18・21 の 18 組、megablast に template（word 11 と 28）、4 スレッド、2 スレッドの outfmt 7、小さな query の 2 つの template、小文字の島の subject と対照、予備の hit list（`-max_target_seqs 3` の outfmt 0）、query の batch と `-subject_besthit`、曖昧な文字、反復、query の分割（`CHUNK_SIZE`）、`BATCH_SIZE=1000` の outfmt 7、NCBI の誤り（word 13、blastn と blastn-short に template、gap 0/0）。blastn-short：既定、4 スレッド、word 4・16、得点と gap、`-evalue 1e3 -dust yes`、`-dust no` の outfmt 7、`-lcase_masking`、`-perc_identity` と `-subject_besthit`、batch、予備の hit list、曖昧な文字、分割、ヒット無しの outfmt 0、`-evalue 0` |
| 同じ script の `dc.div_*`（監査の後） | 22 | 18 の組で NCBI の出力がすべて違う入力の、18 の template の組、megablast の template、outfmt 0、outfmt 7、4 スレッド（下の「template の区別の fixture」）。回帰 fixture の SD の行は合わせて 76 件 |

入力（`LOSAT/tests/fasta/outfmt0/dc_*.fasta`・`short_*.fasta`）は [`make_inputs.py`](make_inputs.py) が repository の genome から決まった方法で切り出す（LC738874 の 2 つの窓に小文字・N・IUPAC、LC738875 の 3 つの窓と乱数の record、primer の長さの query 8 本と CA の反復）。凍結と確かめの記録は [`fixtures-20261003T103816Z/`](fixtures-20261003T103816Z/)（最初の 54 件と manifest の 21 行）と [`fixtures-20261003T164528Z/`](fixtures-20261003T164528Z/)（`dc.div_*` 22 件。既存の 161 件の NCBI の出力は凍結し直しても同じ bytes。ゲートの実行ファイルで 183 件差 0、変更前の実行ファイルで SD の 76 件が違う）。

区別の確認：変更前の実行ファイル（`2a46c2e4…`）は 54 件すべてと manifest の 21 行すべてで違う（task と template の option を拒否する）。

**template の区別の fixture（監査の後）。** 上の入力の 18 の template の組（word 11・12 × 3 種類 × 3 つの長さ）では、NCBI の出力は 9 通りしかなく（例：word 11 の coding_and_optimal 21 は coding 18 と同じ bytes）、走査の変種と 2 つめの template の区別が弱かった。sonnet の agent に、repository の genome の窓から、18 の組で NCBI の出力がすべて違う組を探させた（`~/.cache/losat-web-gui-target/sd-disc-fixtures/`、候補 27 組）。WSSV（AP027280）の 26 の窓と、遠い nimavirus（LC741431・LC738884）の 37 の窓（約 63〜75 % の一致の弱いヒットの周り、12.6 kb と 13.6 kb）で 18 通りになった（`make_inputs.py` の `DIV_*_WINDOWS`、`dc_div_query.fasta`・`dc_div_subject.fasta`）。LOSAT はこの組の 18 の組と megablast の template の行で、outfmt 0・6・7 とも NCBI とバイト一致（`-num_threads 2` は NCBI の警告だけが無い、承認済みの例外）。fixture `dc.div_*` 22 件（18 の組、megablast の template、outfmt 0、outfmt 7、4 スレッド）を足した。

壊した実行ファイルで確かめた（[`mutants/`](mutants/)：`build.sh` が `c0849e604` の写しに 1 つずつ変更を入れて組み、`compare.py` が両方の入力の 19 の行を NCBI と比べる）：

| 壊した所 | 前からの入力 | 新しい入力 |
|---|---|---|
| なし（ゲートの実行ファイル） | 19 行とも一致 | 19 行とも一致 |
| M1：2 つめの template の表を引かない | coding_and_optimal の 6 行のうち 5 行で違う（word 11 の 18 を見逃す） | 6 行すべてで違う |
| M2：dc-megablast の two-hit の窓を 0 に | dc-megablast の 18 行すべてで違う | 18 行すべてで違う |
| M3：11/21 の coding を 11/18 の走査で | `w11_coding_21` で違う | `w11_coding_21` で違う |
| M4：11/18・11/21 の coding を一般の走査で | 19 行とも一致 | 19 行とも一致 |

M4 は同じ出力になる変種である（NCBI の専用の走査は一般の走査と同じ word を同じ順に返す）。新しい入力は M1 の 6 行すべてを区別し、前からの入力が見逃した word 11 の coding_and_optimal 18 も捕まえる。

## sweep（指示書の 4.）

[`sweeps.sh`](sweeps.sh) が E2c の `scoring_sweep.py`・`word_size_sweep.py`（`--tasks dc-megablast,blastn-short`。既定の起動は変わらない）、[`check_inputs.py`](check_inputs.py)（E2g の入力の case に `.dc`・`.short` の変種を足す）、E2f の `batch_sweep.py`、E2c の `slice_sweep.py` を両 task で走らせる。移植の直後の実行（`90c5f0181`、sonnet の agent、[`sweeps/run-90c5f0181/SUMMARY.md`](sweeps/run-90c5f0181/SUMMARY.md)）で入力の case の期待値を決め（下の表の `E2I_EXPECT`）、ゲートと最後のコミットの確かめで同じ結果になった：

| sweep | 件数 | 結果 |
|---|---|---|
| 得点（outfmt 0・6・7） | 各 880（22 の reward・penalty × 10 の gap × 2 task） | 一致 260、同じ誤り 620 |
| word size（4〜28、e-value 10 と 1e5） | 256 | 差 0（dc-megablast の 11・12 以外は NCBI と同じ「word size must be either 11 or 12」） |
| 入力と query の batch（E2g の case と `.dc`・`.short` の変種） | 771 | 予期しない 0。E2g の期待値から変わったもの（`E2I_EXPECT`）：NCBI 自身が誤りにする 14 件（task の gap 5/2 が表に無い得点、word size、gap 0 と greedy 無し）、LOSAT が拒否しなくなった 6 件（greedy の gap の上限は megablast だけ、task の受け付け）、例外 2 の 1 件（blastn-short の e-value 1000 で句読点の題にヒットがあり NCBI が落ちる）。`audit15.small_evalue` の `.short` は、NCBI の blastn-short が全ゲノムで 14 GB を超えるので外した |
| query の batch（`batch_sweep.py`） | 3 × 120（dc-megablast、blastn-short、blastn-short `-evalue 1e-5`） | 差 0 |
| genome の切り出し（`slice_sweep.py`） | 3 × 120 | 差 0 |

## 決めたこと（推奨の案。AUTHORITY.md §H）

1. `-task megablast` に template（NCBI はどの task にも `-template_type`・`-template_length` を当てる）：NCBI と同じに、discontiguous の表と megablast の greedy の拡張と one-hit で動かす（fixture `dc.megablast_template`、`dc.megablast_w11_coding_18`）。blastn と blastn-short では NCBI の検査が拒否し、LOSAT も同じ文言。
2. 明示的な拒否のまま：`-task rmblastn`（行列の得点と masklevel）と、task の既定値が内部で決める `cli.rs` の `is_unported_blastn_arg` の option（`-window_size`、`-off_diagonal_range`、`-xdrop_ungap`、`-xdrop_gap`、`-xdrop_gap_final`、`-no_greedy`、`-ungapped`、`-soft_masking`、`-use_index`、`-index_name`、`-min_raw_gapped_score`）。E2c からの BLASTN の全 task の拒否で、SD で変えていない。
3. ABI v1（`web_api.rs`）は megablast と blastn のまま（TD-1 が v1 を凍結）。ABI v2 は両 task を `describe` し実行する。
4. NCBI の不具合に当たる新しい挙動は無かった（`PD-LOSAT-NCBI-DEFECTS` の新しい適用なし。blastn-short の句読点だけの題は既存の例外 2）。

## ゲート

`b99f42493`（移植と fixture と記録の最後のコミット。エンジンは `c45e57d85` まで）の全ゲート：run [`run-20261003T141051Z/`](run-20261003T141051Z/)（`head.txt`）。script は [`gates/sd_gates.sh`](gates/sd_gates.sh)（S08 の `s08_gates.sh` を写し、E2g の BLASTN の sweep と [`sweeps.sh`](sweeps.sh)、`check_authority.py` を足した。最初に Gate A の字句のパスを `ci_fast_regressions.py` の `stage_lexical_fixtures` で用意する）。成果物のハッシュは `artifacts.sha256`（native `487ac387…`、reactor `ee820996…`・`e8a3dd81…`。独立監査の実行ファイルと同じ）。変更前はセッションの開始の実行ファイル（`~/.cache/losat-web-gui-target/sd/bin/`、S08b のゲートの成果物、native `2a46c2e4…`）。BLASTN の Gate A は capture（S02 の 236 件、Gate A の BLASTN の行を含む）と CI の速い検査の凍結ハッシュで確かめる（E2c・E2f・E2g と同じ）。

| 検査 | 結果 |
|---|---|
| `cargo fmt --check`（LOSAT、adapter）、clippy `-D warnings` 4 構成と adapter の 3 構成 | すべて終了コード 0 |
| `cargo test --all-features`、adapter、wasm32 の web API の試験 | 882 件通過・失敗 0（無視 3）、adapter 7 件、web API 5 件 |
| pure-Rust の境界、`ci_fast_regressions.py` の単体試験、このセッションで足した行の NCBI の参照 | 通過、誤り 0（`verify-refs-session-added.txt`） |
| `check_authority.py` | 45 件、誤り 0 |
| outfmt 0 の fixture（`check_losat.py`、1・2・4 スレッド） | 95 件（SD の `dc.*`・`short.*` 21 行を含む）、差 0（3 つのスレッド数とも。承認済みの例外 2 行は分類どおり） |
| `run_oracle.py`（NCBI の凍結の確かめ）、`precheck_hits.py` | 差 0。`precheck` の差は承認済みの `-db_gencode` の例外の `tblastx.code4.0`・`tblastx.code4.7` だけ |
| BLASTN の回帰 fixture（`blastn_regression_fixtures.py check`） | 161 件、差 0。変更前の実行ファイルは SD の 54 件すべてで違う（`blastn-fixtures-before.tsv`） |
| TBLASTX の回帰 fixture、`CTOOLKIT_COMPATIBLE`（`ctoolkit_compare.py`） | 73 件差 0、279 実行差 0 |
| BLASTN の既存の sweep（megablast と blastn、E2c・E2f・E2g） | S06 の得点 80（一致 46・同じ誤り 34）、E2c の得点 outfmt 0/6/7 各 300 一致・580 同じ誤り、word size 256 差 0、slice 17 の組（40〜300 件ずつ）差 0、題 1023（一致 957・例外 2 が 66、予期しない 0）、batch 9 の組（60〜150 件ずつ）差 0、分割 220 件差 0（分割が NCBI の出力を変える 19 件を含む） |
| BLASTN の入力（E2g の `check_inputs.py`） | 300 件。予期しない 2 件 `audit6.task.dc_megablast`・`audit6.task.blastn_short` は、E2g の期待値（LOSAT の拒否）が SD の前のもので、今は NCBI とバイト一致（`same`）。E2c の期待値を `same` に直した（`4d3c9cc26`、下の最後のコミットの確かめで予期しない 0） |
| SD の sweep（[`sweeps.sh`](sweeps.sh)、dc-megablast と blastn-short） | 得点 outfmt 0/6/7 各 880（一致 260・同じ誤り 620）、word size 256 差 0（dc-megablast の 11・12 以外は NCBI と同じ誤り）、入力 771 件で予期しない 0（`.dc`・`.short` の変種、[`check_inputs.py`](check_inputs.py)）、batch 3 × 120 差 0、slice 3 × 120 差 0（`sd-sweeps/`） |
| CI の速い検査の全件（`ci_fast_regressions.py --all-cases`） | 236 件、失敗 0、許可した既知の不一致 1（`Sakai.MG1655.megablast`） |
| capture（`capture_outputs.py`） | 236 件が S02 の基準とも変更前とも差 0 |
| v1 の WASI の行列（`check_wasm_threading.py`） | 433 の記録、reactor の lifecycle の gate 通過、形式の失敗 0。reactor の記録の S05 との差 21 は S08b と同じファイルと内容（E2g の LOSAT の文言の差だけ）。`v1-requests` は E1d と一致 |
| V-ABI quick | 60 の実行が native の CLI と一致（`dc.default`・`short.default` を含む）、凍結ハッシュ 16/16 |
| V-ABI full | 150 の検索（BLASTN 77 のうち dc-megablast と blastn-short が 18、TBLASTX 36、TBLASTN 28、BLASTP 9）× 4（serial n1、threads n1・n2・n4）= 600 の実行、失敗した部分 0。各実行の 0/6/7 の stream が native の CLI と一致。凍結ハッシュ 860 件中 856 件一致、違う 4 件は既知の `Sakai.MG1655.megablast` outfmt 7。dc-megablast と blastn-short の 72 の実行は凍結の NCBI のハッシュ 80/80 |

### 最後のコミットの確かめ

ゲートの後のコミットは、エンジンの `eeeea4fb2`（使われない定数とその単体試験、注釈。監査の指摘）と、fixture・記録（`c0849e604`・`f26a7d017`・`4d3c9cc26`）だけである。最後のコミット `4d3c9cc26` から成果物を作り直し、変わりうる検査を繰り返した：run [`run-20261003T165813Z/`](run-20261003T165813Z/)（script [`gates/sd_postgate.sh`](gates/sd_postgate.sh)）。

- **成果物はゲートのものとバイト一致**（`artifacts.sha256`：native `487ac387…`、WASI の 4 つ、reactor `ee820996…`・`e8a3dd81…`）。`eeeea4fb2` は実行ファイルを変えない。したがって上のゲートの結果（V-ABI full、v1 の WASI の行列、sweep を含む）は最後のコミットの成果物の結果でもある。
- `cargo fmt --check`、clippy 4 構成と adapter の 3 構成：終了コード 0。`cargo test --all-features` 881 件通過・失敗 0（外した単体試験の 1 件だけ少ない）、adapter 7 件、web API 5 件、CI の Python の検査、このセッションで足した行の NCBI の参照：誤り 0。
- outfmt 0 の fixture 95 件（1・2・4 スレッド）と NCBI の凍結の確かめ：差 0。BLASTN の回帰 fixture 183 件（`dc.div_*` 22 件を含む）差 0、変更前の実行ファイルは SD の 76 件で違う。TBLASTX の回帰 fixture 73 件差 0。`check_authority.py` 45 件。
- E2g の `check_inputs.py` 300 件、予期しない 0（期待値を直した後）。SD の sweep はゲートと同じ結果。
- CI の速い検査の全件 236 件、失敗 0（許可した既知の不一致 `Sakai.MG1655.megablast` だけ）。capture 236 件が S02 の基準とも変更前とも差 0。V-ABI quick 60 の実行が native の CLI と一致、凍結ハッシュ 16/16。

## V-PERF

[`gates/sd_perf.sh`](gates/sd_perf.sh)、V-PERF の lock を取って（アプリ側の S09 は止まる）。結果は最後のコミットの確かめの run [`run-20261003T165813Z/`](run-20261003T165813Z/) にある。変更前はセッションの開始の成果物（S08b のゲートの native `2a46c2e4…` と WASI）、変更後は最後のコミットの成果物（ゲートのものとバイト一致）。case は E2c の [`perf_cases.py`](../losat_web_e2c/perf_cases.py) の BLASTN の 4 つ（`blastn`、`blastn-large`・`blastn-large-fmt0`（Gate A の EDL933 × Sakai の megablast）、`blastn-many`）× native・serial-WASI・threaded-WASI（4 スレッド）、変更前と変更後を交互に（暖機 1 回、3 回の計測）。計算機は静かだった（load 2 前後）。

| 計測 | 結果 |
|---|---|
| `perf-1`（`--repeat 3`、12 組） | 11 組が閾値（×1.05）以内（×0.416〜×1.030）、出力はすべて同じ。超えた 1 組：`blastn` native ×1.095（0.059 → 0.065 秒、起動の時間が大半の case） |
| `perf-2`（`blastn` を `--repeat 5`） | native ×0.981、serial-WASI ×1.013、threaded-WASI ×1.013。出力は同じ。`perf-1` の超過は計測の揺れ（同じ case の native の中央値が `perf-1` と `perf-2` で 0.06 秒と 0.12 秒） |
| dc-megablast と blastn-short の NCBI との比（[`perf_ncbi.py`](perf_ncbi.py)、`perf-ncbi.json`。NCBI BLAST+ 2.17.0 は subject を `makeblastdb` の `-db`、LOSAT は `-subject`。暖機 1 回、交互に 3 回、中央値） | EDL933 × Sakai の dc-megablast：LOSAT 8.73 秒・NCBI 9.16 秒（×0.953、11602 HSP）、4 スレッド LOSAT 6.43 秒（×0.701）。LC738874 × LC738875 の dc-megablast ×0.565、coding_and_optimal 21 ×0.655。blastn-short：LC738874 × LC738875 ×0.781（81055 HSP）、primer 8 本 × LC738875 ×0.382 |

megablast と blastn は非退行（閾値以内）。dc-megablast と blastn-short は新しい経路で、すべての case で NCBI より速い。

## 独立監査（指示書の 7.）

`ncbi_parity_auditor` の役割を、観点ごとに sonnet の agent 4 つが読み取り専用で並行して行った（共通の指示 `~/.cache/losat-web-gui-target/sd-audit/COMMON.md`、観点ごとの指示 `ANGLE_{A,B,C,D}.md`。写しは [`audit/`](audit/)）。対象はゲートの成果物の写し（`sd-audit/bin/`：native `487ac387…`、reactor `ee820996…`・`e8a3dd81…`。`b99f42493` のゲートの `artifacts.sha256` と同じ）。基準は `INVENTORY.tsv`。比べた全件の入力と出力は作業ディレクトリ `~/.cache/losat-web-gui-target/sd-audit/{a,b,c,d}/` にある（大きいので入れない）。(a) と (b) は `FINDINGS.md` の書き込みを道具に断られ、報告を文で返した（下の表に要約した。表は作業ディレクトリの `COVERAGE.tsv`・`CASES.tsv`）。

### 第 1 回（`b99f42493`）：4 観点とも supported

| 観点 | 結論 | 内容 |
|---|---|---|
| (a) 経路の網羅 | supported（低い文書の指摘 2 件、情報 3 件） | NCBI を callgrind の下で 106 の場面（dc 49、short 57）で実行し、届いた 1068 の関数を E2g の監査の 238 の実行と両方の棚卸しに突き合わせた。命令の番地の差（dc 52・short 54 と megablast 72・blastn 90 の基準、268 実行）で分岐ごとにも確かめた。blastn-short の word size 4〜23 × query 60 nt〜200 kb の格子（129 の比較）で `BlastChooseNaLookupTable` が選ぶ表・走査・拡張のすべての変種を通した。NCBI と LOSAT の比較 957 組のうち 927 がバイト一致、30 は承認済みの例外か明示的な拒否（引数の誤りの終了コード、`-num_threads` の警告、`-reward 0`、`BL2SEQ_LEGACY`）。`BATCH_SIZE`・`CHUNK_SIZE`・`OVERLAP_CHUNK_SIZE` の 6 つの組（24 実行）、`PRE_FETCH_SEQS_LIMIT`、`CTOOLKIT_COMPATIBLE` も一致。棚卸しの n/a・exception・GAP・rejected の 82 行の分類はすべて正しい。**指摘：** 1 dc と short が届く 13 の関数に棚卸しの行が無い（option の setter、`ClearFilterOptions` が消す filter、`SWindowMaskerOptions*`。値を運ぶ行はある）、2 C-15 の小さな表のあふれの代わりの表は「起きることが示されていない」とあるが、blastn-short で起きる（16 kb の 8 kb の単位の 2 回の繰り返しと word 7・8。出力は NCBI と一致）。情報：`BL2SEQ_LEGACY` の拒否（E2g）、`s_MBChooseScanSubject` の case 9 の `else` の無さ（出力に影響なし）、命令の差の方法の限界（inline された関数はソースと行で確かめた） |
| (b) 移植の忠実さ | supported（低い 2 件、情報 1 件） | NCBI と LOSAT の 15547 組（stdout・stderr・終了コード、170 件は `2>&1` も）：一致 15062（NCBI の stdout が空でないもの 10285）、承認済みの引数の誤りの例外 411、理由の成り立つ明示的な拒否 74、説明の無い差 0。12 の `DiscontigIndex_*` と `ComputeDiscontiguousIndex` を NCBI のヘッダから gcc で、`disc_lookup.rs` を rustc で組み、12 × 200000 の乱数の accumulator で 240 万行がバイト一致。`s_FillDiscMBTable`（1 から数える位置、曖昧な文字、新しい順の鎖、2 つめの表、`longest_chain`）、4 つの走査（`_11_18_1`・`_11_21_1` の 4 つの phase の式）、`BlastNaWordFinder` と mask された subject、two-hit の窓 40 の対角線の配列と hash（Delta 0、挿入の窓 41、`hit_len` の型、8000 の境界の前後）、task の既定値の組み立て（2352 件）、得点 1050 件、option の検査の順 576 件、query の chunk（1M・5M の境界、`CHUNK_SIZE` 504 件、`BATCH_SIZE` 216 件）、epilog（約 3000 件）を行ごとに読み、差分の実行で確かめた。**指摘：** F1 `coordination.rs` の `scan_range` の注釈が古い（blastn 4）、F2 `constants.rs` の `TWO_HIT_WINDOW`（0）が移植の後に使われていない。F3（情報）対角線の offset の `INT4_MAX/4` の初期化は 536 Mb を越える subject が要り、読むだけで確かめた |
| (c) 拒否の理由 | supported（中 1 件、低 7 件。どれも SD の前からの明示的な拒否か記録した差） | 2561 組。template の option は 4 task × 9 word size × 4 の組 × outfmt 0/6/7 の 432 実行で NCBI とバイト一致（3 つの NCBI の誤りの文言を含む）、448 の値の書き方で受け付けと拒否が一致（`018`・`+18` を受け、`0x12`・`18.0`・`1e1`・` 18`・`17`・`274`・`4294967314` を拒否）。2 つの問題が同時にあるときの誤りの順は 844 実行で一致。greedy の gap の上限は megablast だけ。`-task rmblastn` の拒否は NCBI の reward 0・penalty 0（行列の得点）で成り立つ。NCBI の `-help` の 71 の option はすべて受け付けるか LOSAT の文言で拒否（204 実行）。`-num_threads` と句読点だけの題は承認済みの例外の範囲どおり。**指摘：** F1（中）task が内部で決める option（`-window_size`・`-off_diagonal_range`・`-xdrop_*`・`-no_greedy`）の拒否は避けられる（NCBI は受け付けて出力が変わる。E2c からの全 task の拒否）、F2 NCBI が無視する option（`-use_index`・`-mt_mode`・`-soft_masking true`、outfmt 6/7 の表示の option）、F3 `-reward 0` と 0 でない penalty、F4 HTML の題の検査がヒットの無い subject にも及ぶ（S08 からの引き継ぎ）、F5 LOSAT の上限が NCBI の gap の表の誤りより先、F6 hit list のあふれの拒否がヒットの無い検索にも及ぶ、F7 文言（`-outfmt 13`・`14`・`19`、`-h`、`-penalty -32768` の注釈が primer の長さの query の NCBI の crash を書いていない）、F8 `-num_threads` 65535 以上 |
| (d) ABI v2 と構造化結果 | supported（低 3 件、情報 1 件、どれも SD の前からか記述どおり） | `describe` は serial と threaded で同じ JSON、`-task` の help に 4 つの task、`-template_type`・`-template_length` がある。`validate` は 175 の argv で CLI と同じ受け付けと文言（template の option、task、NCBI の `BLAST query/options error` の 4 つの検査）。62 の BLASTN の検索 × serial n1・threads n1/2/4 の 248 実行で、0・6・7 の stream が native の CLI とバイト一致、stream 3 も一致（7.6 MB・40541 HSP の outfmt 0、250 を越える subject、全 9 の template の組と word 11・12）。manifest の dc と short の 19 行で凍結の NCBI のハッシュ 84/84。HSP の記録を outfmt 0・6・7 と全項目で照合（約 77000 HSP、問題 0）。ABI v1 は凍結のとおり（432 実行が CLI と一致、新しい task と option は前の文言で失敗）。**指摘：** 最後の argv の語の空文字列（ABI の NUL の規則）、LOSAT の拒否の文言に `Error: ` が無い（§4 の記述どおり）、`-subject`・`-query` の無い argv を `validate` が受け付ける、`describe` に `-task`・`-template_type` の `choices` が無い（情報） |

### 第 1 回の指摘への対応

| 指摘 | 対応 |
|---|---|
| (b) F1・F2、(c) F7 の注釈 | `eeeea4fb2`：使われない `TWO_HIT_WINDOW` とその単体試験を外し、`scan_range` の注釈（全 task で `BLAST_SCAN_RANGE_NUCL` 0）と `-penalty -32768` の注釈（primer の長さの query で NCBI が落ちる）を直した。振る舞いは変わらない（下の「最後のコミットの確かめ」） |
| (a) 1・2 | `c0849e604`：棚卸しに B-AU1〜AU4 の 4 行（`GetDefaultsMode`、`SetMBTemplate*`、`ClearFilterOptions` が消す filter、`SWindowMaskerOptions*`。すべて n/a）、C-15 に blastn-short で代わりの表が起きることと出力の一致（`build_inventory.py` の `AUDIT_ROWS`・`AUDIT_NOTES`）。(a) が行が無いとした `CDiscNucleotideOptionsHandle` の構築子は B-3 にある |
| (c) F1〜F8 の残り、(d) の指摘 | どれも SD の前からの BLASTN の明示的な拒否か記録した差、または記述どおりで、完了条件（同じ拒否、バイト一致、明示的な拒否のどれか）を満たす。S08+ の指示書の「SD からの引き継ぎ」に記録した（(c) F1 の task が内部で決める option は、SD の後は値が `TaskConfig` にあるので、移す範囲と順を書いた） |
| fixture の区別（監査の後の自分の確かめ） | 下の「template の区別の fixture」 |

指摘への対応は振る舞いを変えないので、第 2 回の監査は行わない（4 観点とも第 1 回で supported）。

## アプリ側（S09・S12）への注意

- ABI v2 の BLASTN は `-task dc-megablast` と `-task blastn-short` を受け付け、形式 0・6・7 と HSP の記録を出す（`describe` の `-task` の help に 4 つの task、`-template_type`・`-template_length` の 2 つの parameter が増えた）。`describe` は `-task`・`-template_type` の `choices` を持たない（clap の独自の parser のため。選択肢は help の文にだけある）。検索画面の選択肢と既定値は [S12 の指示書](../../losat_web_gui_sessions/session_s12_w3_search_ui.md)の 8. に書いた（task ごとの既定値、template の依存と制約、`validate` が返す NCBI の文言）。
- `-template_type` と `-template_length` は片方だけだと clap の依存の誤り（終了コード 2）。dc-megablast と megablast では word size 11・12 だけ、blastn と blastn-short では template を NCBI の文言で拒否する。
- blastn-short の既定の e-value は 1000 なので、短い query でも大きな報告になりうる（`batch_sweep` の 1 つの case で outfmt 6 が 8321 行）。
- ABI v1 は凍結のまま（megablast と blastn だけ。新しい task と option は前の文言で失敗する）。
- verification の升目：`docs/web/verification_cells.tsv` の最後の 4 行（SD）を `checked` にした。

## 残件と引き継ぎ

- **`main` への PR**：[#112](https://github.com/satoshikawato/LOSAT/pull/112)（SD のコミット `90c5f0181`〜）。merge は保守者。CI の結果はこの記録の次のコミットで書く。
- **S08+**（[指示書](../../losat_web_gui_sessions/session_s08p_e2e_protein_options.md)の「SD からの引き継ぎ」、次のエンジン側のセッション）：SD の監査が記録した BLASTN の明示的な拒否と SD の前からの差。中心は、task が内部で決める option（`-window_size`・`-off_diagonal_range`・`-xdrop_*`・`-no_greedy`。監査 (c) の F1、中程度）で、SD の後は値が `TaskConfig` にあるので、移す範囲と順を書いた。ほかに NCBI が無視する option、`-reward 0`、誤りの順、hit list のあふれの拒否、文言、ABI v2 の `validate` の細部。S08 の引き継ぎ（BLASTN の HTML の題の検査をヒットのある subject に狭める、`-num_threads` 65535 以上）も残る。
- **S12**（[指示書](../../losat_web_gui_sessions/session_s12_w3_search_ui.md)の 8.）：両 task の既定値と template の option の規則を書き足した。
- **保守者の確認**：無し。V-PERF は閾値以内で、SD は保守者の判断の要る NCBI の不具合に当たらなかった。S08b の V-PERF の `tblastx-multi`（NCBI の batch の費用）は今も保守者の確認待ち（E2b のゲート記録）。
