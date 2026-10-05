# LOSAT Web E2d（Session S11）ゲート記録

- 段階：E2d `-query_loc` / `-subject_loc`（[総合計画書](../../losat_web_gui_plan.md) §7 の S11、§4 の G9、計画 DW-9〜DW-11。指示書 [S11](../../losat_web_gui_sessions/session_s11_e2d_query_subject_loc.md)）
- ブランチ：`feature/losat-web-gui`。変更前は E2e の最後のコミット `4fab73fdb`（native の SHA-256 `6f070575…`、E2e の最後のゲートの成果物 `~/.cache/losat-web-gui-target/s08pb2-gate-*`）。変更後は最後のエンジンのコミット `095171eda`（最後のゲートもこのコミットで、native `【G:native】…`。記録は後のコミット）。2026-10-05
- 状態：**作業中（2026-10-05 に中断）**。移植、fixture、sweep、4 観点の監査と再監査 2 回は済み。残りは最後のゲート・Gate A・V-PERF と、この記録の【G】の欄（下の「中断の時点と再開の手順」）

## 中断の時点と再開の手順（2026-10-05）

保守者の指示で、最後のゲートの途中で中断した。エンジンの最後のコミットは `095171eda`（中断の時点で push 済み）。

- 中断したゲート（`095171eda`、記録には使わない。`~/.cache/losat-web-gui-target/s11/paused-run-20261005T074458Z/`）の途中の結果は全て通過：fmt・clippy の 4 構成とアダプタの 3 構成、`cargo test --all-features` 920 件・失敗 0、アダプタと wasm32 の web API の試験、pure-Rust の境界、このセッションで足した行の NCBI の参照の誤り 0、fixture 204 件 × 1・2・4 スレッドで差 0、`run_oracle.py` 204 件差 0、範囲の回帰 fixture 35 件差 0、TBLASTX 84・BLASTN 183 の回帰 fixture 差 0、`CTOOLKIT_COMPATIBLE` 200 行 600 実行で差 0、句読点の定義行と題の sweep で予期しない 0、BLASTN の入力 300 件で予期しない 0、CI の速い検査の全件 236 件で失敗 0（既知の Sakai 1）、範囲の sweep 1860 件で DIFF 0、E2e の BLASTP・TBLASTN の option の sweep は終了コード 0。TBLASTX の option の sweep、capture、v1 の行列、V-ABI は途中で止めた。
- その前の 2 つのゲート（`3d46b4cd8`、`bbceeb330`）は、後でエンジンを直したので記録に使わない（`superseded-run-*`）。`3d46b4cd8` の run の capture の前の検査も全て通過。開発中の capture（`LOSAT-wip2`）は 236 件が S02 の基準と差 0。
- 再開の手順：
  1. `echo docs/evidence/losat_web_e2d/run-$(date -u +%Y%m%dT%H%M%SZ) > ~/.cache/losat-web-gui-target/s11/rundir` で run の directory を決め、`bash docs/evidence/losat_web_e2d/gates/s11_gates.sh`（約 3.5 時間、V-ABI full が最も長い）。
  2. 並行して `bash docs/evidence/losat_web_e2d/gates/s11_gate_a.sh`（出力は `/tmp/claude-1000/` の下。`audit_tblastx_v010.py` は `/tmp` の外を拒否する。中断の前に `OUT` を直した）。
  3. ゲートの後、静かな計算機で `bash docs/evidence/losat_web_e2d/gates/s11_perf.sh`（`REPEAT=3`。閾値を超えた case は `REPEAT=5 SUFFIX=2 CASES=…` で測り直す）。
  4. この記録の【G】の欄、`docs/web/verification_cells.tsv` の S11 の 4 行（`run-<S11 final gate>` を run の directory に、`pending` を `checked` に）、計画の状態の行と §7 の S11 の行、セッション README の表、`evidence.sha256` を書き、コミットして push する。

## 完了条件（計画 §7 の S11 の行と指示書の 5.）

| 条件 | 状態 | 根拠 |
|---|---|---|
| 固定した fixture で outfmt 0/6/7 が NCBI とバイト一致 | 【G】 | outfmt 0/6/7 の manifest の E2d の 53 行（1・2・4 スレッドで差 0、NCBI の凍結の確かめも差 0）、範囲の回帰 fixture 35 件（stdout・stderr・終了コード、NCBI の誤りと誤りの順を含む）。変更前の LOSAT はどの行も「the NCBI BLAST+ option -query_loc is not supported by LOSAT's <PROGRAM>」で拒否する |
| 範囲の全ての書き方と端で NCBI と同じか明示的な拒否 | 【G】 | [`range_sweep.py`](range_sweep.py) の 1860 件で DIFF 0（下の「sweep」）。4 人の監査役の約 11 万の比較でも、一致、同じ誤り、承認済みの例外、明示的な拒否だけ（指摘 A-1 = D-1 は直した） |
| 既存の認証済みの出力に退行なし | 【G】 | fixture 151 件、capture 236 件が S02 の基準と差 0、Gate A、TBLASTX 84・BLASTN 183 の回帰 fixture、CI の速い検査の全件、E2e の option の sweep 3 つ、v1 の WASI の行列 |
| `run_local` と ABI v2 の `validate` から使える | 【G】 | V-ABI full（E2d の fixture の検索を含む）、監査 D の Web の harness（`validate`・`register`・`run` と CLI の約 18,400 件 × 3 形式） |
| `docs/web/verification_cells.tsv` に升目を足す | 満たした | E2d の 4 行（fixture、範囲の回帰 fixture、sweep、V-ABI）を `checked` で足した |
| 独立監査 | 【G】 | 4 観点（BLASTN、BLASTP、TBLASTN・TBLASTX、引数・端・アダプタ）と、指摘の修正の再監査 2 回（下の「独立監査」） |
| V-PERF の非退行 | 【G】 | 下の「V-PERF」 |

## 範囲の意味（要約。詳しくは [`AUTHORITY.md`](AUTHORITY.md)）

NCBI は、範囲を、その役割（query か subject）の入力が読む**すべてのレコード**に当て、各レコードの区間 `[from, min(to, len-1)]`（0 始まり・閉区間）の文字だけを検索し、結果をレコードの座標で報告する。

- 単位：1 始まり・両端を含む。BLASTN・TBLASTX の query と subject、TBLASTN の subject は塩基、BLASTP の query と subject、TBLASTN の query は残基。
- 文法（`ParseSequenceRange`）：`start-stop` の 2 つの 10 進数。0 以下、同じ値、逆向き、2 つでない部分は NCBI の誤り（`BLAST engine error: Invalid specification of query location (…)`、終了コード 3）。`NStr::StringToInt` が読めない部分（空白、`0x`、小数点、2^31 以上など）は LOSAT の明示的な拒否（判断 R1）。
- レコードの端：終わりが長さを越えると黙って切る。始まりが長さ + 1 は文字の無い区間（LOSAT の明示的な拒否、判断 R2。BLASTP の subject は NCBI と同じ警告）。始まりが長さ + 2 以上は、subject なら `BLAST query/options error: Invalid from coordinate (greater than sequence length)`（終了コード 1）、query ならそのレコードを黙って飛ばす（`Query_<n>` の番号は飛ばしたレコードも数え、batch の全てを飛ばすと `BLAST engine error: Empty CBlastQueryVector`、終了コード 3）。
- 検索：統計・DUST・SEG・lowercase の mask・lookup・query の分割・翻訳の frame は区間のもの。e-value と検索空間は区間の長さから。
- 報告：座標はレコードの座標（両 strand・全 frame）。`Length=`・`slen` はレコードの長さ、BLASTP の `qlen` は入力どおりの範囲の長さ（端で切らない）、TBLASTN・TBLASTX の frame はレコードから計算し直した値、database の合計は区間の長さの和。
- NCBI の決まった振る舞いとして再現したもの（`PD-LOSAT-NCBI-DEFECTS` の「決まった結果は再現」）：分割した範囲付きの query の chunk の mask、BLASTN の subject の lowercase の表示、TBLASTX の mask の表示の frame（AUTHORITY §E）。

## NCBI の経路の記録と棚卸し

[`AUTHORITY.md`](AUTHORITY.md)（§A 引数、§B 入力のレコード、§C 検索、§D 報告、§E 再現する NCBI の振る舞い、§F 判断、§G 範囲の外）。ソースは固定 commit `598d8ae6`、振る舞いの確認は NCBI BLAST+ 2.17.0。コードの中の引用は `verify_refs.py` で確かめた（このセッションで足した行の誤り 0、下の「ゲート」）。

[`INVENTORY.tsv`](INVENTORY.tsv)（248 行）。範囲の経路の NCBI の関数を 5 つの範囲（LA 引数 52、QN 核酸の query 33、QP 蛋白の query 48、SB subject 68、RF 結果と報告 47）に分け、読み取り専用の agent（sonnet）が移植の前に行を作った（`inventory/{範囲}.tsv`、規則は [`inventory/COMMON.md`](inventory/COMMON.md)、範囲ごとの指示と覚え書き・oracle の観察は `inventory/{範囲}_brief.md`・`_notes.md`）。列 `status_before` は S11 の前の LOSAT（`4fab73fdb`）。

| | faithful | divergent | unported | rejected | exception | n/a |
|---|---|---|---|---|---|---|
| 移植の前（`status_before`） | 82 | 42 | 80 | 26 | 2 | 16 |

棚卸しの divergent・unported の行を一括で移した（共有の部品 `LOSAT/src/blastinput/seq_range.rs` と 4 つの program の入口・batch・報告）。rejected は判断 R1・R2・R4 と範囲の外（`-db`・`-remote`・`-import_search_strategy`・BLASTX）。移植の後の振る舞いは、fixture、sweep と 4 人の監査で確かめた。

## 移植

| コミット | 内容 |
|---|---|
| `46d93c06d` | 4 つの program への移植。`seq_range.rs`（`parse_sequence_range`、`record_interval`、`cut_subjects`・`cut_queries`、`QueryInput` の batch と飛ばしたレコード、`check_no_empty_interval`、`record_frame`）。BLASTN `run.rs`（範囲の読みの順、batch、`restrict_masks` の chunk の mask、`shift_to_records`、表示の mask）、BLASTP `blast_engine.rs`（`SubjectSet`・`BlastpRanges`、`qlen`・`slen`、epilog）、TBLASTN（`ReportRanges`、`record_coordinates`、frame）、TBLASTX（`QueryRuns`、`shift_to_records`、mask の frame）、`cli.rs` の拒否の表から 2 つを除いた（BLASTX は SX まで拒否）、アダプタの `validate` |
| `d67175bd3` | NCBI BLAST+ 2.17.0 で凍結した fixture（下の「fixture」） |
| `a262c8726` | `verify_refs` で見つけた NCBI の参照の行の範囲 8 か所 |
| `cc65da7b6` | 記録：`AUTHORITY.md`、棚卸し、sweep、ゲートの script |
| `3d46b4cd8` | S12 の指示書に範囲の意味（単位、文法と誤り、レコードの端、結果、`validate` と `run`） |
| `884009a64` | 監査 A-1（= D-1）：option の値の UTF-8 でないバイトを明示的に拒否（4 つの program） |
| `fc9970539` | A-1 の修正の再監査の R-1・R-2：拒否を parser の後に（`-help` と parser の誤りが先、NCBI と同じ順）、`-query=` などの書き方のファイル名はバイトのまま（ファイル名の検査が働く） |
| `bbceeb330` | 監査の報告、`docs/web/abi_v2.md` の `validate` と範囲、S12 の指示書の `validate` と CLI の誤りの順 |
| `095171eda` | 再監査 2 の F-1：実数の引数の読み（`ncbi_string_to_double`）が、4 バイト目をまたぐ複数バイトの文字で panic していた（`-evalue -xx€` など。E2e からの不具合で、`fc9970539` が UTF-8 でない値からも届くようにした）。バイトで比べる |

## fixture

| 集まり | 件数 | 内容 |
|---|---|---|
| `LOSAT/tests/outfmt0_manifest.tsv` の `e2d.*` | 53（BLASTN 20、BLASTP 15、TBLASTN 9、TBLASTX 9。outfmt 0 が 29、6 が 8（BLASTP の custom の field 2 を含む）、7 が 16（同 1 を含む）） | NCBI の stdout を凍結（`run_oracle.py`）。BLASTN：query・subject・両方の範囲、minus の strand、複数のレコードと飛ばすレコード、終わりが長さを越える範囲、subject の範囲と複数の subject、subject と query の lowercase の表示（区間の座標の検査）、DUST の端、dc-megablast、batch の境の飛ばすレコード、`Query_<n>` の番号。BLASTP：`qlen`・`slen` の field、終わりが長さを越える `qlen`、飛ばすレコード、subject の範囲、文字の無い subject の区間、SEG、O の警告、query の分割（区間 25000 残基）、blastp-fast。TBLASTN：subject の 6 つの frame、範囲の端、両方の範囲、lowercase、複数の subject、分割した query の chunk の mask。TBLASTX：両方の範囲、SEG の mask の表示（範囲の始まり 5・31・61）、query の範囲、複数の subject、飛ばすレコード |
| `LOSAT/tests/range_regression_fixtures.py` | 35 | NCBI の stdout・stderr・終了コードを凍結（`LOSAT/tests/fixtures/range_regression/`）。文法の誤りと他の誤りとの順（subject の範囲が先、`-dust` の後、`-evalue 0` の前）、subject の範囲がレコードを越える誤り（4 つの program）、全ての query を飛ばす（`Empty CBlastQueryVector`）、最後の batch を飛ばす（前の batch の報告、epilog なし）、分割した範囲付きの query の mask（`CHUNK_SIZE=20000`） |

入力は `LOSAT/tests/fasta/outfmt0/e2d_*`（約 340 KB）。fixture の数は 204（E2e の終わりの 151 に 53）。

## sweep

[`range_sweep.py`](range_sweep.py) が、4 つの program で、範囲の 33 の書き方（`ParseSequenceRange` と `StringToInt` の読み方の全ての枝）× query・subject、レコードの端（始まり 1・L−1・L・L+1・L+2 × 終わり L−1・L・L+1・2147483647）、2 つの範囲の組、E2d の複数のレコードの入力（短いレコード、飛ばすレコード、subject の誤り）を outfmt 0・6・7 で NCBI BLAST+ 2.17.0 と比べる（stdout、stderr、終了コード）。

| | 件数 | 変更前（`4fab73fdb`） | 最後のゲート |
|---|---|---|---|
| BLASTN | 540 | 全て LOSAT の拒否（option が未移植） | 【G】 |
| BLASTP | 480 | 同上 | 【G】 |
| TBLASTN | 420 | 同上 | 【G】 |
| TBLASTX | 420 | 同上 | 【G】 |

LOSAT の拒否は R1（`StringToInt` が読めない部分）と R2（文字の無い区間）だけ。

## 判断（推奨の案で進め、記録した。保守者の常の指示）

[`AUTHORITY.md`](AUTHORITY.md) §F。

- **R1**：`StringToInt` が読めない範囲の部分（NCBI は build のパスを含む文言で終了コード 255）は明示的な拒否（計画 TD-15 と同じ扱い）。
- **R2**：文字の無い区間（始まりがレコードの長さ + 1）は明示的な拒否（query は 4 つの program、subject は BLASTN・TBLASTN・TBLASTX。BLASTP の subject は既存の空の subject の経路で NCBI と一致）。
- **R3**：LOSAT は subject のファイル全体を読んでから LOSAT の拒否を出す（NCBI はレコードを 1 つずつ読み、範囲の誤りで止まる）。両方が NCBI の誤りを出す組では順が同じ。違いは LOSAT の拒否の側だけ。
- **R4**：`-strand`（BLASTN・TBLASTX）は拒否のまま（範囲と別の option）。
- **A-1**（監査）：option の値の UTF-8 でないバイト（NCBI は値のバイトをそのまま読む）は、引数が parse できた後に明示的な拒否（「the value of -X is not UTF-8; NCBI BLAST+ reads the bytes of an option's value, which is not supported by LOSAT's <PROGRAM>」、終了コード 2）。`-help` と parser の誤り（承認済みの例外 1。数の option は、U+FFFD に置き換えた値を parser が読めないので parser の誤り。NCBI も USAGE）が先。`-query`・`-subject`・`-out` のファイル名は以前からのファイル名の検査。Web の ABI は UTF-8 の argv だけを受けるので、この拒否は CLI だけに現れる。

## ゲート

【G:ゲート】

## Gate A

【G:Gate A】

## V-PERF

【G:V-PERF】

## 独立監査

sonnet の監査役 4 人が読み取り専用で並行して、監査の対象の commit `a262c8726`（エンジンは移植の最後）の実行ファイルを NCBI BLAST+ 2.17.0 と比べた（指示は [`audit/brief/`](audit/brief/)、報告は [`audit/`](audit/)）。主張：4 つの program の `-query Q -subject S` に範囲を付けた検索の stdout・stderr・終了コードが、outfmt 0・6・7 と全てのスレッド数で NCBI と同じ。例外は R1・R2、承認済みの例外、既存の明示的な拒否。

| 観点 | 比較の数 | 結論 | 指摘と対応 |
|---|---|---|---|
| (a) BLASTN | 約 1,500 | unsupported（A-1 だけ） | A-1：範囲の値の UTF-8 でないバイトで LOSAT が parser の誤り（clap、終了コード 2）。NCBI は範囲の誤り（終了コード 3 か 255）。`884009a64`・`fc9970539` で明示的な拒否に |
| (b) BLASTP | 約 12,550 | supported | なし |
| (c) TBLASTN・TBLASTX | 約 10,500 | supported | なし（`-db_gencode` 2・4・12 は範囲の無い場合と同じ差で、承認済みの例外） |
| (d) 引数・端・アダプタ | CLI 約 77,800、Web 約 18,400 × 3 形式 | unsupported（D-1 = A-1 だけ） | D-1 は A-1 と同じ。O-1：`validate` は引数だけを調べるので、引数とレコードの両方に誤りがあると、CLI はレコードの誤り（subject を先に読む）、`validate` は引数の誤りを先に出す。設計どおり（`validate` はレコードを持たない）として `docs/web/abi_v2.md` と S12 の指示書に書いた。O-2：`register` は残基の無いレコードと空の定義行を `run` の前に拒否する（以前からの拒否の類）。`run` の経路は CLI とバイト一致（`register` が通った全件） |

A-1 の修正の再監査：

| 再監査 | 対象 | 比較の数 | 結論 | 指摘と対応 |
|---|---|---|---|---|
| 1（[`audit/REAUDIT1.md`](audit/REAUDIT1.md)） | `884009a64` | 17,662（約 38,000 の実行） | supported | R-1：`-help` の後の UTF-8 でない値で help を出さなくなった（承認済みの例外 1 の外）。R-2：`-query=` などの書き方のファイル名が値の拒否になった。どちらも `fc9970539` で直した |
| 2（[`audit/REAUDIT2.md`](audit/REAUDIT2.md)） | `fc9970539` | 49,271（約 187,000 の実行） | unsupported（F-1 だけ） | R-1・R-2 は解消。F-1：`-evalue`・`-perc_identity`・`-threshold`・`-xdrop_gap`・`-xdrop_gap_final` の値が、U+FFFD に置き換えた後に符号か `in`・`na` で始まり、4 バイト目を複数バイトの文字がまたぐと panic（終了コード 134、`value_parsers.rs:273`）。有効な UTF-8 でも E2e から（`-evalue -xx€`）。`095171eda` で直した |
| 2 の再生（[`audit/REAUDIT2_REPLAY.md`](audit/REAUDIT2_REPLAY.md)） | `095171eda` | 91,072（U+FFFD の対を含む） | supported | 再監査 2 の全ての command line を `095171eda` の実行ファイルで再生し、`fc9970539` の結果と比べた。違うのは F-1 の 694 件だけで、全て parser の誤り（「expected a decimal number (other forms, which NCBI BLAST+ may read, are not supported by LOSAT's <PROGRAM>)」、終了コード 2。U+FFFD に置き換えた同じ命令と同じ）。ほかの 6 件は harness の作業 directory のファイルの状態（前の `-out` が書いたファイル）によるもので、同じ状態では差 0。panic は 0 |

## 保守者に諮ること

推奨の案で進め、記録した。

1. **R1**（`StringToInt` が読めない範囲の部分）：明示的な拒否。推奨：このまま（TD-15 と同じ）。
2. **R2**（文字の無い区間）：明示的な拒否。推奨：このまま（実用が無く、LOSAT は文字の無いレコードも拒否している）。
3. **A-1**（option の値の UTF-8 でないバイト）：parse の後の明示的な拒否。推奨：このまま。
4. **O-1**（`validate` と CLI の誤りの順）：引数とレコードの両方に誤りがある入力で、Web は引数の誤りを先に出す。推奨：ABI v2 の設計として受け入れる（`validate` はレコードを持たない。`run` の結果は CLI と同じ）。
5. V-PERF：【G】

## アプリ側（S12）への注意

- 範囲の意味（単位、文法と誤り、端、結果、`validate` と `run`）は S12 の指示書の 3. に書いた（`3d46b4cd8`、`bbceeb330`）。推奨：フォームは範囲を 1〜レコードの長さに制限する（BLASTP の `qlen` は入力どおりの範囲の長さを出すため、また端の誤りを画面で起こさないため）。
- ABI v2 の `validate` は範囲の文法を CLI と同じ順と文言で返す。レコードとの関係は `run` が返す（`docs/web/abi_v2.md`）。
- web ABI v1 は範囲を拒否のまま（TD-1）。
- verification の升目：`docs/web/verification_cells.tsv` に E2d の 4 行を `checked` で足した。

## 残件と引き継ぎ

- 【G:残件】
- **SX**：BLASTX の `-query_loc` / `-subject_loc`（`is_unported_blastx_arg` は変えていない）。`seq_range.rs` の部品と、この記録の §B〜§D がそのまま使える。
- **次のセッション**：[S12 — 検索画面](../../losat_web_gui_sessions/session_s12_w3_search_ui.md)。
