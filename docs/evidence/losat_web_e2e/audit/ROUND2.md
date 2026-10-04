# E2e（S08+・S08+a・S08+b）独立監査 第 2 回

- 対象：S08+a を merge した後の本線。1 回目の監査は `491292327`（S08+b の 2.1〜2.3 の後。native `d8d18ec0…cb116`、指摘を直す前のゲートの run の native と同じ）、(a) の 2 回目は最後のエンジン `75cc8e565`（最後のゲートの native `4d6c036f…`）。監査の写しは `~/.cache/losat-web-gui-target/s08pb-audit/`（`src/`・`src2/`、`/mnt/c` の外）。
- 方式：Sonnet の監査役が 4 観点を並行して、読み取り専用で NCBI BLAST+ 2.17.0 と LOSAT を同じ引数で実行して比べた（stdout、stderr、終了コード、`2>&1`）。第 1 回の 5 つの報告を (a) BLASTP、(b) TBLASTN、(c) TBLASTX、(d) 報告と入力 の 4 つにまとめ、第 1 回の全ての再現の命令と harness（S08+a が `/mnt/c` の外に直したもの）を繰り返し、S08+a・S08+b の変更の周り（query の分割、巨大な X-drop、BLOSUM45 の同点、hard mask、窓の右端、toolkit の語、BLASTP の同順位）を新しく調べた。指示は [`round2/brief/`](round2/brief/)（`COMMON.md`、`ANGLE_A.md`〜`ANGLE_D.md`、2 回目の `ANGLE_A2.md`）、報告は [`round2/`](round2/)。
- 保守者の判断待ちの項目（D11、D12、D13、D14）は「pending the maintainer」として扱わせた。

## 1 回目の結果（`491292327`）

| 観点 | 結論 | 比べた数 | 新しい指摘 |
|---|---|---|---|
| (a) BLASTP（[`a_blastp.md`](round2/a_blastp.md)） | unsupported | 第 1 回の argv 約 1,443、sweep 約 3,600、無作為 約 10,700、分割 約 2,100、同順位 約 1,900 | R2A-1・R2A-2（abort 2 件） |
| (b) TBLASTN（[`b_tblastn.md`](round2/b_tblastn.md)） | supported | 約 8,000 | R2b-1、R2b-2 |
| (c) TBLASTX（[`c_tblastx.md`](round2/c_tblastx.md)） | supported | 第 1 回の 730、新しい約 900 | R2c-1、R2c-2、R2c-3 |
| (d) 報告と入力（[`d_reports_inputs.md`](round2/d_reports_inputs.md)） | supported | 報告 2,183、入力 6,390、再現 133、新しい 2,552 | R2D-1 |

第 1 回の指摘は 4 観点とも、一致、決めたとおりの明示的な拒否、受け入れのどれかで、差は無かった（S08+a が直した TN-1 の残り・TN-2・TN-4・TN-5・RP-4・`-out -version` を含む）。

## 指摘と対応

| ID | 重さ | 内容 | 対応 |
|---|---|---|---|
| R2A-1 | 高 | BLASTP の compressed の lookup（`-word_size 5`、blastp-fast）で、1〜2 残基の subject があると abort（終了コード 134）。NCBI は走査の範囲を広げて 1 度 prime し、word は subject の後ろの NULLB（compressed の文字が無い）を含むので何も見つけない（aa_ungapped.c:496-500、blast_aascan.c:264-286） | 直した（`0e59e2f64`）：prime の word の subject の外を NULLB として読む（`tblastx/lookup/compressed.rs` の `scan_subject`。BLASTX も同じ関数を使う、SX に記録）。fixture `e2e.blastp.compressed_short_subject` |
| R2A-2 | 中 | BLASTP の one-hit（`-window_size 0`）と低い `-threshold` で、ungapped の Int4 の長さが負になり abort。NCBI は負の長さを `BlastGetStartForGappedAlignment` に `Uint4` で渡し、最初の 11 文字の窓だけを数える（aa_ungapped.c:1054,1083、blast_gapalign.c:3394-3437） | 直した（`0e59e2f64`）：同じ読み方。窓が context か subject の終わりを越えるときは NULLB（`BLAST_SCORE_MIN`）を含むので、和は負で `q_start` になる。fixture `e2e.blastp.one_hit_negative_width`。2.6 MB の出力の realistic な例も NCBI と一致 |
| R2b-1 | 低〜中 | 有限の DBL_MAX の `-evalue`（`1.7976931348623157e308`）で NCBI の tblastn と blastp は SIGSEGV（blast_kappa.c:409,3687 の `best_evalue <= expect_value` が DBL_MAX で真、blast_hits.c:3266）。LOSAT は結果を出した（D12 は無限大だけを拒否） | 直した（`141e49567`）：D12 の拒否を DBL_MAX 以上に広げた（BLASTP・TBLASTN）。`1.7976931348623156e308` は実行し NCBI と一致 |
| R2b-2 | 低（情報） | 小さい `CHUNK_SIZE`（1、7）で LOSAT の TBLASTN は chunk ごとに NCBI の約 10 倍遅い（出力は同じ。chunk ごとの `query_set_setup`） | 残件（実験用の環境変数だけ。既定では chunk は数個） |
| R2c-1 | 中（資源だけ） | TBLASTX の `-threshold` が `INT_MIN` になる値（`+inf` など、完全な近傍）で、LOSAT の記憶が NCBI の最大 2.3 倍（30 kb の query で NCBI 3.8 GB、LOSAT 8.7 GB）。出力は同じ。NCBI は tblastx の長い query を 10002 塩基の chunk に分け、chunk ごとに lookup を作る（LOSAT の TBLASTX は分けない。NCBI の出力は `CHUNK_SIZE` 3000 と 999999 で同じ。`s08pb/txsplit/`） | 残件。最後のゲートの前のゲートの sweep で、この 2 行の LOSAT が並列の実行の中で記憶不足で止まった（単独では outfmt 0・6・7 とも NCBI とバイト一致、LOSAT 8.6 GB、NCBI 5.9 GB。`s08pb/thrinf/`）。最後のゲートの sweep は TBLASTX の並列を 3 にした |
| R2c-2 | 低 | 末尾の `--`：NCBI は位置の引数の始まりとして読み、最後なら何もしない（ncbiargs.cpp:2866-2872）。LOSAT は引数の誤り（終了コード 2） | 直した（`13493774f`）：BLASTN・BLASTP・TBLASTN・TBLASTX で最後の `--` を落とす。後ろに語があれば今までどおり引数の誤り（例外 1）。BLASTX は変えない |
| R2c-3 | 低（文言） | TX-3 の拒否の文言が、標準入力が普通のファイルでも「such as a pipe」と言う | 受け入れ（D10 の拒否の理由は正しい） |
| R2D-1 | 低（注釈） | `blastn/input.rs:418` の `fasta.cpp:966-979`（引用は 967 から） | 直した（`75cc8e565`） |

ゲートの v1 の WASI の行列で見つかったもの（監査の外）：S08+ が BLASTP を NCBI の app の層に移したため、web ABI v1 の BLASTP の振る舞いが E1d の記録から変わっていた（誤りの文言、2 つの誤りの順、知らない tabular の field を誤りにしない）。v1 は凍結（計画 TD-1）なので、v1 の経路で option を先に確かめ、書かない field を以前の文言で拒否するようにした（`efa445b19`。試験 `web_api::tests::blastp_web_pair_keeps_v1_error_order_and_fields`）。option の誤りの文言はエンジンの今の文言（ゲート記録の「v1」）。

## (a) の 2 回目（`75cc8e565`、native `4d6c036f…`）

[`round2/a2_blastp.md`](round2/a2_blastp.md)（指示 [`brief/ANGLE_A2.md`](round2/brief/ANGLE_A2.md)）。結論は unsupported（狭い）。R2A-1・R2A-2 は直った（短い subject の約 4,000 と負の one-hit の長さの約 17,000 の比較で abort なし）。D12 の DBL_MAX と最後の `--` は決めたとおり（4 つの program × 35 の argv）。第 1 回の harness の差は、R2A-1・R2A-2 の abort（sweep の 41）と `cfuzz`・`ofuzz` の 15 が 0 になった。新しい指摘が 2 件：

| ID | 重さ | 内容 | 対応 |
|---|---|---|---|
| R2A2-1 | 中 | `-task blastp-fast`（NCBI は HSP の chaining を有効にする、blast_options_handle.cpp:395-399）で明示的な低い `-threshold`（13 以下）のとき、LOSAT は NCBI が落とす短い HSP を 1 つ多く報告する（約 12,000 の比較のうち 68）。NCBI は init hit list が空の時だけ chaining を飛ばし（blast_gapalign.c:3707-3708）、HSP が 1 つの context もその得点で drop の検査（3628-3636）にかける。LOSAT は list が 1 つ以下で返っていた。第 1 回の実行ファイルも同じ（退行ではない） | 直した（`3b84ce6c4`）：空の時だけ返す。fixture `e2e.blastp.fast_lone_hsp`（outfmt 0）・`fast_lone_hsp_window`（outfmt 7）。監査の harness の blastp-fast の行と差のあった行 7,053 のうち 7,032 が NCBI と一致（R2A2-1 の 68 件を含む）、21 は NCBI が落ちる行（R2A2-2）（[`../s08pb/r2a2_1/`](../s08pb/r2a2_1/)） |
| R2A2-2 | 低（NCBI の不具合） | `-window_size 0` と `-threshold` 1〜3 の一部の one-hit の検索で NCBI が SIGSEGV（5,000 の無作為の組のうち 21）。`BlastGetStartForGappedAlignment`（blast_gapalign.c:3393-3416）が負の長さの 11 文字の窓を配列の外まで読み（valgrind：「0 bytes after a block」）、その値で行列を引く。LOSAT は配列の外を NULLB として読み（`blastp_get_start_for_gapped_alignment_int4_length`）、NCBI を `valgrind -q` の下で走らせた出力とバイト一致（21 件とも） | 保守者に諮る（NCBI の不具合の方針の「確かめられる妥当な結果を承認済みの例外にする」に当たる。推奨は承認済みの例外。それまで LOSAT の振る舞いは変えない） |

## (a) の 3 回目（`3b84ce6c4`、native `6f070575…`）

[`round2/a3_blastp.md`](round2/a3_blastp.md)（指示 [`brief/ANGLE_A3.md`](round2/brief/ANGLE_A3.md)。実行ファイルは最後のゲート `run-20261004T163746Z` の native と同じ bytes）。結論は **supported**。新しい指摘は無い。

- R2A2-1：直った。2 回目の再現と新しい fixture の組 160 件が一致。chaining の広い探索（`w4.py`：短い組、subject ごとに ungapped の整列が 1 つの組、1 つと複数の query（context）を混ぜた batch、100〜1,000 残基、e2e の入力と 300 の subject、query の分割、`-threshold` 1〜30、window・e-value・`-comp_based_stats`・`-max_target_seqs`・outfmt 0/6/7）21,500 件で一致 20,074、LOSAT の拒否 1,424（全て `-comp_based_stats 0`、§K の拒否）、NCBI が落ちる 2（R2A2-2）、差 0。2 回目の実行ファイルと違う 68 件は全て新しい実行ファイルが NCBI と一致。blastp-fast を使わない対照 1,100 件は 2 回目の実行ファイルとバイト一致。blastp-fast の差は約 12,000 件中 68 件から約 28,000 件中 0 件になった。
- R2A2-2：変わらない。NCBI が落ちる 23 件（`w2` の 21、`w4` の 2）で、LOSAT の出力は NCBI の `valgrind -q` の下の出力と 23 件とも一致（保守者の判断待ち、D15）。
- 2 回目の全ての harness を繰り返した：argv 1,443（差 5 は承認済みの `-help` 4 件と両方が timeout の `-task blastp-fast -threshold +inf`。13 行が負荷の下の時間で類が移った）、sweep 3,635（2 回目と行ごとに同じ）、無作為の比較 8,271（差 0）、`w1`・`w1` fast・`g1`・`w2`・`w3`・`w2L`・`g2`・`d12`・`dd` は R2A2-1 の行が一致になった以外は 2 回目と同じ。

## 結論

4 観点とも supported：(b) TBLASTN、(c) TBLASTX、(d) 報告と入力は 1 回目で、(a) BLASTP は 3 回目で（1 回目の R2A-1・R2A-2 と 2 回目の R2A2-1 を直した後）。比べた数は 4 観点で約 12 万（(a) 1 回目 約 19,700、2 回目 約 33,000、3 回目 約 51,000、(b) 約 8,000、(c) 約 2,200、(d) 約 11,000。回をまたぐ繰り返しを含む）。出力に関わる残りの差は無い。

保守者の判断待ち（supported の判断の外）：D11、D12（DBL_MAX 以上を含む）、D13、D14、D15（R2A2-2）。残件（資源・性能・文言）：R2c-1、R2b-2、R2c-3（[`../AUTHORITY.md`](../AUTHORITY.md) §N）。
