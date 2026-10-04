# E2e（S08+）独立監査 第 1 回

- 対象：`f7ce8199f`（S08+ の移植と記録の最後のコミット。監査の実行ファイルは同じ源の `s08p-native-b`）。監査の写しは `~/.cache/losat-web-gui-target/s08p/audit/src`（`/mnt/c` の外。WSL の 9p を避けるため）。
- 方式：5 つの観点を Sonnet の監査役 5 人が並行して、NCBI BLAST+ 2.17.0 と LOSAT を同じ引数で実行して比べた（stdout、stderr、終了コード、`2>&1`）。観点ごとの報告（再現の命令、照合済みの一覧）は [`round1/`](round1/)：BLASTP の app の流れと option（[`blastp.md`](round1/blastp.md)、約 1000 の argv）、TBLASTN（[`tblastn.md`](round1/tblastn.md)、約 1150 の 2 option の組と 700 の無作為の組）、TBLASTX（[`tblastx.md`](round1/tblastx.md)、約 1000）、BLASTP・TBLASTN の報告（[`reports.md`](round1/reports.md)、約 1000）、入力（[`inputs.md`](round1/inputs.md)、約 2100）。
- 結果：指摘 53 件（重複を除くと 46 件）。直したもの、明示的な拒否にしたもの（保守者の確認を待つ判断 D11・D12 を含む）、記録して受け入れたもの、次のセッション（S08+b）に残したもの、に分けた。直した組合せは NCBI で固定した fixture にした（`LOSAT/tests/outfmt0_manifest.tsv` の `e2e.*` の 9 行と、入力を 1 つの query にした `e2e.blastp.threshold_inf`）。

## 指摘と対応

| ID | 重さ | 内容 | 対応 |
|---|---|---|---|
| IN-8 | 高 | BLASTP：O と X の対で同一性を数える（検索は O を X として読むが、NCBI の報告は読んだ文字で比べる） | 直した：最終の HSP で O のある対の同一性・不一致・positive を読んだ文字で数え直す（`blast_engine.rs` の `recount_pyrrolysine_identities`、tabular.cpp:1006-1020）。fixture `e2e.blastp.o_identity`・`o_identity0` |
| IN-2・BP-1 | 高 | BLASTP：全 subject に残基が無いと outfmt 0 の有効な探索空間が query の長さ（NCBI は 0） | 直した（blast_setup.c:729-732）。fixture `e2e.blastp.empty_only` |
| IN-4 | 高 | 最初のレコードが `>` だけの subject を `bio` がレコード無しと読み、NCBI の空の subject の誤りになる（3 つの program） | 明示的な拒否：遅らせた定義行の検査が先に拒否する（`Query is Empty!` は NCBI と同じ） |
| IN-9・TX-3 | 高 | TBLASTX：`-query - -subject -` で `Query is Empty!`（BLASTP・TBLASTN は明示的な拒否） | 直した：BLASTP・TBLASTN と同じく位置の無い標準入力として明示的に拒否 |
| IN-10・IN-14・BP-13 | 高 | アダプタの `register` が BLASTP・TBLASTN の入力の検査をしない。BLASTP の `run_local` に残基の検査と空の query・subject の分岐が無い | 直した：`register` が蛋白の入力に CLI と同じ検査（`check_protein_sequence_lines_of`、`check_protein_input_of`、query の残基の無いレコード）、TBLASTN の subject に核酸の検査をする。`run_local` に空の query・subject と残基の検査 |
| IN-1・RP-1 | 中 | O を含む無効な query の 2 つの警告の順（NCBI は Karlin の警告が先） | 直した：NCBI の `RemoveDuplicates`（blast_aux.cpp:1043-1054）の並べ替え（同じ error id と重さなので文字列の順）。fixture `e2e.blastp.o_invalid_order`・`e2e.tblastn.o_invalid_order`（stderr も固定） |
| IN-3 | 中 | 蛋白の定義行の末尾の空白の前に 50 文字の英字があると LOSAT だけが警告 | 明示的な拒否（蛋白の 50 文字の規則。核酸の 20 文字の規則は蛋白に当てない） |
| IN-5・IN-7b | 中 | 蛋白の配列の行の NBSP などと数字：NCBI は読むときに警告して落とす | 明示的な拒否を読むときに行う（`check_protein_sequence_lines_of`。NCBI が黙って飛ばす空白・`;`・`!`・`#` の行は今までどおり後で） |
| IN-6・RP-5 | 中 | outfmt 0 の `2>&1` で、蛋白の query の題の警告が prolog の前（NCBI は query の batch を読むとき） | 直した：各 batch の最初の query の報告の前に書く（`prepend_batch_title_warnings`、blastp_app.cpp:253-259）。`BATCH_SIZE` で batch を分けた場合も NCBI とバイト一致 |
| IN-7a | 中 | 空の query で LOSAT だけが subject の「no data」の警告を書く | 直した：`Query is Empty!` の後（NCBI の formatter の作成時） |
| IN-11 | 中 | アダプタの `validate` が BLASTP・TBLASTX の `-num_threads` を確かめない | 直した：`check_options` が `validate_threads` を呼ぶ（abi_v2.md：host は実行する argv を検証する） |
| TX-1 | 高 | `-window_size` が 2^31 − query の長さ以上：NCBI は `Int4` が回り込み hit 無し、LOSAT は hit を出す | 明示的な拒否（判断 D11：query の長さ＋window が 2^30 を超える値。BLASTP も） |
| TX-2 | 中 | 同じく (2^30 − 長さ, 2^31 − 長さ) で NCBI も LOSAT も終わらない（D8 が未実装だった） | 明示的な拒否（D8 を実装、`app.rs` の `check_diag_table_window`） |
| TX-4・BP-6 | 中 | `-seg` の locut・hicut の `2.5x`・`1e`・`.` などに LOSAT の拒否の文言 | 直した：NCBI の `StringToDouble`（既定の flag は `strtod` だけ）が読み切らない値は NCBI の誤り（`Invalid input for filtering parameters`）。`strtod` が読み、LOSAT が読まない値（16 進数、無限大、1e400）だけが明示的な拒否 |
| BP-2 | 高 | BLASTP の `-max_target_seqs` 1073741799〜2147483623：NCBI は落ちる、LOSAT は実行 | 明示的な拒否（予備の hit list の大きさを NCBI の `int` で計算し、正でなければ拒否。TBLASTN と同じ） |
| BP-3 | 高 | 同じく 2147483624 以上：NCBI の予備の大きさは 2〜48 に回り込む | 直した（`blastn/hsp.rs` の `get_prelim_hitlist_size` を共有）。fixture `e2e.blastp.mts_wrap` |
| BP-4・TN-6 | 高・中 | 無限大の `-evalue`（`+inf`、`1e999`）：NCBI の blastp・tblastn は入力によって SIGSEGV（300 の subject、compressed の lookup）、LOSAT は結果を出す | 明示的な拒否（判断 D12、BLASTP と TBLASTN。有限の 1e308 は実行し NCBI と一致）。既存の fixture `e2e.blastp.evalue_huge`・`huge7`・`e2e.tblastn.evalue_inf` を `-evalue 1e308` にして固定し直した |
| BP-5 | 高 | `-outfmt "6 sframe"`・`"6 frames"`：NCBI は他に整列を求める field が無いと 0 と 0/0 | 直した（tabular.cpp:900-908 の条件、`blastp_tabular_sets_frames`）。7 つの field の組で一致 |
| BP-8・TN-7 | 低・中 | `-dryrun`（NCBI の隠れた toolkit の option）が LOSAT の不明な option の誤り | 直した：D9 のとおり明示的な拒否（4 つの program） |
| TN-1 | 高 | TBLASTN の `-comp_based_stats 0` で最終の X-drop が小さい（または `INT_MIN` になる巨大な値）と一部の e-value が違う | 直した：NCBI は部分翻訳の左端に、先頭からの翻訳でも fence を置き（blast_hits.c:1224-1228）、traceback が触れると全翻訳でやり直して `stat_length` が全長になる。LOSAT は先頭からの部分翻訳に fence を置かなかった。fixture `e2e.tblastn.fence_xdrop`・`fence_xdrop_huge` |
| TN-3 | 高 | outfmt 0 の gap だけの Sbjct の行に座標 | 直した（showalign.cpp:1605-1627、BLASTX の書き方と同じ） |
| IN-13・BP-12・TX-7・RP-2・RP-3・TN-10・TN-11 | 低 | NCBI の参照の行のずれ、断片の無い注釈 | 直した（このセッションで足したもの 7 件と前からの 4 件） |
| TX-5・BP-7・TN-9 | 中〜低 | LOSAT の明示的な拒否（対応しない outfmt、入力の拒否、スレッドの上限）が NCBI の option の誤りより先に出る | 受け入れ：両方とも終了コード 1 の明示的な失敗（SD の監査 (c) の F5 と同じ扱い） |
| TX-6・BP-9・TN-8 | 低 | スレッドの上限と `NCBI_CONFIG__*` の拒否の文言に program の名が無い（「not supported by LOSAT」） | 受け入れ：全 program と web ABI v1 が共有する文言（TD-1） |
| TX-8・IN-12 | 低 | web ABI v1 の引数の読み方（`inf`・`nan`、0 の個数） | 受け入れ：v1 は凍結（TD-1） |
| TX-9 | 低 | `validate` は入力を見ないので、空の query でも LOSAT の上限を拒否する（CLI は `Query is Empty!`） | 受け入れ：空でない query の実行は同じ拒否になる |
| TX-10 | 低 | AUTHORITY.md の §K と §M の window の記述の食い違い | 直した（§K・§M・§N） |
| BP-10 | 低 | `LOSAT_STARTUP_TRACE` の診断の出力 | 受け入れ：LOSAT の診断の環境変数（abi_v2.md の `LOSAT_TIMING` と同じ扱い） |
| BP-11 | 低 | `--help` の終了コード | 承認済みの例外 1 |
| BP-8（残り） | 低 | `-out -version` で LOSAT は `-version` という名のファイルを作る（NCBI は version を出す）。末尾の `--` | `--` は例外 1（引数の構文）。`-out -version` は S08+b |
| TN-2 | 高 | TBLASTN の `-comp_based_stats 0`、`-lcase_masking -seg no`（hard mask）、少ない `-max_target_seqs`：和の統計の e-value が 1〜4 % 違う | **S08+b**：`-sum_stats false` では一致。NCBI の linking は `Blast_HSPListGetEvalues` を subject の長さ/3 で呼ぶ（trace：`link` 489、`evalue` 163）。LOSAT の Spouge の長さ（`db_length`）か、予備の段で残る subject に依る値を確かめる |
| TN-4 | 高 | 巨大な `-xdrop_gap`・`-xdrop_gap_final`（1e7〜5e8 bit）で LOSAT は数十秒・数 GB（NCBI は 1 秒未満）。出力は同じ | **S08+b**：LOSAT の DP の帯（`blastp/gapalign.rs` の `gap_dp_reserve_*`）は `resize` で全要素を書く。NCBI は `malloc`（触れない page は確保されない）。BLASTP・TBLASTN が共有する経路で V-PERF が要る |
| TN-5 | 高 | BLOSUM45 の組で `-evalue` 5000 以上：同じ得点の HSP の frame・座標が違う（29900 行中 4〜12 行） | **S08+b**：同点の順（blast_hits.c の並べ替え・heap） |
| RP-4 | 低（前から） | 既定の option の TBLASTN：60000 残基の query の 4998 残基の一致で bit score が違う（NCBI 10404、LOSAT 10412。変更前の実行ファイルも同じ） | **S08+b**：composition の行列の調整。入力は [`round1/rp4/`](round1/rp4/)（gzip） |

TN-2・TN-5 の入力は [`round1/tblastn_inputs/`](round1/tblastn_inputs/)。

## 確かめ

直した指摘は、監査の再現の命令を直した実行ファイルで繰り返して、NCBI とのバイト一致か明示的な拒否を確かめた（stdout、stderr、終了コード、`2>&1`）。fixture 137 件（1・4 スレッド）は差 0、変更前の実行ファイル（SD の最後）は新しい fixture で違う。単体試験と CLI の試験を足した（`query_warnings.rs`、`input.rs`、`app.rs`、`value_parsers.rs`、`blast_engine.rs`、`tests/cli_v2.rs`）。
