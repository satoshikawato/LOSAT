# E2h（SF・SFb）ゲート記録：FASTA の読み方（NCBI `CFastaReader` の移植）

状態：**途中**（2026-10-08、SF を区切った）。完了条件（計画 §7 の SF の行）はまだ満たしていない。続きは [SFb の指示書](../../losat_web_gui_sessions/session_sfb_e2h_port.md)。

## SF（2026-10-06 開始、10-06 に `/mnt/c` の障害で中断、10-08 に Linux の clone で再開）で行ったこと

| 手順 | 結果 | コミット |
|---|---|---|
| 0. 準備 | `origin/main`（#113〜#115）を merge（衝突なし）。変更前の成果物（S11 の 8 つ）は `~/.cache/losat-web-gui-target/sf/bin/`、capture の基準は `sf/capture-before/hashes.tsv`（236 件、10-06 に取得） | `00fb6c0e` |
| 取り残された作業の復元 | 10-06 に未コミットだった読み込み器の module と棚卸しの brief を、セッション監査の写し（`stranded-backup-20261008`、SHA-256 を照合）から戻した | `eba45f75` |
| 読み込み器の単体試験 | 初めて流した 9 件のうち 4 件が失敗していた。4 件とも試験の期待が NCBI と違い、読み込み器は NCBI と同じだった。期待を NCBI のソースと oracle の観察に合わせた | `04cc336d` |
| 2. 棚卸し | [`INVENTORY.tsv`](INVENTORY.tsv)（222 行）。LR 行の読み込みと流れ 30、RD レコード 49、BI blastinput の層 55、RP 報告 56、AD アダプタと ABI 32。`status_before`（`f3048ffde`）：rejected-now 74、faithful 52、divergent 46、n/a 26、unported 22、exception 2。範囲ごとの指示と覚え書きは `inventory/`、oracle の入力と出力は `$BUILD_ROOT/sf-e2h/inventory/scratch_*`（追跡しない） | `e1eed966` |
| 1. 権威の記録 | [`AUTHORITY.md`](AUTHORITY.md)（§A〜§L）。§K に NCBI の不具合に見える挙動 12 件（保守者への一括の問い、下記） | `2a903cd7` |
| 3. 移植の計画 | [`PORT_PLAN.md`](PORT_PLAN.md)：段階 S0〜S10、記録の型、危険、未決の問い（判断は下記） | この記録と同じコミット |
| 3. 移植 S0 | 読み込み器の下準備（program への接続は無く、出力は変わらない）。Seq-id の判定を `CSeq_id::Set` の解析に合わせた（`0123`・`ACGT ACGT`・`contig1` は FASTA として読む。広めの拒否は accession の guide の形だけ）。行の読み方の試験（BLASTX の C++ oracle 12180 行、LR の癖）、`from_bytes` と file の一致の性質試験、まとめて読む経路（R1） | `24f56b17`、`905597ba` |
| 5. fixture（NCBI 側） | `LOSAT/tests/fasta_input_fixtures.py` と program ごとの manifest（`fasta_input_fixtures/*.tsv`、528 KB）、入力 683（387 KB）。3932 行：NCBI 2.17.0 を固定した 3908 行と、Seq-id の拒否の 24 行（文言は移植で埋める）。3 回の凍結でハッシュが一致 | `352ad49f` |

### 検証

- `cargo fmt --check`、release build、clippy `--all-targets --all-features -D warnings`、`cargo test --all-features`（`fasta_reader`：20 通過、速度の試験 1 は ignored）。ログは `/home/kawato/losat-baselines/sf-e2h-20261008/logs/`。
- fixture の `check` を変更前の実行ファイル（S11）で流した：same 1131、differs 13（BLASTN の `>` だけの subject：NCBI は `Subject sequence contains no data` の警告と終了コード 3 の `The average subject length is too short`、LOSAT は `Empty CBlastQueryVector`。移植の S4 で直す）、LOSAT の拒否 2764、pending 24。fixture は変更前の拒否と移植を区別する。
- 読み込み器だけの速度（release、3 回の中央値、`bio::io::fasta` との比）：100 Mb の genome（80 桁）1.6 倍、1 行の genome 1.2〜1.3 倍、10 万 × 300 塩基 1.5 倍（まとめて読む経路の前は 4.6〜6.0 倍）。V-PERF は program に接続した後。

## 推奨の案で進めた判断（保守者に委ねられた判断、2026-09-29 の常設の指示、10-07 に再掲）

1. 移植の順は `PORT_PLAN.md` の S0〜S10。`AUTHORITY.md` を S1 の前に書いた。
2. **ABI v1**：v1 は `bio` と今の検査のまま、`FastaRecord::from_bio` で `run_local` に入る（v1 の出力と文言は変わらない）。保守者の判断 4 の目的（TD-1 の凍結）を保つ。判断 4 の前提「v1 が受け付ける入力は bio と NCBI で読み方が同じ」は TBLASTX と BLASTP の v1 で成り立たない（棚卸し AD-24・AD-25）。
3. **Seq-id の行**：`CSeq_id` の解析（BI-15〜19・55）を移した。広めに拒否するのは accession の guide の形（文字と数字の数の 26 の形、先頭の byte が `A-Z`・`_`・`?` のときだけ）。guide の表（1458 の規則）は移さない。`ZZ123456` などは NCBI では FASTA だが LOSAT は拒否する（判断 2 の範囲）。
4. Web の `register` は data loader を有効とし（CLI の既定）、読み込みの誤りと注釈だけの query を NCBI の文言で早く拒否する。
5. 空の `DATA_LOADERS=` は項目なし（data loader は有効、`ncbireg.cpp:984-991`）。`blastdb` も `genbank` も含まない空でない値は両方を無効にする。
6. 読み込み器の API：Seq-id の拒否の後は、拒否した 1 行の次から読める（局所 ID の番号は進まない）。program は拒否で止まる。
7. `from_bytes` は同じ byte の通常の file と同じに読む（NCBI の stream の補充の大きさによる尾の消失を含む）。
8. まとめて読む経路は `FastaStream::bulk` の旗の後ろに置く（program では常に有効。試験では 1 byte ずつの参照の経路と比べる）。
9. fixture の manifest は情報を落とさず小さくした（全長の SHA-256、行の定義は script の生成器、program ごとに 4 file）。CI の速い検査への接続は移植の後（SFb）。TBLASTX と BLASTP の `-lcase_masking` の 34 行は残し、移植か明示的な拒否かを SFb で決める。
10. `PORT_PLAN.md` の Q4〜Q6・Q8・Q10・Q11 は推奨のとおり。BI の覚え書き §10 の推奨（query は batch ごとに読む、subject にも同じ Seq-id の規則、空のレコードと空の batch を 4 program で移植、空の pipe の零 batch、`-parse_deflines` は拒否のまま）も同じ。

## 保守者への一括の問い（`PD-LOSAT-NCBI-DEFECTS`）

`AUTHORITY.md` §K の 12 件。移植は各行の推奨の規則で進める（再現：2〜7・9・11 の切り詰め、拒否：1・8・11 の巨大な gap・12、Linux の意味：10）。保守者の答えで変わるのは主に次の 3 つ：(1) RP-20 tabular の `Subject_` の語と復号できない題で NCBI が終了コード 255 で落ちる → 明示的な拒否、(4) 改行の無い最後の行の尾が届き方（file と pipe）で失われたり残ったりする → file の挙動を再現、(3) 行の途中の単独の CR が行の終わりを失わせる → 再現（BLASTX の移植と同じ）。

## 残件（SFb の最初の作業）

- 移植 S1〜S10（`PORT_PLAN.md`）、fixture の Seq-id の 24 行と CI の速い検査への接続、`-lcase_masking` の判断。
- 指示書の手順 6〜9（sweep、ゲートの全工程、V-PERF、独立監査）。
- アプリ側（`feature/losat-web-gui-app`、S13 の途中、`e200e60b`）の merge は S13 の終了後。
