# SESSION_ID — セッション実行記録

このファイルは未記入のテンプレートである。`SESSION_ID`はS00〜S12の実IDに置換する。`WORK_DIR/reports/SESSION_ID.md`へ保存し、同じsessionの再試行は同ファイルに履歴を追記する。既存repositoryのevidence形式がある場合は項目を対応付け、同じrawデータを複製しない。対象外項目は`N/A: 理由`、未実施は`NOT_RUN: 理由`とし、空欄のまま合格を付けない。

## 1. 結論

実施日・担当・役割: 未記入  
状態: NOT_STARTED  
今回の判断（採用/不採用/適用なし/不足資料）: 未記入  
この判断を支持する最小の証拠: 未記入  
独立レビュー状態: 未実施

## 2. 作業identityと保護した範囲

| 項目 | 値 |
|---|---|
| PLAN_DIR / REPO_ROOT / WORK_DIR | 未記入 |
| 初期固定base | 未記入 |
| 今回の直近採用base | 未記入 |
| candidate SHA / 未commit patch hash | 未記入 |
| base / candidate artifact hash | 未記入 |
| branch・開始時dirty状態 | 未記入 |
| 関連gbdraw SHA / asset hash | 未記入またはN/A |
| toolchain・Cargo.lock・target・features・実flags | 未記入 |
| runtime / browser / OS / CPU / 電源・affinity | 未記入 |
| 他者の変更と今回触っていない範囲 | 未記入 |

## 3. 一つの仮説・適用経路・実装

対象program/task/case、source file/function、実call path: 未記入

計測で確認した負担と、単なる推測を分ける: 未記入

不変条件（入力/encoding、状態遷移、採用順序、出力/error）: 未記入

NCBIの対応source版・file・関数・行範囲・caller: 未記入

今回の変更責任、state所有者と寿命、消すcopy/初期化/不要計算: 未記入

検討した最小の代替と、今回選ばなかった理由: 未記入

変更fileとdiff/patchの実在path: 未記入

## 4. 互換性・安全性の証拠

| ケース / target / thread | 契約・golden所在/hash | 実command・raw出力所在/hash | 結果・最初の差 |
|---|---|---|---|
| 未記入 | 未記入 | 未記入 | NOT_RUN |

Native Gate A/Bを適用した範囲とfingerprint: 未記入またはN/A

scalar/optimized、serial/parallel、error/反復run、境界テスト: 未記入

安全性検査の実施範囲と対象外（SIMD、共有memory、borrow、overflow等）: 未記入

未確認の差異を隠すnormalization/比較器変更がないこと: 未確認

## 5. 性能と費用

S00の測定契約ID・fixture一覧・今回の選択理由: 未記入

計測範囲、cold/warm、準備再利用条件、output sink、順序/warmup/反復数: 未記入

| ケース / mode | base中央値 | candidate中央値 | 差（ms・%） | sample数・min/max | memoryと補助指標 | 判定 |
|---|---:|---:|---:|---|---|---|
| 未記入 | — | — | — | NOT_RUN | NOT_RUN | NOT_RUN |

同一threaded artifactのn1/n2/n4/n8、実spawnと仕事分散: 未記入またはN/A

copy bytes / fill bytes / scratch生成 / redo計算 / early不使用 / heap判定: 未記入またはN/A

初期base比・直近base比・Wasm/Native比の区別: 未記入

whole searchとmicrobenchmarkの区別、前計算を含む収支: 未記入

測定ノイズ、CPU/worker時間の重複、失敗run、未測定platform: 未記入

## 6. 採否と原則レビュー

事前に固定した採用/非回帰/memory基準への適合: 未記入

SRP: preparation/compute/reduction/hostの責任が混ざっていないか: 未記入

OCP/LSP: 既存契約の置換であり別の結果を許していないか: 未記入

ISP/DIP: 必要な入力だけを渡し、kernelがJS/WASI/計測へ依存していないか: 未記入

KISS/DRY/YAGNI: 必要な局所変更か、重複framework/設定/文書を増やしていないか: 未記入

棄却時の差分撤去・保持すべき証拠、採用時の統合条件: 未記入

## 7. 再実行・引き継ぎ・rollback

実際に実行したbuild/test/benchmark commandを、そのまま再実行できる形で記す。S00の共通commandを使った場合は識別子と差分引数、CWD、環境、実artifact pathを記す。秘密情報を含めない。

```text
未記入: このblockはcommandではない。実行時に実在commandへ置換する。
```

生証拠の実在path、入力hash、report間の依存: 未記入

次sessionのID・着手条件・読むfile/関数・最初のcommand: 未記入

BLOCKED/INCONCLUSIVEの場合に必要な具体的資料・環境: 未記入またはN/A

rollbackする採用commit/patchの範囲と方法: 未記入またはN/A

STATUS更新内容: 未記入

## 8. 追記履歴

再試行した場合、日時、base/candidateの変更、追加sample、判断の変更理由を追記する。過去の不都合な結果を消して置き換えない。
