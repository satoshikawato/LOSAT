# 次のセッションへのプロンプト

## 概要
LOSATのTBLASTXでSEGフィルターを実装し、runaway extensionの問題は解決した（最大アラインメント長: 113,828 aa → 2,253 aa）。しかし、**ヒット数がNCBI BLAST+の約1/6**（9,789 vs 62,053）である問題が残っている。

## 現状
- SEGフィルター実装済み（window=12, locut=2.2, hicut=2.5）
- マスク領域のアミノ酸を'X'に置換して拡張を自然に停止
- グラフから見ると、low identity（<80%）と短いヒット（<30 aa）が特に不足

## 次のタスク

### 調査1: SEGマスク率の比較
LOSATは41.33%のアミノ酸をマスクしているが、NCBI BLASTの実際のマスク率と比較する。

```bash
# NCBI segmaskerでマスク率を確認
segmasker -in fasta/NZ_CP006932.fasta -outfmt interval
```

### 調査2: シード検出の比較
`utils.rs:502`で`seed_score < 11`のシードをスキップしている。この閾値を緩和または削除してヒット数を比較。

### 調査3: 診断出力の追加
シード検出数、フィルタリング数、拡張数の詳細統計を出力して、どこでヒットが失われているか特定。

## 参照ファイル
- まとめ: `docs/next_session_tblastx_v10.md`
- SEG実装: `LOSAT/src/utils/seg.rs`
- TBLASTX統合: `LOSAT/src/algorithm/tblastx/utils.rs`
- テスト結果: `tests/losat_out/`, `tests/blast_out/`
- NCBI BLAST: `C:\Users\kawato\Documents\GitHub\ncbi-blast`

## コマンド例
```bash
# LOSATテスト実行
cd LOSAT
time ./target/release/losat tblastx \
  -q ../tests/data/NZ_CP006932.fasta \
  -s ../tests/data/NZ_CP006932.fasta \
  -o tests/losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.out \
  2>&1 | tee tests/losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.log

# シード閾値を無効にしてテスト（utils.rsを編集後）
# seed_score < 11 の行をコメントアウト
```


