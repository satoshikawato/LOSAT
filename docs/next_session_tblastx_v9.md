# TBLASTX セッション v9: SEGフィルター実装

## セッション目標

NCBI BLAST互換のSEGフィルターを実装し、runaway extension問題を解決する。

## 現状の問題

### Runaway Extension
- **NCBI BLAST**: 最大アラインメント長 2,260 aa (NZ_CP006932 self-comparison, genetic code 4)
- **LOSAT**: 最大アラインメント長 113,828 aa（約50倍！）

### センチネル実装の限界
- センチネル実装は完了したが、runaway extension問題は解決していない
- 配列末端のセンチネルは意味がない（拡張は自然にそこで止まる）
- Self-comparisonでは100%一致が続くため、センチネルペナルティだけではX-dropが発動しない

## NCBI BLASTのSEGフィルター

### 基本情報
- **オプション**: `-seg <String>`
- **デフォルト**: `12 2.2 2.5` (window=12, locut=2.2, hicut=2.5)
- **用途**: 低複雑度領域を検出してマスク
- **効果**: マスクされた領域で拡張が停止し、runaway extensionを防ぐ

### 実装場所
- NCBI BLASTソースコード: `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast`
- SEGアルゴリズムの実装を調査してLOSATに移植

### 参考資料
- SEG論文: Wootton & Federhen (1993) "Statistics of local complexity in amino acid sequences and sequence databases"
- NCBI BLASTのSEG実装を調査

## 実装タスク

### 1. SEGアルゴリズムの調査
- [ ] NCBI BLASTのSEG実装を特定
- [ ] SEGアルゴリズムの詳細を理解
- [ ] パラメータ（window, locut, hicut）の意味を確認

### 2. SEGフィルター実装
- [ ] `src/utils/seg.rs` を作成（または既存の低複雑度フィルターを拡張）
- [ ] SEGアルゴリズムを実装
- [ ] アミノ酸シーケンスに適用
- [ ] マスクされた領域を `MaskedInterval` として返す

### 3. TBLASTXへの統合
- [ ] `translation.rs` でSEGフィルターを適用
- [ ] マスクされた領域を `query_masks` に追加
- [ ] 拡張時にマスク領域をチェックして停止

### 4. テスト・検証
- [ ] 単体テストを作成
- [ ] NZ_CP006932 self-comparisonで検証
- [ ] NCBI BLASTと比較（最大アラインメント長が2,260 aa程度になることを確認）

## 参考データ

### テストケース
- **Query/Subject**: `LOSAT/tests/fasta/NZ_CP006932.fasta`
- **Genetic Code**: 4
- **NCBI結果**: `LOSAT/tests/blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out`
  - 最大アラインメント長: 2,260 aa
  - ヒット数: 62,053
- **LOSAT結果（現状）**: `LOSAT/tests/losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.out`
  - 最大アラインメント長: 113,828 aa
  - ヒット数: 50,952

### 解析ツール
- `LOSAT/tests/analyze_tblastx_feature_overlap.py` - GFF3解析スクリプト

## 実装の注意点

1. **既存のDustMaskerとの関係**
   - LOSATには既に `src/utils/dust.rs` がある（DNA用）
   - SEGはアミノ酸用の低複雑度フィルター
   - 両方を適切に使い分ける

2. **マスクの適用タイミング**
   - 翻訳後、アミノ酸シーケンスに対してSEGを適用
   - マスク情報をDNA座標に変換して `query_masks` に追加

3. **パラメータの互換性**
   - NCBI BLASTのデフォルト `12 2.2 2.5` を使用
   - 将来的にはコマンドラインオプションで変更可能にする

## 成功基準

- [ ] NZ_CP006932 self-comparisonで最大アラインメント長が2,500 aa以下になる
- [ ] NCBI BLASTとの最大アラインメント長の差が10%以内
- [ ] 単体テストがすべてパス
- [ ] 既存のテストが壊れない

## 次のセッションへの引き継ぎ

このセッションで完了できなかった場合:
1. SEG実装の進捗状況
2. 残っている問題点
3. 次の優先順位

