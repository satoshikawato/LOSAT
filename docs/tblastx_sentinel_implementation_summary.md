# TBLASTX センチネル実装セッション - まとめ

## 実施日
2025年12月30日

## 目標
NCBI BLAST互換のセンチネル実装により、runaway extension問題を解決する。

## 実施内容

### 1. NCBI BLASTのセンチネル方式の調査

**発見事項:**
- NCBI BLASTは翻訳シーケンスの先頭と末尾にNULLB (0) センチネルを配置
- センチネルは `blast_util.c:BLAST_GetTranslation()` で追加される
- センチネル位置: `prot_seq[0] = NULLB`, `prot_seq[index_prot] = NULLB`
- センチネルは `IS_residue(x)` で `x <= 250` として判定される
- BLOSUM62マトリックスでセンチネル（未知残基）のスコアは `defscore = -4`

**参考コード:**
- `ncbi-blast/c++/src/algo/blast/core/blast_util.c:438-453`
- `ncbi-blast/c++/include/algo/blast/core/blast_util.h:48` (`IS_residue`)

### 2. LOSATへのセンチネル実装

**実装した変更:**

#### 2.1 `translation.rs`
- `translate_sequence()` 関数でシーケンス先頭・末尾に `SENTINEL_BYTE (255)` を追加
- `QueryFrame` 構造体に `aa_len` フィールドを追加（センチネル除く実際のアミノ酸数）
- レイアウト: `[SENTINEL_BYTE, aa_0, aa_1, ..., aa_n-1, SENTINEL_BYTE]`

#### 2.2 `lookup.rs`
- k-merスキャン時にセンチネル位置（0とlen-1）をスキップ
- 論理位置（raw位置-1）をルックアップテーブルに保存
- `encode_aa_kmer()` でセンチネル（255）を検出してNoneを返すように修正

#### 2.3 `extension.rs`
- `get_score()` 関数でセンチネル（255）を検出したら `SENTINEL_PENALTY (-100)` を返す
- これによりX-dropが即座に発動し、拡張が終了

#### 2.4 `utils.rs`
- 座標変換: RAW座標（センチネル込み）とLOGICAL座標（センチネル除く）を適切に使い分け
- 拡張関数にはRAW座標を渡し、結果をLOGICAL座標に変換
- `aa_len` を使用して検索空間計算などを行う

#### 2.5 `constants.rs`
- `SENTINEL_BYTE = 255` を定義
- `SENTINEL_PENALTY = -100` を定義（X-drop=15より十分に大きい）

### 3. テスト結果

**ビルド・単体テスト:**
- ✅ ビルド成功
- ✅ TBLASTX関連の単体テスト67件すべてパス

**実データテスト (NZ_CP006932, genetic code 4):**
- ❌ **問題発見**: LOSATは最大113,828 aaまで拡張（NCBIは2,260 aa）
- センチネル実装は完了したが、runaway extension問題は解決していない

### 4. GFF3解析による原因調査

**実施した解析:**
- NCBI TBLASTX出力とGFF3注釈の重なりを解析
- LOSAT TBLASTX出力とGFF3注釈の重なりを解析

**発見事項:**
- NCBIの2,260 aaヒットは**ORF境界で止まっているわけではない**
- 複数のCDSをまたいでいる（例: rpoC, rpoB, rplL, rplJ）
- LOSATの113,828 aaヒットは**315個のCDSをまたいでいる**
- センチネルだけではrunaway extensionを防げない

**統計:**
- NCBI: 最大アラインメント長 2,260 aa
- LOSAT: 最大アラインメント長 113,828 aa（約50倍）
- NCBI: 62,053 hits
- LOSAT: 50,952 hits

### 5. 問題の本質

**センチネル実装の問題点:**
1. 配列末端のセンチネルは意味がない（拡張は自然にそこで止まる）
2. センチネルは境界到達時にペナルティを適用するが、self-comparisonでは100%一致が続くためX-dropが発動しない
3. NCBIが2,260 aaで止まる理由は**SEGフィルター**による可能性が高い

**NCBI BLASTのSEGフィルター:**
- `-seg` オプションで低複雑度領域をマスク
- デフォルト: `12 2.2 2.5` (window=12, locut=2.2, hicut=2.5)
- マスクされた領域では拡張が止まる

## 未解決の問題

1. **Runaway extension**: LOSATは113,828 aaまで拡張する（NCBIは2,260 aa）
2. **センチネルの効果不足**: センチネル実装だけでは問題解決しない
3. **SEGフィルター未実装**: LOSATにはSEGフィルターがない

## 次のステップ

### 優先度1: SEGフィルター実装
- NCBI BLASTの `-seg` オプション相当を実装
- 低複雑度領域を検出してマスク
- マスクされた領域で拡張を停止

### 優先度2: 拡張終了条件の再検討
- X-drop以外の終了条件を検討
- NCBIの実際の動作をさらに詳細に調査

### 優先度3: ヒット数の差の原因特定
- NCBI: 62,053 hits
- LOSAT: 50,952 hits
- 約18%の差の原因を特定

## 参考資料

- NCBI BLASTソースコード: `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast`
- 解析スクリプト: `LOSAT/tests/analyze_tblastx_feature_overlap.py`
- テストデータ: `LOSAT/tests/fasta/NZ_CP006932.fasta`
- GFF3注釈: `/mnt/d/study/2025-05-04_gbdraw_test/NZ_CP006932.gff3`

## 実装ファイル一覧

- `LOSAT/src/algorithm/tblastx/translation.rs` - センチネル追加
- `LOSAT/src/algorithm/tblastx/lookup.rs` - センチネル位置スキップ
- `LOSAT/src/algorithm/tblastx/extension.rs` - センチネルペナルティ
- `LOSAT/src/algorithm/tblastx/utils.rs` - 座標変換
- `LOSAT/src/algorithm/tblastx/constants.rs` - センチネル定数
- `LOSAT/tests/analyze_tblastx_feature_overlap.py` - GFF3解析スクリプト


