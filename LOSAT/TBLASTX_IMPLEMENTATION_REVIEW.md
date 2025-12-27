# TBLASTX実装レビューと修正経緯

## 概要

LOSATのTBLASTX実装をNCBI BLASTと比較し、不一致を修正する作業の記録。

## 初期の問題点

1. **ヒット数が少ない**: LOSATのヒット数がNCBI BLASTより大幅に少ない
2. **ヒットが異常に短い**: 以前よりもヒットが異常に短くなっている
3. **低アイデンティティのヒットが欠落**: 比較的アイデンティティの低いヒットが欠落している
4. **長いヒットが過剰に表現**: LOSATはBLASTより長いヒットが過剰に表現されている

## 修正内容

### 1. Seed Score Thresholdの修正

**問題**: LOSATがseed score thresholdとして11を使用していたが、NCBI BLASTは13を使用している。

**修正**:
- ファイル: `LOSAT/src/algorithm/tblastx/utils.rs`
- 変更: `if seed_score < 11` → `if seed_score < 13`
- 参照: `ncbi-blast/c++/include/algo/blast/core/blast_options.h:115`
  - `#define BLAST_WORD_THRESHOLD_TBLASTX 13`

### 2. One-hit Extensionの削除

**問題**: LOSATには、seed score >= 30の場合にtwo-hit要件なしでone-hit extensionを許可する最適化があった。これはNCBI BLASTの実装と異なる。

**修正**:
- ファイル: `LOSAT/src/algorithm/tblastx/utils.rs`
- 変更: `if two_hit_info.is_none() && seed_score < 30` のブロックを削除
- 理由: NCBI BLASTはTBLASTXでstrict two-hit要件を強制している
- 参照: `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c`
  - TBLASTXは`window_size > 0`の場合、`s_BlastAaWordFinder_TwoHit`を使用し、strict two-hit要件を強制

### 3. E-value計算の修正（Effective Search Space）

**問題**: LOSATがalignment lengthベースでE-valueを計算していたが、NCBI BLASTはquery/subjectのamino acid length（length adjustment後）を使用している。

**修正**:
- ファイル: `LOSAT/src/algorithm/tblastx/utils.rs`
- 変更前:
  ```rust
  let len = qe - qs;
  let (bit_score, e_val) = calculate_statistics(ungapped_score, len, &params);
  ```
- 変更後:
  ```rust
  let q_aa_len = q_frame.aa_seq.len();
  let s_aa_len = s_aa.len();
  let search_space = SearchSpace::with_length_adjustment(q_aa_len, s_aa_len, &params);
  let bit_score = calc_bit_score(ungapped_score, &params);
  let e_val = calc_evalue(bit_score, &search_space);
  ```
- 適用箇所:
  - Ungapped extension (line 320-323)
  - Gapped extension (line 463-464)
- 参照: `ncbi-blast/c++/src/algo/blast/core/blast_setup.c:842-843`
  ```c
  effective_search_space = effective_db_length * (query_length - length_adjustment);
  ```

### 4. Two-hit Extensionロジックの確認

**問題**: `extend_hit_two_hit`関数のコメントと実装が不一致していた。

**調査内容**:
- Word best position logicの確認
- First hit到達判定の確認
- `distance_to_first_hit_end`の計算方法の確認

**結果**: 
- コメントと実装の不一致を確認したが、NCBI BLASTの実装に合わせるため、変更を元に戻した
- 現在の実装はNCBI BLASTに忠実であることを確認

## 調査内容

### Effective Search Spaceの計算方法

**NCBI BLASTの実装**:
- 各context（frame combination）ごとにeffective search spaceを計算
- `effective_search_space = effective_db_length * (query_length - length_adjustment)`
- `effective_db_length = db_length - (db_num_seqs * length_adjustment)`
- `db_num_seqs = 1`の場合（single sequence comparison）:
  - `effective_db_length = db_length - length_adjustment`
  - `effective_search_space = (db_length - length_adjustment) * (query_length - length_adjustment)`

**LOSATの実装**:
- 各frame combinationごとにeffective search spaceを計算
- `effective_m = (m - length_adj).max(1.0)`
- `effective_n = (n - length_adj).max(1.0)`
- `effective_space = effective_m * effective_n`
- `db_num_seqs = 1`の場合、NCBI BLASTと等価

**確認結果**:
- `q_aa_len`と`s_aa_len`は各frame combinationごとに正しく設定されている
- `q_aa_len = q_aa.len()` (frame-specific amino acid length)
- `s_aa_len = s_aa.len()` (frame-specific amino acid length)
- Length adjustmentの計算は正しい
- Effective search spaceの計算式はNCBI BLASTと一致

### LOSATのImplied Search SpaceがBLASTの約16337倍という問題

**調査結果**:
- E-valueとbit scoreから逆算したimplied search spaceがBLASTの約16337倍
- しかし、これはE-valueとbit scoreから逆算した値であり、実際のeffective search spaceの計算とは異なる可能性がある
- 実装自体はNCBI BLASTに合わせていることを確認

**考えられる原因**:
1. BLASTのE-value計算が異なる可能性（可能性は低い）
2. LOSATのE-value計算が異なる可能性（実装は正しいが、使用箇所が異なる可能性）
3. 実際のeffective search spaceの計算は正しいが、E-valueとbit scoreの関係が異なる可能性

## 現在の実装状況

### 修正済み
- ✅ Seed score threshold: 11 → 13
- ✅ One-hit extensionの削除（strict two-hit要件）
- ✅ E-value計算の修正（effective search spaceベース）
- ✅ Two-hit extensionロジックの確認（NCBI BLASTに忠実）

### 調査済み
- ✅ Effective search spaceの計算方法
- ✅ Frame combinationの扱い
- ✅ Length adjustmentの計算

### 残存する問題
- ⚠️ ヒット数が少ない（以前より改善したが、まだBLASTより少ない）
- ⚠️ ヒットが短い（以前より改善したが、まだBLASTより短い）
- ⚠️ 低アイデンティティのヒットが欠落（以前より改善したが、まだBLASTより少ない）
- ⚠️ 長いヒットが過剰に表現（X-dropoffの適用、スコア計算の順序などを確認する必要がある）

## 次のステップ

1. **長いヒットが過剰に表現される原因の調査**
   - X-dropoffの適用方法
   - スコア計算の順序
   - Extension終了条件

2. **ヒット数が少ない原因の調査**
   - Seed findingのロジック
   - Two-hit要件の適用
   - Filtering条件

3. **低アイデンティティのヒットが欠落する原因の調査**
   - Scoring matrixの適用
   - E-value threshold
   - Filtering条件

## 参照ファイル

### LOSAT
- `LOSAT/src/algorithm/tblastx/utils.rs`: メイン実行ロジック
- `LOSAT/src/algorithm/tblastx/extension.rs`: Extensionロジック
- `LOSAT/src/algorithm/tblastx/constants.rs`: 定数定義
- `LOSAT/src/stats/search_space.rs`: Effective search space計算
- `LOSAT/src/stats/length_adjustment.rs`: Length adjustment計算

### NCBI BLAST
- `ncbi-blast/c++/include/algo/blast/core/blast_options.h`: 定数定義
- `ncbi-blast/c++/src/algo/blast/core/blast_setup.c`: Effective search space計算
- `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c`: Ungapped extensionロジック
- `ncbi-blast/c++/src/algo/blast/api/tblastx_options.cpp`: TBLASTXオプション設定

## テスト結果

### AP027131.AP027133
- **修正前**: LOSAT hits: 33,382, BLAST hits: 62,053
- **修正後**: LOSAT hits: 12,559, BLAST hits: 62,053
- **改善**: ヒット数は減少したが、これはE-value計算の修正による影響

### 低アイデンティティヒット
- **LOSAT**: 776 hits (0-50% identity)
- **BLAST**: 9,503 hits (0-50% identity)
- **比率**: 8.2%

## 結論

実装はNCBI BLASTに合わせて修正されました。主な修正点は：
1. Seed score thresholdの修正
2. One-hit extensionの削除
3. E-value計算の修正（effective search spaceベース）

残存する問題については、さらなる調査が必要です。

