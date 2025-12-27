# 次のセッション用プロンプト

## 背景

LOSATのTBLASTX実装をNCBI BLASTと比較し、不一致を修正する作業を継続中です。

## 現在の状況

### 実装済みの修正

1. ✅ **Two-hit要件を満たさない場合でもextensionを実行**
   - 以前: Two-hit要件を満たさない場合は`continue`でスキップ
   - 修正後: Two-hit要件を満たさない場合でも`extend_hit_two_hit`を呼び出し
   - `s_left_off`を現在のseed位置に設定することで、left extensionのみが実行される
   - 診断情報: `ungapped_one_hit_extensions`と`ungapped_two_hit_extensions`を追加

2. ✅ **cutoff_scoreとgap_triggerの実装**
   - `CUTOFF_E_TBLASTX = 1e-300`からcutoff_scoreを計算
   - `gap_trigger`（`BLAST_GAP_TRIGGER_PROT = 22.0` bit scoreから計算）で制限
   - `cutoff_score = MIN(cutoff_from_evalue, gap_trigger)`

3. ✅ **diag_maskとdiag_coordの実装**
   - `diag_mask = diag_array_length - 1`（2のべき乗 - 1）
   - `diag_coord = (query_offset - subject_offset) & diag_mask`
   - `mask_key = (q_idx, q_f_idx, diag_coord)`

### テスト結果（修正後）

| テストケース | LOSAT hits | NCBI hits | 比率 | 改善 |
|------------|-----------|----------|------|------|
| NZ_CP006932 self | 32,865 | 62,053 | 53.0% | 2.4倍 |
| AP027132 vs NZ_CP006932 | 32,045 | 62,580 | 51.2% | 2.3倍 |
| AP027078 vs AP027131 | 11,563 | 30,180 | 38.3% | 4.4倍 |
| AP027131 vs AP027133 | 11,967 | 14,876 | 80.5% | 3.9倍 |
| AP027133 vs AP027132 | 34,334 | 54,869 | 62.6% | 2.3倍 |

### 診断情報（AP027132 vs NZ_CP006932）

- K-mer matches: 59,947,933
- Seeds filtered (low score): 17,055,652 (28.5%)
- Seeds passed to extension: 40,760,053
  - One-hit extensions: 40,691,988 (99.8%)
  - Two-hit extensions: 68,065 (0.2%)
- Total ungapped-only: 40,760,053
- Filtered (cutoff_score): 37,410,306 (91.8%)
- E-value passed: 32,045
- E-value failed: 3,317,702 (9.2%)

### 残存する問題

1. **ヒット数がまだ少ない**: 38-80%の比率（平均約53%）
2. **低アイデンティティヒットが少ない**: 21-28%の比率（NCBI BLAST: 42-63%）
3. **cutoff_scoreで91.8%がフィルタリング**: `gap_trigger=46`が制限している可能性

### Identity分布の比較

| テストケース | <50% (LOSAT/NCBI) | 50-70% (LOSAT/NCBI) | 70-90% (LOSAT/NCBI) | >=90% (LOSAT/NCBI) |
|------------|------------------|-------------------|-------------------|------------------|
| NZ_CP006932 self | 6,357 / 26,070 (24.4%) | 14,468 / 16,442 (88.0%) | 7,182 / 3,892 (184.5%) | 4,858 / 15,649 (31.0%) |
| AP027132 vs NZ_CP006932 | 5,989 / 27,325 (21.9%) | 14,101 / 17,020 (82.8%) | 7,110 / 4,515 (157.5%) | 4,845 / 13,715 (35.3%) |
| AP027078 vs AP027131 | 4,986 / 17,361 (28.7%) | 4,647 / 10,276 (45.2%) | 1,734 / 2,467 (70.3%) | 196 / 71 (276.1%) |
| AP027131 vs AP027133 | 4,748 / 9,503 (50.0%) | 5,177 / 4,516 (114.6%) | 1,830 / 747 (245.0%) | 212 / 105 (201.9%) |
| AP027133 vs AP027132 | 6,521 / 25,049 (26.0%) | 17,364 / 21,979 (79.0%) | 8,953 / 5,074 (176.4%) | 1,496 / 2,762 (54.2%) |

## NCBI BLAST実装の確認済み事項

詳細は`NCBI_BLAST_IMPLEMENTATION_DETAILS.md`を参照。

### 重要な確認事項

1. ✅ **Two-hit要件を満たさない場合**: NCBI BLASTでは`continue`でスキップされ、extensionは実行されない
   - 参照: `aa_ungapped.c:538-543`, `aa_ungapped.c:549-551`, `aa_ungapped.c:566-573`
   - **注意**: 現在のLOSATの実装は、Two-hit要件を満たさない場合でもextensionを実行している（これはNCBI BLASTと異なる）

2. ✅ **cutoff_scoreの計算**:
   - `CUTOFF_E_TBLASTX = 1e-300`から計算
   - Search space: `MIN(subj_length, query_length) * subj_length`
   - `gap_trigger`で制限: `new_cutoff = MIN(new_cutoff, gap_trigger)`
   - 参照: `blast_parameters.c:348-374`

3. ✅ **cutoff_scoreとE-valueのsearch spaceの違い**:
   - cutoff_score: `MIN(subj_length, query_length) * subj_length`（事前フィルタリング用）
   - E-value: `effective_db_length * (query_length - length_adjustment)`（最終統計評価用）
   - **結論**: 意図的な違い

4. ✅ **Extensionの実装**:
   - Left extensionは常に実行される
   - Right extensionはleft extensionがfirst hitに到達した場合のみ実行
   - 戻り値: `MAX(left_score, right_score)`
   - 参照: `aa_ungapped.c:1088-1158`

## 次の調査項目

### 優先度: 高

1. **Two-hit要件を満たさない場合の処理を再検討**
   - **問題**: 現在のLOSATは、Two-hit要件を満たさない場合でもextensionを実行している
   - **NCBI BLAST**: Two-hit要件を満たさない場合は`continue`でスキップ
   - **調査**: NCBI BLASTの実装に合わせて、Two-hit要件を満たさない場合はextensionを実行しないように修正すべきか？
   - **影響**: 現在の実装により、ヒット数が2.3-4.4倍に改善しているが、これはNCBI BLASTと異なる動作

2. **cutoff_scoreでフィルタリングされているヒットのスコア分布を分析**
   - **問題**: 91.8%のヒットが`cutoff_score`でフィルタリングされている
   - **調査**: 
     - `cutoff_score`でフィルタリングされているヒットのスコア範囲を確認
     - 低アイデンティティヒットのスコア分布を分析
     - `gap_trigger=46`が適切かどうかを判断

3. **低アイデンティティヒットが少ない原因を特定**
   - **問題**: LOSATは21-28%の比率（NCBI BLAST: 42-63%）
   - **調査**:
     - 低アイデンティティヒットのスコア分布
     - `cutoff_score`でフィルタリングされている低アイデンティティヒットの数
     - E-valueフィルタリングの影響

### 優先度: 中

4. **gap_trigger=46の適切性を確認**
   - **問題**: `gap_trigger=46`が`cutoff_score`を制限している
   - **調査**:
     - NCBI BLASTの実際の`gap_trigger`値を確認
     - 他のテストケースでの`gap_trigger`値を確認
     - `gap_trigger`が低すぎる可能性を検証

5. **cutoff_scoreの計算方法を再確認**
   - **調査**: `BLAST_Cutoffs`関数の実装を詳細に確認
   - Search spaceの計算が正しいか確認

## コードの場所

### LOSAT実装

- `LOSAT/src/algorithm/tblastx/utils.rs`: メイン実行ロジック
  - Line 361-394: Two-hit要件のチェックとextension呼び出し
  - Line 425-502: cutoff_scoreの計算とフィルタリング
  - Line 504-565: E-value計算とHSP保存

- `LOSAT/src/algorithm/tblastx/extension.rs`: Extensionロジック
  - Line 170-275: `extend_hit_two_hit`関数

- `LOSAT/src/algorithm/tblastx/constants.rs`: 定数定義
  - `CUTOFF_E_TBLASTX = 1e-300`
  - `GAP_TRIGGER_BIT_SCORE = 22.0`

- `LOSAT/src/algorithm/common/diagnostics.rs`: 診断情報
  - `ungapped_one_hit_extensions`: One-hit extensionのカウント
  - `ungapped_two_hit_extensions`: Two-hit extensionのカウント

### NCBI BLAST参照

- `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c`: Two-hit要件とextension
- `ncbi-blast/c++/src/algo/blast/core/blast_parameters.c`: cutoff_scoreとgap_trigger
- `ncbi-blast/c++/src/algo/blast/core/blast_setup.c`: Effective search space

## テストケース

### テストコマンド

```bash
export LOSAT_BIN=./LOSAT/target/release/losat

# NZ_CP006932 self
(time $LOSAT_BIN tblastx -q ./LOSAT/tests/fasta/NZ_CP006932.fasta -s ./LOSAT/tests/fasta/NZ_CP006932.fasta -o ./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1 )&>./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.log

# AP027132 vs NZ_CP006932
(time $LOSAT_BIN tblastx -q ./LOSAT/tests/fasta/AP027132.fasta -s ./LOSAT/tests/fasta/NZ_CP006932.fasta -o ./losat_out/AP027132.NZ_CP006932.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1 )&>./losat_out/AP027132.NZ_CP006932.tlosatx.n1.log

# AP027078 vs AP027131
(time $LOSAT_BIN tblastx -q ./LOSAT/tests/fasta/AP027078.fasta -s ./LOSAT/tests/fasta/AP027131.fasta -o ./losat_out/AP027078.AP027131.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1 )&>./losat_out/AP027078.AP027131.tlosatx.n1.log

# AP027131 vs AP027133
(time $LOSAT_BIN tblastx -q ./LOSAT/tests/fasta/AP027131.fasta -s ./LOSAT/tests/fasta/AP027133.fasta -o ./losat_out/AP027131.AP027133.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1 )&>./losat_out/AP027131.AP027133.tlosatx.n1.log

# AP027133 vs AP027132
(time $LOSAT_BIN tblastx -q ./LOSAT/tests/fasta/AP027133.fasta -s ./LOSAT/tests/fasta/AP027132.fasta -o ./losat_out/AP027133.AP027132.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1 )&>./losat_out/AP027133.AP027132.tlosatx.n1.log
```

### NCBI BLAST結果ファイル

- `./LOSAT/tests/blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out`
- `./LOSAT/tests/blast_out/AP027132.NZ_CP006932.tblastx.n1.out`
- `./LOSAT/tests/blast_out/AP027078.AP027131.tblastx.n1.out`
- `./LOSAT/tests/blast_out/AP027131.AP027133.tblastx.n1.out`
- `./LOSAT/tests/blast_out/AP027133.AP027132.tblastx.n1.out`

## 重要なドキュメント

- `NCBI_BLAST_IMPLEMENTATION_DETAILS.md`: NCBI BLASTの実装詳細（確認済み事項と疑問点）
- `TBLASTX_IMPLEMENTATION_REVIEW.md`: これまでの修正履歴と調査結果

## 次のセッションで最初に行うこと

1. **Two-hit要件を満たさない場合の処理を再検討**
   - NCBI BLASTの実装に合わせて、Two-hit要件を満たさない場合はextensionを実行しないように修正するか判断
   - 現在の実装（extensionを実行）とNCBI BLASTの実装（スキップ）のどちらが正しいか確認

2. **cutoff_scoreでフィルタリングされているヒットのスコア分布を分析**
   - 診断情報から、`cutoff_score`でフィルタリングされているヒットのスコア範囲を確認
   - 低アイデンティティヒットが多く失われているか確認

3. **低アイデンティティヒットの検出を改善**
   - `cutoff_score`や`gap_trigger`の調整を検討
   - E-value thresholdの調整を検討

## 注意事項

- 現在のLOSATの実装は、Two-hit要件を満たさない場合でもextensionを実行しているが、これはNCBI BLASTと異なる動作
- この実装により、ヒット数が2.3-4.4倍に改善しているが、NCBI BLASTとの互換性を優先するか、感度を優先するか判断が必要
- `cutoff_score`で91.8%がフィルタリングされているため、低アイデンティティヒットが多く失われている可能性が高い

