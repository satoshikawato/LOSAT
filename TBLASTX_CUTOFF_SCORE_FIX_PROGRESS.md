# TBLASTX cutoff_score修正の進捗と挫折

## 作業日時
2024-12-28

## 問題の概要

LOSATのTBLASTX出力で、特に低同一性（50%未満）や低ビットスコアのヒットがNCBI BLASTと比較して大幅に少ない。
- LOSAT: 21-28%の比率（NCBI BLAST比）
- NCBI BLAST: 42-63%の比率
- 最小ビットスコア: LOSAT 32.3 vs NCBI BLAST 22.1

## 調査結果

### 診断ログからの発見

```
cutoff_from_evalue=2668  (CUTOFF_E_TBLASTX=1e-300から計算)
gap_trigger=46            (22.0 bit scoreから計算されたraw score)
cutoff_score_max=72       (args.evalue=10.0から計算、eff_searchsp使用)
final_cutoff_score=46     (最終的なcutoff_score)
limiting_factor=gap_trigger
```

**問題点**: `cutoff_score_max=72`が`gap_trigger=46`より大きいにもかかわらず、`gap_trigger`が制限要因になっていた。

## 実施した修正

### 1. 診断ログの追加 ✅
- `cutoff_score_max`、`gap_trigger`、`final_cutoff_score`の値をログ出力
- `limiting_factor`を追加して、どの値が制限要因かを識別

### 2. `cutoff_score_max`と`gap_trigger`の適用順序の修正 ✅
- **問題**: `cutoff_score_max >= gap_trigger`の場合でも、`gap_trigger`が先に適用されていた
- **修正**: `cutoff_score_max >= gap_trigger`の場合、`gap_trigger`を適用せず、`cutoff_score_max`を直接使用
- **結果**: `limiting_factor=cutoff_score_max`になったが、`cutoff_score_max=72`は依然として高すぎる

### 3. `cutoff_score_max`の計算方法の修正 ✅
- **問題**: `cutoff_score_max`が`evalue=10.0`から計算され、`gap_trigger=46`より大きくなっていた
- **修正**: `cutoff_score_max = gap_trigger`に設定（「evalue < 10 OR bitscore > threshold」のロジックを実装）
- **結果**: `cutoff_score_max=46`（`gap_trigger`と同じ）になったが、ヒット数は変わらず

## 現在の状態

### 修正後の値
```
cutoff_from_evalue=2668
gap_trigger=46
cutoff_score_max=46  (gap_triggerと同じ)
final_cutoff_score=46
limiting_factor=cutoff_score_max(gap_trigger)
```

### テスト結果
- ヒット数: 13,850（修正前と同じ）
- 最小ビットスコア: 32.3（修正前と同じ）
- NCBI BLAST: 62,053ヒット、最小ビットスコア22.1

## 挫折と未解決の問題

### 1. `cutoff_score_max`の計算方法が正しくない可能性
- **問題**: `cutoff_score_max=46`（`gap_trigger`と同じ）に設定したが、ヒット数が増えていない
- **原因の可能性**:
  - `cutoff_score_max`の計算に使用する`evalue`が異なる
  - `cutoff_score_max`の計算方法自体が間違っている
  - NCBI BLASTでは`cutoff_score_max`が`gap_trigger`より大きい場合の処理が異なる

### 2. 実装とNCBI BLASTの実装の不一致
- **NCBI BLASTの実装**: 常に`gap_trigger`を先に適用してから`cutoff_score_max`を適用
  ```c
  new_cutoff = MIN(new_cutoff, gap_trigger);
  new_cutoff *= scale_factor;
  new_cutoff = MIN(new_cutoff, cutoff_score_max);
  ```
- **LOSATの修正**: `cutoff_score_max >= gap_trigger`の場合、`gap_trigger`を適用しない
- **問題**: この修正が正しいかどうか不明

### 3. `cutoff_score_max`の値が高すぎる
- **問題**: `cutoff_score_max=72`（32.3 bit score）は、NCBI BLASTが検出する22.1 bit score（約45 raw score）のヒットを検出できない
- **試行**: `cutoff_score_max = gap_trigger`に設定したが、ヒット数が増えていない
- **原因の可能性**: `cutoff_score_max`の計算に使用する`evalue`や`search_space`が異なる

## 次のステップ（未実施）

### 1. NCBI BLASTの実装の再確認
- `cutoff_score_max`が`gap_trigger`より大きい場合の実際の動作を確認
- `cutoff_score_max`の計算に使用する`evalue`を確認
- `cutoff_score_max`が`gap_trigger`より大きい場合でも、より低スコアのヒットが検出される理由を調査

### 2. `cutoff_score_max`の計算方法の再検討
- `cutoff_score_max`を`gap_trigger`に対応するevalueから計算する
- `cutoff_score_max`を`gap_trigger`より大きく設定する（例: `gap_trigger + margin`）
- `cutoff_score_max`の計算に使用する`search_space`を再確認

### 3. 低スコアヒットの分析
- 22.1 bit score付近のヒットがなぜ検出されないか、より詳細に分析
- `cutoff_score`でフィルタリングされているヒットのスコア分布を確認
- E-valueフィルタリングの影響を確認

## 学んだこと

1. **NCBI BLASTの実装は複雑**: `cutoff_score_max`と`gap_trigger`の関係は単純ではない
2. **診断ログの重要性**: 実際の値を確認することで、問題の原因を特定できる
3. **実装の順序が重要**: `cutoff_score_max`と`gap_trigger`の適用順序が結果に大きく影響する
4. **計算式の確認**: 計算式は正しくても、使用するパラメータ（`evalue`、`search_space`）が異なると結果が変わる

## 参考資料

- NCBI BLAST実装: `ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:348-374, 920-1000`
- LOSAT実装: `LOSAT/src/algorithm/tblastx/utils.rs:450-540`
- 診断ログ: `losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.*.log`

