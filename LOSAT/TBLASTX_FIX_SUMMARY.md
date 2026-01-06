# TBLASTX Score Calculation Parity Fix Summary

**修正日**: 2026-01-06

## Summary of Changes

### Fixed Issues:

#### 1. Karlinパラメータ使用に関するコメント修正
- **問題**: コメントに「cutoff score search space calculation で gapped params を使用」と誤って記載
- **修正内容**:
  - `utils.rs` のコメントを修正: tblastxはすべての計算で ungapped params を使用することを明記
  - NCBIコード参照 (`blast_setup.c:768`) を追加
  - `ncbi_cutoffs.rs` のパラメータ名を `gapped_params` → `karlin_params` に変更（汎用性を明確化）

#### 2. パラメータ使用の検証
以下の計算すべてで ungapped params を使用していることを確認:
- ✅ eff_searchsp計算: ungapped params使用
- ✅ cutoff_score_max計算: ungapped params使用  
- ✅ per-subject cutoff更新: ungapped params使用
- ✅ bit score/E-value計算: ungapped params使用

**NCBI実装確認**:
- `blast_setup.c:768`: `kbp_ptr = (scoring_options->gapped_calculation ? sbp->kbp_gap_std : sbp->kbp)`
- tblastxでは `gapped_calculation = FALSE` → `kbp_ptr = sbp->kbp` (ungapped)

### Remaining Issue:

#### ヒット数が2倍の問題
- **現象**: LOSAT 29,766 vs NCBI 14,877 (約2倍)
- **Bitscore**: LOSATがわずかに高い
- **テストケース**: AP027131.AP027133 (長い配列、600kb+)

#### 推定原因:
1. **Extension logic**: 過剰なHSP生成
2. **Cutoff application**: カットオフ適用/フィルタリングの問題
3. **Sum-statistics E-value**: E-value計算の問題

### Next Steps:

追加調査が必要な項目:
1. Extension X-drop終了条件の再検証
2. Extensionでのカットオフ適用の確認
3. Sum-statistics linkingのE-value計算の詳細確認

### 修正ファイル:
- `LOSAT/src/algorithm/tblastx/ncbi_cutoffs.rs`
- `LOSAT/src/algorithm/tblastx/utils.rs`

### 関連ドキュメント:
- `TBLASTX_NCBI_PARITY_STATUS.md` - 詳細なステータスレポート
- `TBLASTX_NCBI_PARITY_FIX_PLAN.md` - 修正計画

