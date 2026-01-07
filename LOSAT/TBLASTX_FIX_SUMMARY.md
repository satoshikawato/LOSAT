# TBLASTX Score Calculation Parity Fix Summary

**修正日**: 2026-01-06  
**更新日**: 2026-01-06 (Long Sequence Excess HSP Fix Plan追加)

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

### Critical Issue Identified:

#### 長い配列での過剰なHSP生成（7倍の差）
- **現象**: LOSAT 338,859 HSPs vs NCBI 30,000-45,000 HSPs（約7倍）
- **特に低ビットスコア範囲（<30 bits）で過剰なHSPが生成される**
- **テストケース**: 長い配列（600kb+）

#### 根本原因分析:

1. **Cutoff計算の問題**:
   - 長い配列（600kb+）での`searchsp`計算: `MIN(query_aa, subject_nucl) * subject_nucl = 120 billion`
   - `CUTOFF_E_TBLASTX = 1e-300`とgap decay (0.5)により、実効的なcutoffが非常に低くなる
   - `cutoff_score_max`（ユーザーのE-value=10.0から計算）が41未満の場合、それが制限要因になる
   - **参照**: `ncbi_cutoffs.rs:399-433` (`cutoff_score_for_update_tblastx`)
   - **NCBI参照**: `blast_parameters.c:348-374` (`BlastInitialWordParametersUpdate`)

2. **HSP Filteringのタイミング**:
   - Extension直後のcutoffチェックは正しく実装されている
   - しかし、計算されたcutoff値が長い配列で低すぎる可能性がある
   - **参照**: `utils.rs:1950-1958` (Extension後のcutoffチェック)
   - **NCBI参照**: `aa_ungapped.c:575-591` (Extension後のcutoffチェック)

3. **Reevaluationの効果**:
   - 座標変換は正しく実装されている（`get_ungapped_hsp_list`）
   - Reevaluationでcutoffチェックが実行されている
   - しかし、低すぎるcutoffにより、低スコアHSPが通過してしまう可能性がある
   - **参照**: `utils.rs:708-757` (`reevaluate_ungapped_hsp_list`)
   - **NCBI参照**: `blast_hits.c:675-733` (`Blast_HSPReevaluateWithAmbiguitiesUngapped`)

### Next Steps:

#### Phase 1: Cutoff計算の検証と修正（最優先）
1. **Cutoff計算のデバッグ出力追加** (`ncbi_cutoffs.rs`)
   - 長い配列（600kb+）でのcutoff計算過程を記録
   - `searchsp`, `update_cutoff`, `gap_trigger`, `cutoff_score_max`, `final_cutoff`を出力
   - **参照**: `TBLASTX_LONG_SEQUENCE_HSP_FIX_PLAN.md` Phase 1.1

2. **Cutoff適用の検証** (`utils.rs`)
   - Extension直後のcutoffチェックが正しく動作しているか確認
   - 長い配列でのcutoff適用統計を追加
   - **参照**: `TBLASTX_LONG_SEQUENCE_HSP_FIX_PLAN.md` Phase 1.2

#### Phase 2: Reevaluationの検証と修正
3. **Reevaluationのcutoff適用確認** (`utils.rs`, `reevaluate.rs`)
   - Reevaluationで使用されるcutoffが正しいか確認
   - Reevaluationで削除されるHSPの統計を追加
   - **参照**: `TBLASTX_LONG_SEQUENCE_HSP_FIX_PLAN.md` Phase 2.1

4. **座標変換の再確認** (`utils.rs`)
   - `get_ungapped_hsp_list`の座標変換が正しく動作しているか確認
   - NCBIの`s_AdjustInitialHSPOffsets`と完全に一致しているか検証
   - **参照**: `TBLASTX_LONG_SEQUENCE_HSP_FIX_PLAN.md` Phase 2.2

#### Phase 3: 統計とデバッグ出力の追加
5. **HSP生成統計の追加** (`utils.rs`)
   - 長い配列でのHSP生成統計を追加
   - Cutoff適用前後のHSP数、Reevaluation前後のHSP数、スコア分布を記録
   - **参照**: `TBLASTX_LONG_SEQUENCE_HSP_FIX_PLAN.md` Phase 3.1

6. **デバッグ出力の追加** (`utils.rs`)
   - 長い配列での処理過程をデバッグ出力
   - Cutoff値、HSP数、スコア分布を出力
   - **参照**: `TBLASTX_LONG_SEQUENCE_HSP_FIX_PLAN.md` Phase 3.2

### 修正ファイル:
- `LOSAT/src/algorithm/tblastx/ncbi_cutoffs.rs`: Cutoff計算のデバッグ出力追加
- `LOSAT/src/algorithm/tblastx/utils.rs`: HSP filtering, reevaluation, 統計追加
- `LOSAT/src/algorithm/tblastx/reevaluate.rs`: Reevaluationのcutoff適用確認

### 関連ドキュメント:
- `TBLASTX_NCBI_PARITY_STATUS.md` - 詳細なステータスレポート
- `TBLASTX_NCBI_PARITY_FIX_PLAN.md` - 修正計画
- **`TBLASTX_LONG_SEQUENCE_HSP_FIX_PLAN.md`** - **長い配列での過剰なHSP生成の修正計画（新規）**

