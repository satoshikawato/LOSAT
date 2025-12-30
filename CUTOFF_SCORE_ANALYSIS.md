# Cutoff Score 分析結果

## 診断ログからの発見

### 実際のcutoff_scoreの値

診断ログから、以下の値が確認されました：

```
cutoff_from_evalue=2668  (CUTOFF_E_TBLASTX=1e-300から計算)
gap_trigger=46            (22.0 bit scoreから計算されたraw score)
cutoff_score_max=72       (args.evalue=10.0から計算、eff_searchsp使用)
cutoff_before_max=46      (gap_triggerで制限された後の値)
final_cutoff_score=46     (最終的なcutoff_score)
limiting_factor=gap_trigger
```

### 問題の分析

1. **gap_triggerが制限要因**
   - `cutoff_score_max=72` は `gap_trigger=46` よりも大きい
   - しかし、`cutoff_before_max=46` なので、`MIN(46, 72) = 46` となり、最終的に `gap_trigger=46` が制限要因
   - **すべてのcontextで `limiting_factor=gap_trigger`**

2. **cutoff_score_maxの役割**
   - `cutoff_score_max=72` は `gap_trigger=46` よりも大きいため、実際には制限要因になっていない
   - これは、`cutoff_score_max` が `gap_trigger` よりも小さい場合にのみ有効になることを意味する

3. **低スコアヒットのフィルタリング**
   - 現在の `final_cutoff_score=46` (raw score) は約22.0 bit scoreに相当
   - NCBI BLASTは22.1 bit scoreのヒットを検出しているが、LOSATは32.3 bit score以上のみ検出
   - これは、`cutoff_score=46` が高すぎることを示している

### NCBI BLASTの実装との比較

NCBI BLASTの実装（`blast_parameters.c:348-374`）:
```c
if (program_number != eBlastTypeBlastn)  
    new_cutoff = MIN(new_cutoff, gap_trigger);
new_cutoff *= (Int4)sbp->scale_factor;
new_cutoff = MIN(new_cutoff, 
                 hit_params->cutoffs[context].cutoff_score_max);
```

現在のLOSAT実装:
```rust
let mut cutoff = cutoff_from_evalue.min(gap_trigger);  // 46
cutoff = (cutoff as f64 * scale_factor) as i32;        // 46 * 1.0 = 46
cutoff = cutoff.min(cutoff_score_max);                 // MIN(46, 72) = 46
```

**実装は正しい**が、`cutoff_score_max` が `gap_trigger` よりも大きい場合、`gap_trigger` が制限要因になる。

### 解決策の検討

1. **cutoff_score_maxがgap_triggerよりも小さい場合の確認**
   - より小さなsearch spaceやより大きなevalue thresholdの場合、`cutoff_score_max` が `gap_trigger` よりも小さくなる可能性がある
   - その場合、`cutoff_score_max` が制限要因となり、より低スコアのヒットが通過できる

2. **gap_triggerの適用条件の再検討**
   - NCBI BLASTでは、`cutoff_score_max` が `gap_trigger` よりも小さい場合、`cutoff_score_max` を使用する
   - しかし、現在の実装では、`gap_trigger` が先に適用されるため、`cutoff_score_max` が小さくても効果がない

3. **実際の問題**
   - 現在のテストケースでは、`cutoff_score_max=72 > gap_trigger=46` のため、`gap_trigger` が制限要因
   - より小さなsearch spaceやより大きなevalue thresholdを使用すれば、`cutoff_score_max` が小さくなり、より多くのヒットが通過できる可能性がある

### 次のステップ

1. **cutoff_score_maxの計算を確認**
   - `eff_searchsp` を使用した計算が正しいか確認
   - より小さなsearch spaceでの動作を確認

2. **gap_triggerの適用順序の再検討**
   - NCBI BLASTの実装を再確認し、`cutoff_score_max` が `gap_trigger` よりも小さい場合の動作を確認

3. **低スコアヒットの分析**
   - 22.1 bit score付近のヒットがなぜ検出されないか、より詳細に分析


