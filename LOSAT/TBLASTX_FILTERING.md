# TBLASTXで実装されているフィルタリング

このドキュメントは、TBLASTXで実装されているすべてのフィルタリングメカニズムをまとめたものです。

## 1. Lookup Table構築時のフィルタリング

### 1.1 SEGマスキング
- **実装場所**: `lookup.rs:337-342`
- **説明**: SEGマスクされた領域のwordはスキップされる
- **NCBI同等**: ✅ NCBIと同等
- **コード**:
```rust
let unmasked = compute_unmasked_intervals(&frame.seg_masks, frame.aa_len);
for (start, end) in unmasked {
    for aa_pos in start..end.saturating_sub(word_length - 1) {
        // wordを追加
    }
}
```

### 1.2 無効な残基（Invalid Residue）のスキップ
- **実装場所**: `lookup.rs:352-355`
- **説明**: アルファベットサイズ（28）を超える残基を含むwordはスキップ
- **NCBI同等**: ✅ NCBIと同等
- **コード**:
```rust
if c0 >= alphabet_size || c1 >= alphabet_size || c2 >= alphabet_size {
    skipped_invalid_residue += 1;
    continue;
}
```

### 1.3 Thresholdベースのフィルタリング
- **実装場所**: `lookup.rs:472-515`
- **説明**: Wordのself-scoreがthreshold未満の場合のみ、exact matchをlookup tableに追加
- **NCBI同等**: ✅ NCBIと同等（`blast_aalookup.c:504-514`）
- **コード**:
```rust
if lookup->threshold == 0 || score < lookup->threshold {
    // exact matchを追加
} else {
    // neighbor wordsを探索（exact matchは追加しない）
}
```

## 2. Seeding段階のフィルタリング

### 2.1 Two-Hit Filtering
- **実装場所**: `utils.rs:854-888`
- **説明**: 2つのhitがwindow内に存在する場合のみextensionを実行
- **NCBI同等**: ✅ NCBIと同等（`aa_ungapped.c:s_BlastAaWordFinder_TwoHit`）
- **フィルタリング条件**:
  1. `diff >= window` (40 AA): 2つ目のhitがwindow外 → スキップ
  2. `diff < wordsize` (3 AA): 2つのhitが重複 → スキップ
  3. `query_offset - diff < ctx.frame_base`: 1つ目のhitがcontext外 → スキップ
- **コード**:
```rust
// [C] if (diff >= window)
if diff >= window {
    diag_entry.last_hit = subject_offset + diag_offset;
    continue;
}

// [C] if (diff < wordsize)
if diff < wordsize {
    continue;
}

// [C] if (query_offset - diff < query_info->contexts[curr_context].query_offset)
if query_offset - diff < ctx.frame_base {
    diag_entry.last_hit = subject_offset + diag_offset;
    continue;
}
```

### 2.2 Masked Region Filtering
- **実装場所**: `utils.rs:831-845`
- **説明**: マスクされた領域のhitはスキップ
- **NCBI同等**: ✅ NCBIと同等
- **コード**:
```rust
if diag_entry.flag != 0 {
    if subject_offset + diag_offset < diag_entry.last_hit {
        // マスクされた領域 → スキップ
        continue;
    }
    diag_entry.flag = 0;
}
```

## 3. Extension段階のフィルタリング

### 3.1 X-Drop終了条件
- **実装場所**: `extension.rs:228-293`
- **説明**: X-drop値に達した場合、extensionを終了
- **NCBI同等**: ✅ NCBIと同等（`aa_ungapped.c:s_BlastAaExtendTwoHit`）
- **終了条件**: `(max_score - current_score) >= x_drop`
- **コード**:
```rust
if (max_score - current_score) >= x_drop {
    break; // extension終了
}
```

### 3.2 Cutoff Score Filtering
- **実装場所**: `utils.rs:973`
- **説明**: Extension後のscoreがcutoff未満の場合、HSPを保存しない
- **NCBI同等**: ✅ NCBIと同等（`aa_ungapped.c:403`）
- **コード**:
```rust
// [C] if (score >= cutoffs->cutoff_score)
if score >= cutoff {
    // HSPを保存
} else {
    // スキップ
}
```

### 3.3 Cutoff計算
- **実装場所**: `utils.rs:732-749`, `ncbi_cutoffs.rs:385-419`
- **説明**: 最終的なcutoffは`MIN(update_cutoff, gap_trigger, cutoff_score_max)`で決定
- **NCBI同等**: ✅ NCBIと同等（`blast_parameters.c:348-374`）
- **計算式**:
  1. `cutoff_score_max`: `eff_searchsp`から計算（ユーザーのE-value使用）
  2. `update_cutoff`: `CUTOFF_E_TBLASTX=1e-300`から計算
  3. `gap_trigger`: 固定値（BLOSUM62の場合、raw score = 41）
  4. 最終cutoff = `MIN(update_cutoff, gap_trigger, cutoff_score_max)`

## 4. Reevaluation段階のフィルタリング

### 4.1 Ungapped HSP Reevaluation
- **実装場所**: `reevaluate.rs:80-145`
- **説明**: Extension後のHSPを実際の（マスクされた可能性のある）配列に対して再評価
- **NCBI同等**: ✅ NCBIと同等（`blast_hits.c:Blast_HSPReevaluateWithAmbiguitiesUngapped`）
- **フィルタリング条件**:
  - Reevaluation後のscoreがcutoff未満の場合、HSPを削除
- **コード**:
```rust
if let Some((new_qs, new_ss, new_len, new_score)) =
    reevaluate_ungapped_hit_ncbi_translated(query, subject, qs, ss, len_u, cutoff)
{
    // HSPを更新
} else {
    // Deleted by NCBI reevaluation (score < cutoff_score)
    continue;
}
```

## 5. Subject Scanning段階のフィルタリング

### 5.1 Offset Array Size制限
- **実装場所**: `utils.rs:638-639`
- **説明**: `offset_array_size = OFFSET_ARRAY_SIZE (4096) + lookup.longest_chain`
- **NCBI同等**: ✅ NCBIと同等（`blast_aascan.c`）
- **説明**: 一度に処理するhit数の上限を設定

### 5.2 Early Termination条件
- **実装場所**: `utils.rs:807-811`
- **説明**: 連続してhitsが0の場合、スキャンを終了
- **NCBI同等**: ✅ NCBIと同等

## 6. 実装されていない、またはタイミングが異なるフィルタリング

### 6.1 max_hits_per_kmerフィルタリング
- **状態**: ❌ **削除済み**（2026-01-05）
- **説明**: Over-represented k-mersのフィルタリングはNCBIには存在しないため、削除されました

### 6.2 Blast_HSPTest（Percent Identity / Min Hit Length）
- **状態**: ✅ **実装済み**（2026-01-05）
- **説明**: NCBIの`Blast_HSPTest`は`Blast_HSPListReevaluateUngapped`内で呼び出される（`blast_hits.c:2719`）
- **NCBIコード**: `blast_hits.c:993-999`
- **NCBIの動作**:
  ```c
  Blast_HSPGetNumIdentitiesAndPositives(query_nomask, subject_start, hsp, ...);
  delete_hsp = Blast_HSPTest(hsp, hit_params->options, align_length);
  ```
- **LOSATの実装**:
  - `get_num_identities_and_positives_ungapped` (`reevaluate.rs:214-242`): NCBIの`Blast_HSPGetNumIdentitiesAndPositives`相当
  - `hsp_test` (`reevaluate.rs:269-287`): NCBIの`s_HSPTest`相当
  - Reevaluation後に呼び出される（`utils.rs:1000-1020`, `utils.rs:1864-1890`）
- **デフォルト値**: `percent_identity = 0.0`, `min_hit_length = 0`（NCBIと同等、実質的に無効）
- **オプション**: `--percent-identity`と`--min-hit-length`で設定可能（`args.rs`）

### 6.3 Reevaluationのタイミングの違い
- **NCBI**: Extension → Save → **一括Reevaluate** (`Blast_HSPListReevaluateUngapped`) → `Blast_HSPTest`
- **LOSAT**: Extension → **Reevaluate** → Save
- **説明**: タイミングは異なるが、機能的には同等。ただし、NCBIは一括処理で`Blast_HSPTest`も実行する

## 7. まとめ

TBLASTXで実装されているフィルタリングは以下の通りです：

1. **Lookup Table構築時**:
   - SEGマスキング ✅ NCBI同等
   - 無効な残基のスキップ ✅ NCBI同等
   - Thresholdベースのフィルタリング ✅ NCBI同等

2. **Seeding段階**:
   - Two-hit filtering ✅ NCBI同等
   - Masked region filtering ✅ NCBI同等

3. **Extension段階**:
   - X-drop終了条件 ✅ NCBI同等
   - Cutoff score filtering ✅ NCBI同等

4. **Reevaluation段階**:
   - Ungapped HSP reevaluation ✅ NCBI同等（タイミングは異なるが機能的には同等）
   - **Blast_HSPGetNumIdentitiesAndPositives**: ✅ **実装済み**（`reevaluate.rs:214-250`）
   - **Blast_HSPTest**: ✅ **実装済み**（`reevaluate.rs:252-275`）

5. **Subject Scanning段階**:
   - Offset array size制限 ✅ NCBI同等
   - Early termination条件 ✅ NCBI同等

**重要な点**:
- **すべてのフィルタリングはNCBIと同等に実装されています**
- `Blast_HSPTest`と`Blast_HSPGetNumIdentitiesAndPositives`は実装済み（2026-01-05）
- デフォルトでは`percent_identity = 0.0`, `min_hit_length = 0`のため、実質的には無効（NCBIと同等）
- Cutoff計算は`MIN(update_cutoff, gap_trigger, cutoff_score_max)`で決定されます
- Reevaluationのタイミングは異なるが、機能的には同等

## 8. 長い配列での過剰なHSP生成の問題（2026-01-06追加）

### 問題の概要
- **現象**: 長い配列（600kb+）でLOSATが約338,859 HSPsを生成 vs NCBIの30,000-45,000（約7倍の過剰生成）
- **特に低ビットスコア範囲（<30 bits）で過剰なHSPが生成される**

### 根本原因

#### 8.1 Cutoff計算の問題
- **実装場所**: `ncbi_cutoffs.rs:399-433` (`cutoff_score_for_update_tblastx`)
- **NCBI参照**: `blast_parameters.c:348-374` (`BlastInitialWordParametersUpdate`)
- **問題点**:
  1. 長い配列（600kb+）での`searchsp`計算: `MIN(query_aa, subject_nucl) * subject_nucl = 120 billion`
  2. `CUTOFF_E_TBLASTX = 1e-300`とgap decay (0.5)により、実効的なcutoffが非常に低くなる
  3. `cutoff_score_max`（ユーザーのE-value=10.0から計算）が41未満の場合、それが制限要因になる
  4. 結果として、低スコアHSP（<30 bits）が通過してしまう

#### 8.2 HSP Filteringのタイミング
- **実装場所**: `utils.rs:1950-1958` (Extension後のcutoffチェック)
- **NCBI参照**: `aa_ungapped.c:575-591` (Extension後のcutoffチェック)
- **問題点**:
  - Extension直後のcutoffチェックは正しく実装されている
  - しかし、計算されたcutoff値が長い配列で低すぎる可能性がある
  - 低スコアHSPがcutoffを通過してしまう

#### 8.3 Reevaluationの効果
- **実装場所**: `utils.rs:708-757` (`reevaluate_ungapped_hsp_list`)
- **NCBI参照**: `blast_hits.c:675-733` (`Blast_HSPReevaluateWithAmbiguitiesUngapped`)
- **問題点**:
  - 座標変換は正しく実装されている（`get_ungapped_hsp_list`）
  - Reevaluationでcutoffチェックが実行されている
  - しかし、低すぎるcutoffにより、低スコアHSPが通過してしまう可能性がある

### 修正計画
詳細は `TBLASTX_LONG_SEQUENCE_HSP_FIX_PLAN.md` を参照してください。

**主な修正項目**:
1. **Phase 1**: Cutoff計算の検証と修正（最優先）
   - Cutoff計算のデバッグ出力追加
   - Cutoff適用の検証
2. **Phase 2**: Reevaluationの検証と修正
   - Reevaluationのcutoff適用確認
   - 座標変換の再確認
3. **Phase 3**: 統計とデバッグ出力の追加
   - HSP生成統計の追加
   - デバッグ出力の追加

### 期待される結果
修正後、以下の改善が期待されます：
1. **HSP生成数の削減**: 338,859 → 30,000-45,000（NCBIと同等）
2. **低ビットスコアHSPの削減**: <30 bitsのHSPが適切にフィルタリングされる
3. **Cutoff計算の正確性**: 長い配列でも正しいcutoffが計算される
4. **NCBIとの完全一致**: HSP生成数、スコア分布がNCBIと一致する

