# TBLASTXセッション v12 まとめ - cutoff_score問題とsum statistics修正

## 1. 実施した修正

### 1.1 cutoff_score計算の修正 (`utils.rs` 758-786行)

**問題**: `cutoff_score_max`が`gap_trigger`に固定されていた

**修正内容**: NCBIと同じロジックに修正
- NCBI `BlastHitSavingParametersUpdate` (blast_parameters.c:943-946)を参照
- `cutoff_score_max`を`expect_value`（デフォルト10.0）から計算

```rust
// 修正後のコード
let cutoff_score_max = raw_score_from_evalue(args.evalue, &params, &search_space);

// Apply NCBI logic:
// 1. MIN(cutoff_from_evalue, gap_trigger)
let mut cutoff = cutoff_from_evalue.min(gap_trigger as i32);
// 2. *= scale_factor
cutoff = (cutoff as f64 * scale_factor) as i32;
// 3. MIN(cutoff, cutoff_score_max)
cutoff = cutoff.min(cutoff_score_max);
```

### 1.2 sum statisticsのソート順修正 (`sum_stats_linking.rs` 135-145行)

**問題**: HSPを昇順でソートしていた

**修正内容**: NCBIと同じ降順ソートに修正
- NCBI `s_RevCompareHSPsTbx` (link_hsps.c:331-376)を参照

```rust
// 修正後：降順ソート（大きいpositionが先）
sorted_indices.sort_by(|&a, &b| {
    let ha = &result_hits[a];
    let hb = &result_hits[b];
    hb.q_aa_start.cmp(&ha.q_aa_start)
        .then(hb.q_aa_end.cmp(&ha.q_aa_end))
        .then(hb.s_aa_start.cmp(&ha.s_aa_start))
        .then(hb.s_aa_end.cmp(&ha.s_aa_end))
});
```

### 1.3 リンク条件へのwindow_size制約追加 (`sum_stats_linking.rs` 355-375行)

**問題**: window_size制約がなく、遠く離れたHSPもリンク候補になっていた

**修正内容**: NCBIと同じ制約を追加
- NCBI link_hsps.c:720-736を参照

```rust
// 修正後：b4, b5制約を追加
let h_qe_gap = h_qe + WINDOW_SIZE as i32;
let h_se_gap = h_se + WINDOW_SIZE as i32;

let b1 = qo <= h_qe;     // overlap in query
let b2 = so <= h_se;     // overlap in subject
let b4 = qo > h_qe_gap;  // beyond window in query
let b5 = so > h_se_gap;  // beyond window in subject

if !(b0 || b1 || b2 || b4 || b5) {
    // Valid link candidate
}
```

### 1.4 `ignore_small_gaps`フラグの実装 (`sum_stats_linking.rs` 27-74行)

**問題**: 検索空間が小さい場合でも常にsmall gap計算を実行していた

**修正内容**: NCBIと同じ条件分岐を追加
- NCBI blast_parameters.c:1062-1078を参照

```rust
fn calculate_link_hsp_cutoffs(...) -> (i32, bool) {
    let ignore_small_gaps = search_sp <= 8.0 * window_sq;
    
    if !ignore_small_gaps {
        x_variable /= 1.0 - gap_prob + EPSILON;
    }
    
    (cutoff, ignore_small_gaps)
}
```

## 2. 現在の状況

### 2.1 ヒット数比較

| 条件 | ヒット数 | 最低bit score |
|------|---------|---------------|
| NCBI BLAST | 62,053 | 22.1 |
| LOSAT (SEGあり) | 9,838 | 22.0 |
| LOSAT (SEGなし) | 50,227 | 22.0 |

### 2.2 改善点

- **最低bit score**: 32 → 22 に改善（NCBIと同等）
- **SEGなしでのヒット数**: NCBIの約80%に到達

### 2.3 SEGフィルタの影響

```
SEG masked 543175 aa (41.33%) of 1314198 total aa across all frames
Combined DNA masks: 577369 bases (87.87%)
```

LOSATのSEGフィルタが41%のアミノ酸をマスクしており、これがヒット数減少の主因。

## 3. 残りの問題点

### 3.1 SEGフィルタの過剰マスク

- LOSATでSEGを無効にすると50,227件に増加
- NCBIもデフォルトでSEGを使用（'12 2.2 2.5'）
- しかしNCBIは62,053件を出力
- LOSATのSEGマスク実装がNCBIと異なる可能性

### 3.2 sum statisticsでのリンク形成

- リンク自体は形成されている（約96,000件のリンク）
- しかし、チェーン化されたHSPのE-valueがNCBIより大きい
- sum計算における`score - cutoff`の処理が影響している可能性

### 3.3 E-value計算の詳細な差異

同じposition、同じbit scoreのHSPで：
- NCBI: E-value 0.0（リンク後）
- LOSAT: E-value > 10（フィルタで除外）

## 4. 今後の方針

### 4.1 優先度1: SEGフィルタの調査

1. LOSATのSEGフィルタ実装（`utils.rs` 170-248行）をNCBIと比較
2. マスク率の差異の原因を特定
3. NCBIの`seg.c`との実装差を修正

### 4.2 優先度2: sum statistics E-value計算の精査

1. `BLAST_LargeGapSumE`（blast_stat.c:4532-4573）との比較
2. `blast_sum_p`計算の検証
3. リンクされたチェーンのE-value割り当て処理の確認

### 4.3 優先度3: HSP検出数の改善

1. シードマッチング処理の比較
2. 拡張処理の比較
3. フレーム処理の比較

## 5. 参照ファイル

### LOSAT
- `LOSAT/src/algorithm/tblastx/utils.rs` - メイン処理
- `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs` - sum statistics実装
- `LOSAT/src/stats/sum_statistics.rs` - E-value計算関数
C:\Users\kawato\Documents\GitHub\ncbi-blast
### NCBI BLAST
- `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_parameters.c` - cutoff計算
- `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c` - HSPリンク処理
- `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_stat.c` - E-value計算
- `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_setup.c` - 検索空間計算

