# TBLASTX NCBI BLAST パリティ改修計画

## 調査日: 2026-01-06

## 概要

blastnの指摘（off-diagonal checkingの欠如、HSP cullingの欠如）を参考に、LOSAT TBLASTXとNCBI TBLASTXの実装上の違いを網羅的に確認し、NCBI BLASTに完全に一致させるための改修計画をまとめる。

## ⚠️ 最優先改修項目

**Post-Processing HSP Filtering (HSP Culling)** を最優先で実装する必要がある。

---

## 1. Post-Processing HSP Filtering (HSP Culling) の欠如 ⚠️ 最優先

### 問題点

**NCBI BLAST** (`hspfilter_culling.c:78-120`):
- Interval treeベースのHSP cullingアルゴリズムを実装
- `s_DominateTest`関数で支配判定を行う
- スコア/長さのトレードオフ公式: `4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2`
- 50%以上のオーバーラップが必要
- Query座標のみを考慮（subject座標は考慮しない）
- Interval treeを使用して効率的に実装（O(log n)検索）

**LOSAT TBLASTX**:
- 現在、HSP cullingは実装されていない
- `chaining.rs:262-265`でdomination filterをスキップしている
- コメント: "NCBI BLAST does not apply HSP culling for TBLASTX by default"
- **しかし、これは誤解**: NCBI BLASTは`--culling_limit`オプションでcullingを有効化可能
- デフォルトでは無効だが、実装自体は存在する

### NCBI実装詳細

**Dominance Test** (`hspfilter_culling.c:78-120`):
```c
// NCBI reference: hspfilter_culling.c:78-120
static Boolean s_DominateTest(LinkedHSP *p, LinkedHSP *y) {
    Int8 b1 = p->begin;
    Int8 b2 = y->begin;
    Int8 e1 = p->end;
    Int8 e2 = y->end;
    Int8 s1 = p->hsp->score;
    Int8 s2 = y->hsp->score;
    Int8 l1 = e1 - b1;
    Int8 l2 = e2 - b2;
    Int8 overlap = MIN(e1,e2) - MAX(b1,b2);
    Int8 d = 0;

    // If not overlap by more than 50%
    if(2 *overlap < l2) {
    	return FALSE;
    }

    /* the main criterion:
       2 * (%diff in score) + 1 * (%diff in length) */
    d  = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2;
    
    // If identical, use oid as tie breaker
    if(((s1 == s2) && (b1==b2) && (l1 == l2)) || (d == 0)) {
    	if(s1 != s2) {
    		return (s1>s2);
    	}
    	if(p->sid != y->sid) {
    		return (p->sid < y->sid);
    	}
    	if(p->hsp->subject.offset > y->hsp->subject.offset) {
    		return FALSE;
    	}
    	return TRUE;
    }

   	if (d < 0) {
   		return FALSE;
    }

    return TRUE;
}
```

**Interval Tree実装** (`hspfilter_culling.c:200-470`):
- `CTreeNode`構造体でinterval treeを実装
- Query座標範囲（`begin`, `end`）で分割
- `s_SaveHSP`関数でHSPを追加し、meritを計算
- `s_ProcessCTree`関数で再帰的に支配関係を更新

### 改修内容

1. **HSP Cullingモジュールの実装** (`src/algorithm/tblastx/hsp_culling.rs`):
   - `LinkedHSP`構造体の実装
   - `dominate_test`関数の実装（NCBIの`s_DominateTest`と同等）
   - Interval tree (`CTreeNode`)の実装
   - `save_hsp`関数の実装（NCBIの`s_SaveHSP`と同等）
   - `process_ctree`関数の実装（NCBIの`s_ProcessCTree`と同等）

2. **統合** (`utils.rs`):
   - Sum-stats linkingの後にHSP cullingを適用
   - `--culling_limit`オプションの追加（デフォルト: 0 = 無効）
   - NCBI BLASTと同様に、デフォルトでは無効だが実装は存在

3. **NCBIコード参照の追加**:
   - `hspfilter_culling.c:78-120`をコメントとして引用
   - `hspfilter_culling.c:200-470`をコメントとして引用

### 影響範囲

- `src/algorithm/tblastx/hsp_culling.rs`: 新規作成（HSP culling実装）
- `src/algorithm/tblastx/utils.rs`: HSP cullingの統合
- `src/algorithm/tblastx/args.rs`: `--culling_limit`オプションの追加
- `src/algorithm/tblastx/chaining.rs`: Domination filterの削除（cullingに置き換え）

---

## 2. Karlinパラメータ使用の明確化 ✅ 完了 (2026-01-06)

### 問題点

**旧LOSAT**:
- コメントに「cutoff score search space calculation で gapped params を使用」と誤って記載
- パラメータ名が `gapped_params` で誤解を招く可能性

**NCBI実装**:
- `blast_setup.c:768`: `kbp_ptr = (scoring_options->gapped_calculation ? sbp->kbp_gap_std : sbp->kbp)`
- tblastxでは `gapped_calculation = FALSE` のため、**すべての計算で ungapped params (`sbp->kbp`) を使用**

### 修正内容

1. **コメント修正** (`utils.rs`):
   - tblastxはすべての計算（eff_searchsp, cutoff, bit score, E-value）で ungapped params を使用することを明記
   - NCBIコード参照 (`blast_setup.c:768`) を追加

2. **パラメータ名の明確化** (`ncbi_cutoffs.rs`):
   - `gapped_params` → `karlin_params` に変更
   - 関数シグネチャとコメントを更新して、パラメータが gapped/ungapped のどちらでも使用可能であることを明記

### 検証結果

- eff_searchsp計算: ungapped params使用 ✅
- cutoff_score_max計算: ungapped params使用 ✅
- per-subject cutoff更新: ungapped params使用 ✅
- bit score/E-value計算: ungapped params使用 ✅

### 残存問題

- ヒット数が2倍（29,766 vs 14,877）の問題は依然として未解決
- 原因は拡張ロジック、カットオフ適用、またはSum-statistics E-value計算にある可能性

### 実装詳細

**Dominance Test Formula**:
```
d = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2
```
- `s1`, `s2`: HSPスコア
- `l1`, `l2`: HSP長さ（query座標での長さ）
- `d >= 0`: pがyを支配する
- `d < 0`: pがyを支配しない

**Overlap Check**:
```
overlap = MIN(e1, e2) - MAX(b1, b2)
if (2 * overlap < l2) return FALSE;  // 50%未満のオーバーラップ
```

**Interval Tree Structure**:
- Query座標範囲で分割
- 各ノードにHSPリストを保持
- Merit（支配される回数）を追跡
- Meritが0になったHSPを削除

---

## 2. One-Hit Mode実装の欠如

### 問題点

**NCBI BLAST** (`aa_ungapped.c:200-236`):
- `BlastAaWordFinder`は`ewp->diag_table->multiple_hits`フラグでone-hit/two-hitを切り替える
- One-hit mode: `s_BlastAaWordFinder_OneHit` (lines 712-812)
- Two-hit mode: `s_BlastAaWordFinder_TwoHit` (lines 439-619)

**LOSAT TBLASTX**:
- Two-hit modeのみ実装 (`utils.rs:2800-3000`)
- One-hit mode (`window_size = 0`)が実装されていない
- `--window-size 0`オプションが存在するが、実際にはtwo-hitロジックが実行される

### NCBI実装詳細

**One-Hit Mode** (`aa_ungapped.c:712-812`):
```c
// NCBI reference: aa_ungapped.c:776-802
diag_coord = (subject_offset - query_offset) & diag_mask;
diff = subject_offset - (diag_array[diag_coord].last_hit - diag_offset);

/* do an extension, but only if we have not already extended this far */
if (diff >= 0) {
    Int4 curr_context = BSearchContextInfo(query_offset, query_info);
    BlastUngappedCutoffs *cutoffs = word_params->cutoffs + curr_context;
    score = s_BlastAaExtendOneHit(matrix, subject, query,
                                  subject_offset, query_offset,
                                  cutoffs->x_dropoff, 
                                  &hsp_q, &hsp_s, &hsp_len,
                                  wordsize, use_pssm, &s_last_off);

    /* if the hsp meets the score threshold, report it */
    if (score >= cutoffs->cutoff_score) {
        BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                         query_offset, subject_offset, hsp_len,
                         score);
    }
    diag_array[diag_coord].last_hit =
        s_last_off - (wordsize - 1) + diag_offset;
    ++hits_extended;
}
```

**重要なポイント**:
- `diff >= 0`の時のみextensionを実行（既にextensionした範囲を超えている場合のみ）
- Two-hitチェックなし（window_sizeチェックなし）
- `s_BlastAaExtendOneHit`を使用（two-hit extensionではない）

### 改修内容

1. **One-Hit Mode実装の追加** (`utils.rs`):
   - `window_size == 0`の場合、one-hit modeロジックを実行
   - `diff >= 0`チェックのみ（two-hitチェックなし）
   - `extend_hit_ungapped`を使用（two-hit extensionではない）

2. **NCBIコード参照の追加**:
   - `aa_ungapped.c:712-812`をコメントとして引用

### 影響範囲

- `src/algorithm/tblastx/utils.rs`: `run_with_neighbor_map`関数内のseeding/extensionロジック
- `src/algorithm/tblastx/extension.rs`: `extend_hit_ungapped`がone-hit extensionとして正しく動作するか確認

---

## 3. Off-Diagonal Checking（TBLASTXには不要）

### 確認結果

**blastn** (`na_ungapped.c:685-717`):
- Off-diagonal checkingが存在（single wordの場合）
- `Delta = MIN(scan_range, window_size - word_length)`で範囲を計算
- 近接する対角線上のhitをチェックしてextensionをトリガー

**TBLASTX** (`aa_ungapped.c`):
- **Off-diagonal checkingは存在しない**
- Protein searchでは不要（blastn特有の機能）

### 結論

TBLASTXではoff-diagonal checkingは実装不要。これは正しい。

---

## 4. Diagonal Trackingの実装確認

### NCBI実装

**Two-Hit Mode** (`aa_ungapped.c:518-606`):
```c
// NCBI reference: aa_ungapped.c:518-530
/* If the reset bit is set, an extension just happened. */
if (diag_array[diag_coord].flag) {
    /* If we've already extended past this hit, skip it. */
    if ((Int4) (subject_offset + diag_offset) <
        diag_array[diag_coord].last_hit) {
        continue;
    }
    /* Otherwise, start a new hit. */
    else {
        diag_array[diag_coord].last_hit =
            subject_offset + diag_offset;
        diag_array[diag_coord].flag = 0;
    }
}
/* If the reset bit is cleared, try to start an extension. */
else {
    /* find the distance to the last hit on this diagonal */
    last_hit = diag_array[diag_coord].last_hit - diag_offset;
    diff = subject_offset - last_hit;

    if (diff >= window) {
        /* We are beyond the window for this diagonal; start a new hit */
        diag_array[diag_coord].last_hit =
            subject_offset + diag_offset;
        continue;
    }

    if (diff < wordsize) {
        /* If the difference is less than the wordsize (i.e. last hit and this hit overlap), give up */
        continue;
    }

    // ... two-hit extension ...
}
```

**One-Hit Mode** (`aa_ungapped.c:776-802`):
```c
// NCBI reference: aa_ungapped.c:776-802
diag_coord = (subject_offset - query_offset) & diag_mask;
diff = subject_offset - (diag_array[diag_coord].last_hit - diag_offset);

/* do an extension, but only if we have not already extended this far */
if (diff >= 0) {
    // ... one-hit extension ...
    diag_array[diag_coord].last_hit =
        s_last_off - (wordsize - 1) + diag_offset;
}
```

### LOSAT実装確認

**Two-Hit Mode** (`utils.rs:2829-2917`):
- NCBIと一致している（コメントで参照も記載済み）

**One-Hit Mode**:
- **実装されていない**

### 改修内容

1. **One-Hit Modeのdiagonal tracking実装**:
   - `diff >= 0`チェックのみ
   - `last_hit`更新: `s_last_off - (wordsize - 1) + diag_offset`
   - Two-hitチェック（window, wordsize overlap）をスキップ

---

## 5. Extension関数の確認

### NCBI実装

**One-Hit Extension** (`aa_ungapped.c:151-162`):
```c
// NCBI reference: aa_ungapped.c:151-162
static Int4 s_BlastAaExtendOneHit(Int4 ** matrix,
                         const BLAST_SequenceBlk* subject,
                         const BLAST_SequenceBlk* query,
                         Int4 s_off,
                         Int4 q_off,
                         Int4 dropoff,
                         Int4* hsp_q,
                         Int4* hsp_s,
                         Int4* hsp_len,
                         Int4 word_size,
                         Boolean use_pssm,
                         Int4* s_last_off);
```

**Two-Hit Extension** (`aa_ungapped.c:184-197`):
```c
// NCBI reference: aa_ungapped.c:184-197
static Int4 s_BlastAaExtendTwoHit(Int4 ** matrix,
                         const BLAST_SequenceBlk* subject,
                         const BLAST_SequenceBlk* query,
                         Int4 s_left_off,
                         Int4 s_right_off,
                         Int4 q_right_off,
                         Int4 dropoff,
                         Int4* hsp_q,
                         Int4* hsp_s,
                         Int4* hsp_len,
                         Boolean use_pssm,
                         Int4 word_size,
                         Boolean *right_extend,
                         Int4* s_last_off);
```

### LOSAT実装確認

**`extend_hit_ungapped`** (`extension.rs:409-500`):
- One-hit extensionとして実装されている
- NCBIの`s_BlastAaExtendOneHit`と同等

**`extend_hit_two_hit`** (`extension.rs:502-580`):
- Two-hit extensionとして実装されている
- NCBIの`s_BlastAaExtendTwoHit`と同等

### 結論

Extension関数は正しく実装されている。One-hit modeのseedingロジックのみ追加すればよい。

---

## 6. Context Managementの確認

### NCBI実装

**Context検索** (`aa_ungapped.c:560, 783`):
```c
// NCBI reference: aa_ungapped.c:560, 783
Int4 curr_context = BSearchContextInfo(query_offset, query_info);
BlastUngappedCutoffs *cutoffs = word_params->cutoffs + curr_context;
```

### LOSAT実装確認

**Neighbor-Map Mode** (`utils.rs:2814`):
- `ctx_flat = ctx_base[q_idx as usize] + q_f_idx as usize;`
- 各`(q_idx, q_f_idx)`が独立したcontextとして扱われている
- NCBIの`BSearchContextInfo`と同等の動作

### 結論

Context managementは正しく実装されている。

---

## 7. Cutoff Score適用の確認

### NCBI実装

**Two-Hit Mode** (`aa_ungapped.c:588-591`):
```c
// NCBI reference: aa_ungapped.c:588-591
if (score >= cutoffs->cutoff_score)
    BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                     query_offset, subject_offset, hsp_len,
                     score);
```

**One-Hit Mode** (`aa_ungapped.c:794-797`):
```c
// NCBI reference: aa_ungapped.c:794-797
if (score >= cutoffs->cutoff_score) {
    BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                     query_offset, subject_offset, hsp_len,
                     score);
}
```

### LOSAT実装確認

**Two-Hit Mode** (`utils.rs:2922-2925`):
```rust
let cutoff = cutoff_scores_ref[ctx_flat][s_f_idx];
if score < cutoff {
    continue;
}
```

### 結論

Cutoff score適用は正しく実装されている。One-hit modeでも同じロジックを使用すればよい。

---

## 8. 改修実装計画

### Phase 1: Post-Processing HSP Filtering (HSP Culling) 実装 ⚠️ 最優先

**ファイル**: `src/algorithm/tblastx/utils.rs`

**変更箇所**: `run_with_neighbor_map`関数内のseeding/extensionロジック

**実装内容**:

1. **Window sizeチェックの追加**:
   ```rust
   let window = args.window_size as i32;
   let is_one_hit_mode = window == 0;
   ```

2. **One-Hit Modeロジックの実装**:
   ```rust
   if is_one_hit_mode {
       // NCBI reference: aa_ungapped.c:776-802 (s_BlastAaWordFinder_OneHit)
       // One-hit mode: extend immediately if diff >= 0
       let last_hit = diag_entry.last_hit - diag_offset;
       let diff = subject_offset - last_hit;
       
       // NCBI: "do an extension, but only if we have not already extended this far"
       if diff >= 0 {
           // ... one-hit extension ...
       }
   } else {
       // Existing two-hit mode logic
       // ...
   }
   ```

3. **NCBIコード参照の追加**:
   - `aa_ungapped.c:712-812` (s_BlastAaWordFinder_OneHit)をコメントとして引用
   - `aa_ungapped.c:776-802` (one-hit extension logic)をコメントとして引用

### Phase 2: One-Hit Mode実装

**ファイル**: `src/algorithm/tblastx/utils.rs`

**変更箇所**: `run_with_neighbor_map`関数内のseeding/extensionロジック

**実装内容**:

1. **Window sizeチェックの追加**:
   ```rust
   let window = args.window_size as i32;
   let is_one_hit_mode = window == 0;
   ```

2. **One-Hit Modeロジックの実装**:
   ```rust
   if is_one_hit_mode {
       // NCBI reference: aa_ungapped.c:776-802 (s_BlastAaWordFinder_OneHit)
       // One-hit mode: extend immediately if diff >= 0
       let last_hit = diag_entry.last_hit - diag_offset;
       let diff = subject_offset - last_hit;
       
       // NCBI: "do an extension, but only if we have not already extended this far"
       if diff >= 0 {
           // ... one-hit extension ...
       }
   } else {
       // Existing two-hit mode logic
       // ...
   }
   ```

3. **NCBIコード参照の追加**:
   - `aa_ungapped.c:712-812`をコメントとして引用
   - `aa_ungapped.c:776-802`をコメントとして引用

### Phase 3: テスト

1. **One-Hit Modeテスト**:
   - `--window-size 0`で実行
   - NCBI BLAST+の`-window_size 0`と出力を比較

2. **Two-Hit Mode回帰テスト**:
   - デフォルト（`window_size=40`）で実行
   - 既存のテストが通ることを確認

### Phase 4: ドキュメント更新

1. **`TBLASTX_NCBI_PARITY_STATUS.md`の更新**:
   - One-hit mode実装完了を記録
   - NCBIコード参照を追加

---

## 9. 実装詳細

### HSP Culling実装詳細

**新規ファイル**: `src/algorithm/tblastx/hsp_culling.rs`

```rust
// NCBI reference: hspfilter_culling.c:52-60
struct LinkedHSP {
    hsp: UngappedHit,
    ctx_idx: usize,
    s_idx: u32,
    begin: i32,  // query offset in plus strand
    end: i32,    // query end in plus strand
    merit: i32,  // how many other hsps in the tree dominates me?
    next: Option<Box<LinkedHSP>>,
}

// NCBI reference: hspfilter_culling.c:78-120
fn dominate_test(p: &LinkedHSP, y: &LinkedHSP) -> bool {
    let b1 = p.begin;
    let b2 = y.begin;
    let e1 = p.end;
    let e2 = y.end;
    let s1 = p.hsp.raw_score as i64;
    let s2 = y.hsp.raw_score as i64;
    let l1 = e1 - b1;
    let l2 = e2 - b2;
    let overlap = (e1.min(e2) - b1.max(b2)).max(0);
    
    // NCBI: "If not overlap by more than 50%"
    if 2 * overlap < l2 {
        return false;
    }
    
    // NCBI: "the main criterion: 2 * (%diff in score) + 1 * (%diff in length)"
    // d = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2
    let d = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2;
    
    // NCBI: "If identical, use oid as tie breaker"
    if (s1 == s2 && b1 == b2 && l1 == l2) || d == 0 {
        if s1 != s2 {
            return s1 > s2;
        }
        if p.s_idx != y.s_idx {
            return p.s_idx < y.s_idx;
        }
        if p.hsp.s_aa_start > y.hsp.s_aa_start {
            return false;
        }
        return true;
    }
    
    if d < 0 {
        return false;
    }
    
    true
}

// NCBI reference: hspfilter_culling.c:200-207
struct CTreeNode {
    begin: i32,   // left endpoint
    end: i32,     // right endpoint
    left: Option<Box<CTreeNode>>,
    right: Option<Box<CTreeNode>>,
    hsplist: Option<Box<LinkedHSP>>,  // hsps belong to this node, start with low merits
}

// NCBI reference: hspfilter_culling.c:430-470
fn save_hsp(tree: &mut CTreeNode, a: &LinkedHSP) -> bool {
    // Descend the tree and check domination
    // If valid, insert and update merit
    // ...
}
```

**統合** (`utils.rs`):
```rust
// After sum-stats linking
let linked_hits = apply_sum_stats_even_gap_linking(...);

// Apply HSP culling if enabled
let culled_hits = if args.culling_limit > 0 {
    apply_hsp_culling(linked_hits, args.culling_limit, query_contexts)
} else {
    linked_hits
};
```

### One-Hit Mode Seeding Logic

```rust
// NCBI reference: aa_ungapped.c:712-812 (s_BlastAaWordFinder_OneHit)
// NCBI reference: aa_ungapped.c:776-802 (one-hit extension logic)
if is_one_hit_mode {
    // NCBI: diag_coord = (subject_offset - query_offset) & diag_mask
    // We use direct index: diag_coord = (query_offset - subject_offset + diag_offset) as usize
    let diag_coord = (query_offset - subject_offset + diag_offset) as usize;
    
    if diag_coord >= diag_array_len {
        continue;
    }
    
    let diag_entry = &mut diag_array[ctx_flat][diag_coord];
    
    // NCBI: diff = subject_offset - (diag_array[diag_coord].last_hit - diag_offset)
    let last_hit = diag_entry.last_hit - diag_offset;
    let diff = subject_offset - last_hit;
    
    // NCBI: "do an extension, but only if we have not already extended this far"
    // NCBI reference: aa_ungapped.c:782
    if diff >= 0 {
        // Convert to raw positions (+1 for sentinel)
        let q_raw = (query_offset + 1) as usize;
        let s_raw = (subject_offset + 1) as usize;
        
        if q_raw + 2 >= q_aa.len() || s_raw + 2 >= s_aa.len() {
            continue;
        }
        
        // NCBI: s_BlastAaExtendOneHit
        // NCBI reference: aa_ungapped.c:787-791
        let x_dropoff = x_dropoff_per_context_ref[ctx_flat];
        let (hsp_q, hsp_qe, hsp_s, _hsp_se, score, _right_extend, s_last_off) =
            extend_hit_ungapped(
                q_aa,
                s_aa,
                q_raw,
                s_raw,
                0, // seed_score (not used in one-hit mode)
                x_dropoff,
            );
        
        let hsp_len = hsp_qe.saturating_sub(hsp_q);
        
        // NCBI: Update last_hit after extension
        // NCBI reference: aa_ungapped.c:799-800
        diag_entry.last_hit = (s_last_off as i32) - (wordsize - 1) + diag_offset;
        
        // NCBI: Check score threshold
        // NCBI reference: aa_ungapped.c:794-797
        let cutoff = cutoff_scores_ref[ctx_flat][s_f_idx];
        if score < cutoff {
            continue;
        }
        
        // NCBI: BlastSaveInitHsp equivalent
        // Store HSP with absolute coordinates
        let frame_base = 0i32;
        let hsp_q_absolute = frame_base + (hsp_q as i32) - 1;
        let hsp_qe_absolute = frame_base + ((hsp_q + hsp_len) as i32) - 1;
        
        init_hsps.push(InitHSP {
            q_start_absolute: hsp_q_absolute,
            q_end_absolute: hsp_qe_absolute,
            s_start: hsp_s as i32,
            s_end: (hsp_s + hsp_len) as i32,
            score,
            ctx_idx: ctx_flat,
            s_f_idx,
            q_idx,
            s_idx: s_idx as u32,
            q_frame: q_frame.frame,
            s_frame: s_frame.frame,
            q_orig_len: q_frame.orig_len,
            s_orig_len: s_len,
        });
    }
    // If diff < 0, skip (already extended past this position)
    continue;
}
```

---

## 10. 検証項目

### 機能テスト

1. **One-Hit Mode動作確認**:
   - `--window-size 0`で実行
   - Extensionが正しく実行されることを確認
   - HSP数がtwo-hit modeと異なることを確認

2. **Two-Hit Mode回帰テスト**:
   - デフォルト（`window_size=40`）で実行
   - 既存のテストが通ることを確認

3. **NCBI BLAST+比較**:
   - 同じクエリ/サブジェクトでNCBI BLAST+と比較
   - `-window_size 0`と`--window-size 0`で出力を比較
   - HSP数、スコア、座標が一致することを確認

### パフォーマンステスト

1. **One-Hit Modeパフォーマンス**:
   - 長い配列（600kb+）で実行
   - Extension数がtwo-hit modeより多いことを確認（期待される動作）

---

## 11. まとめ

### 必須改修項目

1. **One-Hit Mode実装** (Priority: HIGH)
   - `window_size == 0`の場合の処理を追加
   - `diff >= 0`チェックのみ（two-hitチェックなし）
   - NCBIコード参照を追加

### 確認済み（改修不要）

1. **Off-Diagonal Checking**: TBLASTXには不要（blastn特有）
2. **Extension関数**: 正しく実装済み
3. **Context Management**: 正しく実装済み
4. **Cutoff Score適用**: 正しく実装済み
5. **Two-Hit Mode**: 正しく実装済み

### 実装優先度

1. **Phase 1**: One-Hit Mode実装（必須）
2. **Phase 2**: テストと検証
3. **Phase 3**: ドキュメント更新

---

## 12. NCBIコード参照

### 主要参照ファイル

- `ncbi-blast/c++/src/algo/blast/core/hspfilter_culling.c` ⚠️ 最優先
  - Lines 52-60: `LinkedHSP`構造体定義
  - Lines 78-120: `s_DominateTest`関数（支配判定）
  - Lines 122-133: `s_FullPass`関数（merit更新）
  - Lines 135-163: `s_ProcessHSPList`関数（HSPリスト処理）
  - Lines 200-207: `CTreeNode`構造体定義（interval tree）
  - Lines 228-247: `s_CTreeNodeNew`関数（ノード作成）
  - Lines 332-370: `s_ProcessCTree`関数（tree処理）
  - Lines 430-470: `s_SaveHSP`関数（HSP追加）

- `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c`
  - Lines 200-236: `BlastAaWordFinder` (one-hit/two-hit切り替え)
  - Lines 712-812: `s_BlastAaWordFinder_OneHit` (one-hit mode実装)
  - Lines 439-619: `s_BlastAaWordFinder_TwoHit` (two-hit mode実装)
  - Lines 151-162: `s_BlastAaExtendOneHit` (one-hit extension)
  - Lines 184-197: `s_BlastAaExtendTwoHit` (two-hit extension)

### 実装時の注意事項

1. **NCBIコードを必ず参照すること**
2. **コメントにNCBIコードの行番号を記載すること**
3. **推測せず、NCBI実装を忠実に移植すること**
4. **テストでNCBI BLAST+と出力を比較すること**

