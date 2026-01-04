# TBLASTX Seeding/Extension段階でのNCBIとLOSATの違い

## 調査日: 2026-01-04

## 概要

TBLASTXのseeding/extension段階で、NCBI BLAST+とLOSATの間に以下の重要な違いが発見されました。これらの違いが、長い配列（600kb+）での過剰なHSP生成（338,859 vs 30,000-45,000）の原因である可能性が高いです。

---

## 1. Reevaluationのタイミングの違い

### NCBIの処理フロー

```
Extension → Cutoffチェック → BlastSaveInitHsp（保存）
  ↓
init_hitlistに保存（絶対座標）
  ↓
BLAST_GetUngappedHSPList（座標変換：絶対→context内相対）
  ↓
hsp_listに変換
  ↓
一括でBlast_HSPListReevaluateUngapped（全HSPをreevaluate）
```

**コード参照**:
- `aa_ungapped.c:588-591`: `BlastSaveInitHsp`で保存
- `blast_engine.c:560-563`: `BLAST_GetUngappedHSPList`で変換
- `blast_engine.c:1492-1497`: `Blast_HSPListReevaluateUngapped`で一括reevaluation

### LOSATの処理フロー

```
Extension → Cutoffチェック → Reevaluation（個別）→ 保存
```

**コード参照**:
- `utils.rs:991-1016`: Extension直後に個別にreevaluation

### 影響

タイミングの違い自体は問題ではありませんが、**座標変換のタイミング**が異なることが問題です。

---

## 2. 座標変換の欠如（最重要）

### NCBIの座標変換フロー

#### ステップ1: `BlastSaveInitHsp` (`aa_ungapped.c:588-591`)

```c
BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                 query_offset, subject_offset, hsp_len,
                 score);
```

- `hsp_q`, `hsp_s`: extension結果の**絶対座標**（`query->sequence`と`subject->sequence`に対する、sentinel含む）
- `query_offset`, `subject_offset`: 元のhitのoffset（concatenated buffer内の**絶対座標**）
- `ungapped_data->q_start = hsp_q`（**絶対座標**）
- `ungapped_data->s_start = hsp_s`（**絶対座標**）

#### ステップ2: `BLAST_GetUngappedHSPList` (`blast_gapalign.c:4756-4758`)

```c
context = s_GetUngappedHSPContext(query_info, init_hsp);
s_AdjustInitialHSPOffsets(init_hsp,
                          query_info->contexts[context].query_offset);
```

**`s_AdjustInitialHSPOffsets`** (`blast_gapalign.c:2384-2392`):
```c
static NCBI_INLINE void
s_AdjustInitialHSPOffsets(BlastInitHSP* init_hsp, Int4 query_start)
{
    init_hsp->offsets.qs_offsets.q_off -= query_start;
    if (init_hsp->ungapped_data) {
        init_hsp->ungapped_data->q_start -= query_start;  // ★座標変換
    }
    ASSERT(init_hsp->ungapped_data == NULL ||
           init_hsp->ungapped_data->q_start >= 0);
}
```

- **絶対座標 → context内相対座標に変換**
- `ungapped_data->q_start`がcontext開始位置からの相対座標になる

#### ステップ3: `Blast_HSPInit` (`blast_gapalign.c:4760-4767`)

```c
Blast_HSPInit(ungapped_data->q_start,  // context内相対座標
              ungapped_data->length+ungapped_data->q_start,
              ungapped_data->s_start,
              ungapped_data->length+ungapped_data->s_start,
              ...)
```

**`Blast_HSPInit`** (`blast_hits.c:170-171`):
```c
new_hsp->query.offset = query_start;  // ungapped_data->q_start（context内相対座標）
new_hsp->subject.offset = subject_start;
```

- `hsp->query.offset`は**context内相対座標**

#### ステップ4: `Blast_HSPListReevaluateUngapped` (`blast_hits.c:2696-2704`)

```c
context = hsp->context;
query_start = query_blk->sequence + query_info->contexts[context].query_offset;
query = query_start + hsp->query.offset;  // context内相対座標
```

- `hsp->query.offset`は**context内相対座標**（`s_AdjustInitialHSPOffsets`で変換済み）
- `query_start`はcontext開始位置のポインタ
- `query = query_start + hsp->query.offset`で正しい位置を取得

### LOSATの座標処理

#### Extension結果 (`utils.rs:922-935`)

```rust
let q_raw = (query_offset - ctx.frame_base) as usize;  // frame内相対座標
let query = &ctx.aa_seq;  // frame内の配列（sentinel含む）
let (hsp_q_u, hsp_qe_u, hsp_s_u, _hsp_se_u, score, right_extend, s_last_off_u) =
    extend_hit_two_hit(query, subject, ...);
let mut qs = hsp_q as usize;  // frame内相対座標（sentinel含む）
```

- `qs`は**frame内相対座標**（sentinel含む）のまま

#### Reevaluation (`utils.rs:1004`)

```rust
reevaluate_ungapped_hit_ncbi_translated(query, subject, qs, ss, len_u, cutoff)
```

- `qs`は**frame内相対座標**のまま使用
- **座標変換が行われていない**

### 問題点

1. **NCBI**: `BlastSaveInitHsp`で絶対座標保存 → `BLAST_GetUngappedHSPList`でcontext内相対座標に変換 → `Blast_HSPListReevaluateUngapped`で使用
2. **LOSAT**: frame内相対座標のままreevaluationに使用（**座標変換なし**）

### 座標の基準点の違い

#### NCBI (`Blast_HSPListReevaluateUngapped`)

```c
query_start = query_blk->sequence + query_info->contexts[context].query_offset;
query = query_start + hsp->query.offset;
```

- `query_start`: context開始位置のポインタ（concatenated buffer内）
- `hsp->query.offset`: context開始位置からの**相対座標**
- `query`: 正しい位置のポインタ

#### LOSAT (`reevaluate_ungapped_hit_ncbi_translated`)

```rust
query = &ctx.aa_seq;  // frame内の配列（sentinel含む）
qs = hsp_q as usize;  // frame内相対座標（sentinel含む）
```

- `query`: frame内の配列（sentinel含む）
- `qs`: frame内相対座標（sentinel含む）
- **座標の基準点が異なる**

### 影響

座標の基準点が異なるため、reevaluation時の配列アクセス位置が異なり、結果が変わる可能性があります。これがHSP生成数の差の原因である可能性が高いです。

---

## 3. 座標変換の実装箇所の確認

### NCBIの実装箇所

1. **`BLAST_GetUngappedHSPList`** (`blast_gapalign.c:4719-4775`)
   - ungapped-onlyの場合に呼ばれる
   - `blast_engine.c:560-563`から呼び出される

2. **`s_AdjustInitialHSPOffsets`** (`blast_gapalign.c:2384-2392`)
   - 座標変換を実行
   - `init_hsp->ungapped_data->q_start -= query_start`

3. **`s_GetUngappedHSPContext`** (`blast_gapalign.c:2370-2375`)
   - contextを取得
   - `BSearchContextInfo(init_hsp->offsets.qs_offsets.q_off, query_info)`

### LOSATの実装状況

- **座標変換が実装されていない**
- `utils.rs:999-1016`でreevaluationを実行するが、座標変換なし

---

## 4. 処理フローの完全な比較

### NCBI (ungapped-only)

```
1. aa_ungapped.c: s_BlastAaWordFinder_TwoHit
   └─ Extension実行
   └─ aa_ungapped.c:588-591: BlastSaveInitHsp
      └─ ungapped_data->q_start = hsp_q (絶対座標)
      └─ init_hitlistに保存

2. blast_engine.c:560-563: BLAST_GetUngappedHSPList
   └─ blast_gapalign.c:4756-4758: 座標変換
      └─ s_AdjustInitialHSPOffsets
         └─ ungapped_data->q_start -= query_offset (context内相対座標に変換)
   └─ blast_gapalign.c:4760-4767: Blast_HSPInit
      └─ hsp->query.offset = ungapped_data->q_start (context内相対座標)
   └─ hsp_listに保存

3. blast_engine.c:1492-1497: Blast_HSPListReevaluateUngapped
   └─ blast_hits.c:2696-2704: 座標計算
      └─ query_start = query_blk->sequence + query_info->contexts[context].query_offset
      └─ query = query_start + hsp->query.offset (context内相対座標)
   └─ blast_hits.c:2702-2705: Blast_HSPReevaluateWithAmbiguitiesUngapped
      └─ reevaluation実行
```

### LOSAT

```
1. utils.rs:922-930: extend_hit_two_hit
   └─ Extension実行
   └─ hsp_q, hsp_s取得（frame内相対座標）

2. utils.rs:991-1016: Cutoffチェック + Reevaluation
   └─ utils.rs:1004: reevaluate_ungapped_hit_ncbi_translated
      └─ qs, ss使用（frame内相対座標、座標変換なし）
   └─ 保存
```

---

## 5. 修正が必要な箇所

### 修正案

LOSATでも`BLAST_GetUngappedHSPList`相当の座標変換を実装する必要があります。

1. **座標変換の実装**
   - `BLAST_GetUngappedHSPList`相当の関数を実装
   - `s_AdjustInitialHSPOffsets`相当の座標変換を実行
   - context内相対座標に変換

2. **Reevaluationのタイミング変更**
   - Extension直後のreevaluationを削除
   - 座標変換後に一括でreevaluationを実行

3. **座標の基準点の統一**
   - NCBIと同様に、context開始位置からの相対座標を使用
   - `query_start + hsp->query.offset`の形式でアクセス

---

## 6. コード参照

### NCBI

- `aa_ungapped.c:588-591`: `BlastSaveInitHsp`の呼び出し
- `blast_extend.c:360-375`: `BlastSaveInitHsp`の実装
- `blast_engine.c:560-563`: `BLAST_GetUngappedHSPList`の呼び出し
- `blast_gapalign.c:4719-4775`: `BLAST_GetUngappedHSPList`の実装
- `blast_gapalign.c:2384-2392`: `s_AdjustInitialHSPOffsets`の実装
- `blast_gapalign.c:2370-2375`: `s_GetUngappedHSPContext`の実装
- `blast_hits.c:150-189`: `Blast_HSPInit`の実装
- `blast_engine.c:1492-1497`: `Blast_HSPListReevaluateUngapped`の呼び出し
- `blast_hits.c:2609-2737`: `Blast_HSPListReevaluateUngapped`の実装
- `blast_hits.c:2696-2704`: 座標計算部分
- `blast_hits.c:675-733`: `Blast_HSPReevaluateWithAmbiguitiesUngapped`の実装

### LOSAT

- `utils.rs:918-930`: Extension実行
- `utils.rs:991-1016`: Cutoffチェック + Reevaluation
- `utils.rs:1004`: `reevaluate_ungapped_hit_ncbi_translated`の呼び出し
- `reevaluate.rs:80-145`: `reevaluate_ungapped_hit_ncbi_translated`の実装

---

## 7. 結論

**最も重要な違い**: 座標変換の欠如

- NCBIは`BLAST_GetUngappedHSPList`内で絶対座標からcontext内相対座標に変換
- LOSATは座標変換を行わず、frame内相対座標のままreevaluationに使用
- これにより、reevaluation時の配列アクセス位置が異なり、結果が変わる可能性がある

この違いが、長い配列での過剰なHSP生成（338,859 vs 30,000-45,000）の根本原因である可能性が高いです。

---

## 8. 次のステップ

1. `BLAST_GetUngappedHSPList`相当の関数を実装
2. 座標変換ロジックを追加
3. Reevaluationのタイミングを変更（Extension直後 → 座標変換後）
4. テストでHSP生成数の差が解消されるか確認

