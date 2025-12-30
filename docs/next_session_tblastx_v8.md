# LOSAT TBLASTX 改善セッション引継ぎ（v8）

## 今回のセッションで発見・修正した問題

### 1. X_DROP_UNGAPPED の単位問題（重要・修正済み）

**問題発見:**
NCBIの`BLAST_UNGAPPED_X_DROPOFF_PROT = 7`は**bit score**であり、実行時に**raw score**に変換される。

**NCBIの変換式** (`blast_parameters.c:220-221`):
```c
x_dropoff_raw = x_dropoff_bits * ln(2) / lambda
             = 7 * 0.693 / 0.3176
             ≈ 15
```

**修正前のLOSAT:**
```rust
pub const X_DROP_UNGAPPED: i32 = 7;  // raw scoreとして使用していた（間違い！）
```

**修正後のLOSAT:**
```rust
pub const X_DROP_UNGAPPED_BITS: f64 = 7.0;
pub const X_DROP_UNGAPPED: i32 = 15;  // NCBIと同じ変換済み値
```

**効果:**
- 低identity比較: 6,488 → 8,020 hits（約24%増加）
- しかしまだNCBI（14,877）の約54%

## 残存する問題

### Self-comparisonでのrunaway extension

**症状:**
- LOSAT: 最長68,354 aa（100% identity）、120個の1000+ hits
- NCBI: 最長2,260 aa、6個の1000+ hits

**原因:**
- 100%一致の場合、各アミノ酸で+4〜+11のスコア
- Stop codon (`*-*` = -4) では累積スコアが15以上下がらない
- X-drop条件 `(maxscore - score) >= 15` を満たせない

**NCBIの謎:**
- NCBIでも`*-*` = +1（我々は-4に変更した）
- なぜNCBIは2,260 aaで止まるのか不明

## NCBIとの徹底比較結果

### 確認済み（同じ動作）

| 項目 | 値 | 確認方法 |
|------|-----|---------|
| Two-hit window | 40 | `blast_options.h:57` |
| cutoff_score | gap_trigger ≈ 42 raw | `blast_parameters.c:369` |
| 対角線座標 | `(q - s) & mask` | `aa_ungapped.c:516` |
| ワード内スコア | 累積最大位置 | `aa_ungapped.c:1108-1119` |
| 左拡張開始 | score = 0 | `aa_ungapped.c:893` |
| 右拡張開始 | score = left_score | `aa_ungapped.c:838` |
| 最終スコア | MAX(left, right) | `aa_ungapped.c:1157` |
| 拡張後の更新 | 常に実行 | `aa_ungapped.c:596-606` |

### NCBIコードの重要な参照箇所

```
aa_ungapped.c:
  - s_BlastAaWordFinder_TwoHit: 439-619 (メインループ)
  - s_BlastAaExtendTwoHit: 1088-1158 (拡張関数)
  - s_BlastAaExtendLeft: 886-921
  - s_BlastAaExtendRight: 831-866

blast_parameters.c:
  - x_dropoff変換: 219-221
  - cutoff_score計算: 348-374

blast_aalookup.c:
  - neighboring words: s_AddWordHitsCore (546-606)
```

## 現在のヒット数比較

### 低identity (AP027131 vs AP027133)

| バージョン | ヒット数 | NCBIとの比率 |
|-----------|---------|-------------|
| NCBI BLAST | 14,877 | 100% |
| LOSAT (X_DROP=7) | 6,488 | 43.6% |
| LOSAT (X_DROP=15) | 8,020 | 53.9% |

### Self-comparison (NZ_CP006932)

| バージョン | ヒット数 | 1000+ hits | 最長 |
|-----------|---------|-----------|------|
| NCBI BLAST | 62,059 | 6 | 2,260 aa |
| LOSAT (X_DROP=15) | 37,091 | 120 | 68,354 aa |

## 次のセッションでの優先作業

### 優先度1: NCBIのself-comparison動作の解明

NCBIで2,260 aaで止まる理由を特定：
1. `score <= 0`条件がどこかで満たされている？
2. 何か別の終了条件がある？
3. データに特殊文字がある？

### 優先度2: Runaway extension対策

考えられる解決策：
1. `*-*`スコアをさらに負にする（例: -8）
2. 拡張の最大長を制限する
3. 連続正スコアのカウントで終了

### 優先度3: ヒット数の差

残り約6,800ヒットの差の原因：
1. 拡張スコアの詳細比較
2. Neighboring wordsの確認
3. Two-hit条件の微細な違い

## 変更したファイル

| ファイル | 変更内容 |
|---------|---------|
| `constants.rs` | `X_DROP_UNGAPPED = 15` に変更 |
| `next_session_tblastx_v8.md` | このドキュメント |

## ビルド状態

```bash
cargo build --release  # ✓ 成功
cargo test --release   # ✓ 成功（一部DNA関連テストは除く）
```
