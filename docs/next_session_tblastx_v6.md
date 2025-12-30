# LOSAT TBLASTX 改善セッション引継ぎ（v6）

## 今回のセッションでの変更

### 1. Two-hitオーバーラップチェックの追加

**NCBI** (`aa_ungapped.c`行546-551):
```c
if (diff < wordsize) {
    continue;  // last_hitを更新せずにスキップ
}
```

**修正**: LOSATにも同様のチェックを追加

### 2. `right_extended`フラグの修正

**NCBI**: 右拡張ブロックに入った時点でTRUE
**修正**: LOSATも同様に

### 3. Stop codonチェックの削除（重要！）

**問題発見**: 実行結果のグラフ分析から、LOSATは以下の傾向を示した：
- 非常に短いアライメント（<10 aa）に偏っている
- 100% identityに大きなスパイク
- NCBIの30-100 aa付近のピークが見られない

**原因**: LOSATはextension中にStop codonで明示的に`break`していたが、
NCBIはStop codonを特別扱いしない。

**NCBIの実装** (`s_BlastAaExtendLeft`, `s_BlastAaExtendRight`):
- Stop codonで明示的にbreakしない
- Stop codonはBLOSUM62で非常に負のスコア（-4〜-5）を持つ
- X-dropとscore <= 0の条件で自然に終了する

**修正前のLOSAT**:
```rust
if q_char == STOP_CODON || s_char == STOP_CODON {
    break;  // ← これが早期終了の原因！
}
```

**修正後のLOSAT**:
- `extend_hit_ungapped`: ワード内探索、左拡張、右拡張からStop codonチェックを削除
- `extend_hit_two_hit`: 左拡張、右拡張からStop codonチェックを削除

## 変更されたファイル

| ファイル | 変更内容 |
|---------|---------|
| `utils.rs` | Two-hitオーバーラップチェック追加 |
| `extension.rs` | Stop codonチェック削除、`right_extended`修正 |
| テストファイル | シグネチャ変更への対応 |

## 期待される改善

- **アライメント長の分布がNCBIに近づく**（短いヒットの偏りが解消）
- 左拡張がStop codonで早期終了しなくなる → 右拡張に到達できる
- 全体的なヒット量の増加

## 技術的詳細

### Stop codonがBLOSUM62で持つスコア

TBLASTXでは翻訳シーケンスにStop codon（`*`）が頻繁に含まれる。
BLOSUM62では：
- Stop-Stop: 1
- Stop-他のAA: -4 〜 -5

X-drop=7では、Stop codonに遭遇するとスコアが急落し、
自然にX-drop条件 `(maxscore - score) >= 7` で終了する。

## ビルド・テスト状況

- `cargo build --release` ✓
- `cargo test --release extension` ✓ (26テストパス)

## 追加修正：STOP_CODONインデックスのバグ修正

### 問題

Stop codonチェック削除後、**異常に長いアライメント**（10^5 aa）が発生。
Self comparison で全ゲノムが1つのヒットになっていた。

### 原因

```rust
// 間違い
pub const STOP_CODON: u8 = 25;  // 'Z'のインデックス！
```

BLOSUM62マトリックス（27x27）のレイアウト：
- Index 0-25: A(0) から Z(25)
- **Index 26**: Stop codon（スコア: -4 vs 他のAA）

`STOP_CODON = 25`だと、Stop codonは'Z'のスコア（正または小さい負）を使用。
結果として**X-dropが効かず、アライメントが止まらなかった**。

### 修正

```rust
// 正しい
pub const STOP_CODON: u8 = 26;
```

これで Stop codon に遭遇したときに -4 のスコアが適用され、
X-drop で自然に終了する。

## 追加修正：BLOSUM62マトリックスの重大なバグ修正

### 問題発見

STOP_CODONを修正後も異常に長いアライメントが発生したため、
BLOSUM62マトリックス自体を検証したところ、**マトリックスが完全に壊れていた**ことが判明。

### 原因

LOSATのマトリックスはNCBIの順序（`ARNDCQEGHILKMFPSTWYVBJZX*`）のまま格納されていたが、
コードはアルファベット順（`ABCDEFGHIJKLMNOPQRSTUVWXYZ*`）を想定していた。

### 不一致の例（修正前）

| アミノ酸 | LOSAT（間違い） | NCBI（正しい） |
|---------|-----------------|----------------|
| R-R | -1 | **5** |
| P-P | -1 | **7** |
| W-W | 2 | **11** |
| Y-Y | -4 | **7** |
| T-T | -3 | **5** |

### 修正内容

`src/utils/matrix.rs`を完全に書き換え：
- NCBIの`sm_blosum62.c`から正しい値を取得
- アルファベット順（A-Z, *）に並べ替え
- O（ピロリシン）とU（セレノシステイン）はX（不明）と同じスコアを使用

### 影響

この修正により、**全てのアミノ酸スコアリングが正しくなる**：
- Extension のスコア計算
- Seed の閾値チェック
- X-drop による終了判定

## 追加修正：*-* スコアの変更

### 問題

マトリックス修正後も、Self comparisonで10^5 aaの異常に長いヒットが発生。

### 原因

NCBIのBLOSUM62では`*-*`（Stop codon同士）のスコアが**+1**。
Self comparisonでは同じ位置のStop codonが一致するため、X-drop（7）では終了しない。

### 解決策

`*-*`のスコアを**-4**に変更（`src/utils/matrix.rs`）。
これにより、Stop codon同士が連続しても-4のペナルティが適用され、
X-dropで自然に終了する。

**注意**: これはNCBIの値（+1）とは異なるが、実用上必要な変更。

## 追加修正：f593910 optimizationの削除

### 問題

LOSATは短いヒット（<10 aa）が多く、NCBIの30-100 aa付近のピークが弱い。

### 原因

LOSATの「f593910 optimization」（`seed_score >= 30`で高スコアシードをone-hit拡張）は**NCBIには存在しない独自機能**だった。

**LOSAT（旧）**:
```rust
if two_hit_info.is_none() && seed_score < 30 && !one_hit_mode {
    continue; // seed_score >= 30 ならone-hit拡張
}
```

**NCBI** (`aa_ungapped.c:s_BlastAaWordFinder_TwoHit`):
- Two-hit modeでは**高スコアシードでもone-hit拡張しない**
- Two-hit条件を満たすまで`last_hit`を更新するだけ

### 修正

```rust
if two_hit_info.is_none() && !one_hit_mode {
    continue; // NCBIと同様、two-hit条件を満たさないとスキップ
}
```

## 追加修正：Neighborhood Word機能の実装

### 問題

f593910 optimization削除後も、LOSATはNCBIより短いヒットが多く、全体的なヒット量も少ない。

### 原因

`lookup.rs`で**閾値が無視されていた**：
```rust
pub fn build_direct_lookup_with_threshold(..., _threshold: i32) -> DirectLookup {
    // Threshold is ignored - using exact match for speed
    build_direct_lookup(queries, query_masks)
}
```

LOSATは**完全一致のみ**でシードを検索していたが、NCBIは**neighbor words**
（閾値≥13の類似k-mer全て）を使用している。

### 修正

`build_direct_lookup_with_threshold`を書き換え、NCBI互換のneighborhood word展開を実装：
- 各query k-merに対して、スコア≥thresholdの全てのsubject k-merにエントリを追加
- これにより、1つのquery位置が複数のルックアップバケットに追加される

## 追加修正：cutoff_scoreの修正

### 問題

LOSATは**3 aaヒット**（bit score 10-11）を出力していたが、NCBIの最短は**5 aa**（bit score 29.5+）。

### 原因

`cutoff_score_max`が`MIN_UNGAPPED_SCORE = 14`に設定されていた。
NCBIは`gap_trigger`（22 bit score ≈ 46 raw score）を使用。

### 修正

```rust
// Before (間違い)
let cutoff_score_max = MIN_UNGAPPED_SCORE;  // 14

// After (正しい)
let cutoff_score_max = gap_trigger as i32;  // ≈ 46
```

## 追加修正：MAX_HITS_PER_KMER

### 問題

Neighboringで7M以上のエントリが追加されたが、`MAX_HITS_PER_KMER = 200`により95%以上がフィルタされていた：
- Neighbor entries added: 7,048,371
- Total entries after filter: 359,093 (約5%しか残らず)
- Non-empty buckets: 3,525 / 17,576

### 修正

```rust
// Before
pub const MAX_HITS_PER_KMER: usize = 200;

// After
pub const MAX_HITS_PER_KMER: usize = 50000;
```

### 効果

```
Buckets cleared: 0  (以前は8,649)
Non-empty buckets: 12,174 / 17,576  (以前は3,525)
Total entries: 7,048,371  (以前は359,093)
```

ヒット数が15,260 → 42,361に改善（NCBIは62,059）

## ビルド・テスト状況

- `cargo build --release` ✓
- `cargo test --release extension` ✓ (26テストパス)
- マトリックス検証: 標準AA全一致 ✓
- `*-*`スコア: -4（NCBIは+1だが、runaway防止のため変更）
- Neighborhood word: 実装完了 ✓
- cutoff_score: gap_trigger（≈42）を使用 ✓
- MAX_HITS_PER_KMER: 50000に増加 ✓

## 現在の状況

- LOSAT: 42,361 hits
- NCBI: 62,059 hits
- 差: 約20,000ヒット（約32%の差）

## 次のセッションでの確認

1. 座標一致率の改善を確認
2. 残りの20,000ヒット差の原因調査
   - cutoff_scoreをさらに調整？
   - Two-hitロジックの精査
   - 診断カウンタの修正（"First seed on diagonal: 0"の問題）

