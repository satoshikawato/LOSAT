# 次セッションプロンプト: sum_stats_linking のパフォーマンス改善

## 前提

- **NCBIの実装がground truth**。出力に影響する"簡略化"は禁止。
- 必要なら構造を壊してでもNCBIに合わせる。
- 調査に逃げるな。コードを直せ。

## 前回セッションで修正済み

### 1. best 選択ループの修正
`sum_stats_linking.rs` lines 287-302 で `linked_to != -1000` チェックを追加。
NCBI は linked list で処理済み HSP が物理的に削除されるが、LOSAT は array なので明示的な除外が必要。

### 2. 抽出ループの安全チェック
lines 649-686 で処理済み HSP への遷移時に break を追加。無限ループと二重カウントを防止。

### 3. can_skip の安全化
lines 489-504 で singleton (`prev_link.link[1].is_none()`) のみ skip を許可。
処理済み HSP を含むチェーンの再利用バグを修正。

## 現在の問題: パフォーマンス

### 症状

```
NZ_CP006932 self-comparison (neighbor-map mode):
- 232,359 raw ungapped hits
- 4 groups (q_strand, s_strand)
- 各グループ 50-70K HSPs
- 60秒のタイムアウトでも完了しない
```

### 根本原因

1. **O(n²) の inner loop**: INDEX 0/1 ループが各 HSP に対して全先行 HSP を走査
2. **use_current_max が効かない**: dense な入力では path_changed が常に true
3. **can_skip が効かない**: 安全版は singleton のみ許可なので、ほとんどスキップできない
4. **多数のイテレーション**: 各イテレーションで 1-10 HSPs しか抽出しない

### NCBI との構造的差異

| 項目 | NCBI | LOSAT |
|------|------|-------|
| データ構造 | linked list | array |
| 処理済み HSP | リストから削除 (O(1)) | マークのみ (フィルタに O(n)) |
| can_skip | なし (prev_link.sum-1 をヒント) | singleton のみ |
| inner loop | linked list 走査 | array 走査 |

## 改善案

### 案1: linked list への移行（大規模変更）（推奨）

NCBI と同じ linked list 構造に変更:
- 処理済み HSP を O(1) で削除
- max 探索が自動的にアクティブ HSP のみを走査
- can_skip が不要になる

## テストコマンド

```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT
cargo build --release
timeout 300 ./target/release/LOSAT tblastx \
  -q tests/fasta/NZ_CP006932.fasta \
  -s tests/fasta/NZ_CP006932.fasta \
  --db-gencode 4 --query-gencode 4 --neighbor-map \
  -o /tmp/losat_output.out 2>&1 | head -50
```


### 案2: NCBI 準拠の can_skip（推奨）

NCBI は can_skip しない。代わりに `h_sum = prev_link.sum - 1` をヒントとして inner loop を高速化:

```rust
// NCBI lines 812-824: prev_link をヒントとして使用（スキップではない）
if !first_pass {
    if let Some(pl) = prev_link {
        if hsp_links[pl].linked_to >= 0 {
            h_sum = hsp_links[pl].sum[1] - 1; // 開始点を設定
        }
    }
}
// inner loop は必ず実行（next_larger jump で高速化）
```

**現状**: この実装は既に入っている。問題は inner loop 自体のコスト。



### 案3: active_hsps の再利用

現在、`!use_current_max` のたびに `active_hsps` を再構築:
```rust
let active_hsps: Vec<usize> = (0..n)
    .filter(|&i| hsp_links[i].linked_to != -1000)
    .collect();
```

**改善**: `active_hsps` を `Vec` ではなく `BitSet` で管理し、抽出時に O(1) で削除。

### 案4: 早期終了

E-value が閾値（例: 10.0）を超えるチェーンしか残っていない場合、早期終了:

```rust
if prob[0] > 10.0 && prob[1] > 10.0 {
    // 残りの HSP はフィルタで落ちるので、処理を打ち切る
    for i in 0..n {
        if hsp_links[i].linked_to >= 0 {
            hsp_links[i].linked_to = -1000;
        }
    }
    break;
}
```

**注意**: 出力同等性を維持するため、打ち切り前に singleton E-value を設定する必要あり。


## 成功基準

| 指標 | 現状 | 目標 |
|------|------|------|
| 実行時間 (neighbor-map) | タイムアウト (>300秒) | < 60秒 |
| Final hits | 未完了 | 62K (±10%) |
| 出力同等性 | - | NCBI と一致 |

## 関連ファイル

- `src/algorithm/tblastx/sum_stats_linking.rs` — 主要修正対象
- `docs/TBLASTX_NCBI_PARITY_DEVLOG.md` — 進捗記録

## NCBI リファレンス

- `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c` — メインアルゴリズム
- lines 607-625: max 探索
- lines 659-686: lh_helper 構築
- lines 812-824: prev_link ヒント
- lines 827-861: INDEX 1 inner loop (next_larger jump)

## 優先順位
1. **案1 (linked list)** — 大規模変更
2. **案2 (早期終了)** — 最小変更で効果大

3. **案3 (BitSet)** — 中程度の変更で O(n) フィルタを削減
4. **案4 (can_skip 検証)** — 現状で正しく動作しているか確認


