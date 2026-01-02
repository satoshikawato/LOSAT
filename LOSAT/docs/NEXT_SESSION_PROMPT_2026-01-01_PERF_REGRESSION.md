# 次セッションプロンプト: sum_stats_linking パフォーマンス低下の調査と修正

## 前提

- **NCBIの実装がground truth**。出力に影響する"簡略化"は禁止。
- 必要なら構造を壊してでもNCBIに合わせる。
- **テストは全パッチ完了後に一度だけ**。無駄なテストはしない。

## 前回セッションで実施した変更

### 1. 侵入型双方向リンクリスト
`HspLink` に `next_active`, `prev_active` を追加し、O(1) でのアクティブ HSP 削除を実現。

### 2. `LhHelper` 構造体の拡張
- `hsp_idx`: 対応する HSP インデックス（NCBI の `ptr` フィールド相当）
- `maxsum1`: フレーム符号ごとの running max（NCBI lines 669-671, 874）

### 3. `next_larger` ジャンプ修正
デクリメント**前**の位置から読むように修正（NCBI lines 831-843）:
```rust
let current_idx = j_lh_idx;
let next_larger = lh_helpers[current_idx].next_larger;
j_lh_idx -= 1;
if b0 {
    j_lh_idx = next_larger;
}
```

### 4. `ordering_method` 追加
`HspLink.ordering_method` を追加し、チェーン抽出時に設定（NCBI line 973）。

### 5. デバッグ出力ゲート
全 `eprintln!("[DEBUG...]")` を `LOSAT_DIAGNOSTICS` でゲート。

## 現在の問題

**パフォーマンスが低下した**。以前は 33秒 で完了していたが、修正後は 90秒+ でもタイムアウト。

## 調査すべき箇所

### 1. `maxsum1` のオーバーヘッド
NCBI line 850:
```c
if(0) if(H2_helper->maxsum1<=H_hsp_sum)break;
```
`if(0)` で無効化されている。計算だけして使わないのはオーバーヘッド。**削除を検討**。

### 2. `next_larger` ジャンプが機能しているか
`sum_stats_linking.rs`:
- 初期構築: lines 418-426
- INDEX 1 ループ後更新: lines 622-629
- 内部ループでのジャンプ: lines 570-586

ジャンプが正しく機能していれば O(n) に近づくはず。機能していなければ O(n²)。

### 3. `use_current_max` 最適化が効いているか
`sum_stats_linking.rs` lines 299-375:
- `path_changed` の設定
- チェーン検証ループ

dense な入力では `path_changed` が常に true で `use_current_max = false` になる可能性。

### 4. 構造体サイズ増加
`LhHelper` と `HspLink` に新フィールドを追加したことでキャッシュ効率が低下した可能性。

## 修正方針

1. **`maxsum1` 計算を削除** - NCBI で使われていない
2. **`next_larger` ジャンプのデバッグ** - 本当にジャンプしているか確認
3. **プロファイリング** - どこで時間がかかっているか特定

## テストコマンド

```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT
cargo build --release
/usr/bin/time -p timeout 120 ./target/release/LOSAT tblastx \
  -q tests/fasta/NZ_CP006932.fasta \
  -s tests/fasta/NZ_CP006932.fasta \
  --db-gencode 4 --query-gencode 4 --neighbor-map --seg \
  -o /tmp/losat_output.out 2>&1
```

## 成功基準

| 指標 | 目標 |
|------|------|
| 実行時間 | < 60秒 |
| ヒット数 | ~62,000 (NCBI ±10%) |

## NCBI 参照ファイル

- `ncbi-blast/c++/src/algo/blast/core/link_hsps.c`
  - `s_BlastEvenGapLinkHSPs()`: lines 400-1000
  - INDEX 1 inner loop: lines 827-861
  - `next_larger` computation: lines 675-684
  - `maxsum1` (disabled): line 850

## 関連ファイル

- `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs`
- `LOSAT/docs/TBLASTX_NCBI_PARITY_DEVLOG.md`

