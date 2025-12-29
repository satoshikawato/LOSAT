# TBLASTX 今後の方針

## 優先度高

### 1. ヒット数差異の調査

現在のヒット数がf593910/NCBIより少ない問題を調査。

**現状の診断結果:**
- Extensions: ほぼ同じ（1,968,086 vs 1,892,458）
- cutoff_scoreフィルタ: 多く弾かれている（582,681 vs 466,503）

**原因候補:**
1. `extension.rs`のstop codonチェック
   - f593910: 4箇所
   - 現在: 9箇所
   - extensionが短くなり、スコアが低くなる可能性

2. extension内部ロジックの差異

**調査方法:**
```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT
git diff f593910 HEAD -- LOSAT/src/algorithm/tblastx/extension.rs | less
```

### 2. NCBIとの互換性向上

- neighboring words threshold (13) の正しい実装方法を再検討
- sum-statisticsリンキングの再実装（現在は無効化）

---

## 優先度中

### 3. 他のテストケースでの検証

- AP027132.NZ_CP006932.tblastx
- 小さいゲノム同士の比較

### 4. E-value計算の検証

- 20-29 bit score帯のヒットがNCBIと同じE-valueになるか確認

---

## 優先度低

### 5. neighboring words再実装

現在は完全一致のみ（f593910方式）で高速。
NCBIはthreshold=13でneighboring wordsを使用。

**以前の試み:**
- 構築時に全neighbor展開 → 遅すぎた
- スキャン時にneighborテーブル使用 → 3分かかった

**代替案:**
- 事前計算neighbor展開の高速化
- より効率的なデータ構造の検討

### 6. sum-statistics完全実装

`sum_stats_linking.rs`は現在パススルー。
NCBIスタイルのeven-gap linkingを正しく実装すると、E-value計算が改善される可能性。

参照:
- NCBI: `ncbi-blast/c++/src/algo/blast/core/link_hsps.c`
- LOSAT: `LOSAT/src/stats/sum_statistics.rs`（関数は実装済み）

