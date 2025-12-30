# セッション11 プロンプト

## コンテキスト
LOSATのTBLASTX実装をNCBI BLAST+と一致させる作業を継続中。

## 現在の最重要課題
**sum_stats_linkingが遅すぎる（2分以上）**

### 原因
- strand-basedグループ分け（4グループ）で1グループ約11万HSP
- O(N²)アルゴリズムで12兆回の比較
- NCBIの`lh_helper` + `next_larger`最適化が不完全

### 現状
- next_larger最適化を実装したが、whileループ内で毎回再計算しているため効果なし
- NCBIは初回のみ計算し、`changed`フラグで差分更新

## 作業指示

### 必須タスク
1. NCBIの`link_hsps.c`を読み、以下を完全実装せよ：
   - `changed`フラグによる差分更新（再計算スキップ）
   - `next_larger`の効率的な更新
   - `use_current_max`による最大値キャッシュ

2. 参照ファイル：
   - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c`
   - 特に589-900行のwhileループ

### 目標
- 処理時間: 10秒以内（現在2分以上）
- ヒット数: 62,053に近づける（現在9,077）

### 禁止事項
- 「実装が複雑すぎるので簡略化」は死刑
- NCBIの実装を忠実に再現せよ

## 参照ドキュメント
`@docs/next_session_tblastx_v11.md`

