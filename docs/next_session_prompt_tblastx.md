# 次のセッション用プロンプト

以下のプロンプトをコピーして次のセッションで使用してください：

---

## プロンプト

LOSATのTBLASTX実装がNCBI BLASTと比較してヒット数が大幅に少ない問題の調査を継続してください。

### 背景
- LOSAT: 1,595ヒット
- NCBI BLAST: 62,053ヒット（40倍の差）
- 主な問題: Two-hitフィルタが97%のシードを除外
- 低スコアヒット（22-29ビットスコア）: LOSAT 0個 vs NCBI 27,985個

### 前セッションで実装済み
1. 非ギャップKarlin-Altschulパラメータ（λ=0.3176, K=0.134）
2. Two-hitフラグメカニズム（diag_extended）
3. 近傍ワード生成（閾値11）
4. K-merエンコーディング修正（base 26）
5. One-hitモードサポート（--window-size 0）

### 優先タスク
1. ルックアップテーブルにエントリ数診断を追加して、LOSATとNCBI BLASTのエントリ数を比較
2. NCBI BLASTの特定の低スコアヒット（例: 座標450063-449974対1135-1224、bit score 24.9）をLOSATでトレースし、なぜ見つからないか調査
3. 対角線ごとのk-mer一致数を出力し、off-diagonalでの一致パターンを分析

### 重要ファイル
- LOSAT: `LOSAT/src/algorithm/tblastx/utils.rs`, `lookup.rs`, `extension.rs`
- NCBI BLAST参照: `C:\Users\kawato\Documents\GitHub\ncbi-blast\c++\src\algo\blast\core\aa_ungapped.c`
- 引継ぎ文書: `docs/session_handoff_tblastx_comparison.md`

### テストコマンド
```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT
cargo build --release
LOSAT_DIAGNOSTICS=1 ./target/release/losat tblastx --query tests/fasta/NZ_CP006932.fasta --subject tests/fasta/NZ_CP006932.fasta --query-gencode 4 --db-gencode 4
```

### 目標
NCBIのアルゴリズムを正確に再現し、ヒット数の乖離を解消すること。

---

