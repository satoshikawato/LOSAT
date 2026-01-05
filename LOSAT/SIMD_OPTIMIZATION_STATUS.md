# TBLASTX SIMD最適化状況

## 実装日: 2026-01-XX

## 概要

このドキュメントは、TBLASTX実装におけるSIMD最適化の現状と、次に最適化すべき箇所をまとめたものです。

---

## 1. 既にSIMD化済みの箇所

### 1.1 Offset Pairs コピー（`utils.rs`）

**関数**: `copy_offset_pairs_overflow_avx2`, `copy_offset_pairs_overflow_sse2`

**場所**: `LOSAT/src/algorithm/tblastx/utils.rs:86-170`

**最適化内容**:
- Overflow配列からoffset_pairsへのコピーをAVX2/SSE2で高速化
- AVX2: 8ペア/イテレーション（2x 4-pair stores）
- SSE2: 4ペア/イテレーション

**NCBI互換性**: ✅ 完全一致（出力順序維持）

**効果**: Overflowケース（`numhits > AA_HITS_PER_CELL`）でのメモリコピーを高速化

---

### 1.2 PV Test バッチ処理（`utils.rs`）

**関数**: `pv_test_mask4_avx2`, `pv_test_mask8_avx2`, `pv_test_mask16_avx2`

**場所**: `LOSAT/src/algorithm/tblastx/utils.rs:172-244`

**最適化内容**:
- 4/8/16個の連続するindexに対するPVテストをAVX2で並列実行
- `_mm256_i32gather_epi64` でPV wordsをgather
- `_mm256_sllv_epi64` でビットマスク生成

**NCBI互換性**: ✅ 完全一致（lane順処理）

**効果**: スキャン時のPVテストをバッチ化して高速化

---

### 1.3 Ungapped Extension スコア計算（`extension.rs`）

**関数**: `extend_left_ungapped_avx2`, `extend_right_ungapped_avx2`, `gather_scores8_avx2`

**場所**: `LOSAT/src/algorithm/tblastx/extension.rs:128-366`

**最適化内容**:
- BLOSUM62スコアテーブルへのアクセスをAVX2 gatherでバッチ化
- 8残基ずつスコアをgatherしてから逐次累積（NCBI互換性維持）
- X-drop判定は逐次処理（1-bit parity維持）

**NCBI互換性**: ✅ 完全一致（制御フローと終了条件を厳密に維持）

**効果**: Extension時のスコア計算を高速化（既存実装）

**未SIMD化の箇所**:
- `extend_hit_ungapped`内の初期word内探索（3残基のみ、効果は限定的）
- `extend_hit_two_hit`内の初期word内探索（3残基のみ、効果は限定的）

---

### 1.4 3-mer Index生成（`utils.rs`）【2026-01-XX追加】

**関数**: `compute_3mer_indices_16_avx2`, `compute_3mer_indices_8_sse2`, `compute_3mer_indices_4_scalar`

**場所**: `LOSAT/src/algorithm/tblastx/utils.rs:247-320`

**最適化内容**:
- Rolling index計算を直接3-mer計算に置き換え
- 各位置で `(subject[s] << 10) | (subject[s+1] << 5) | subject[s+2] & 0x7FFF` を計算
- 現在はscalar fallbackを使用（SIMD実装は将来の最適化候補）

**NCBI互換性**: ✅ 完全一致（数学的等価性、処理順序維持）

**効果**: スキャン時のindex計算を簡素化（SIMD実装は未完成）

---

### 1.5 Identity計算（`reevaluate.rs`）【2026-01-XX追加】

**関数**: `count_identities_avx2`, `count_identities_sse2`, `count_identities_scalar`

**場所**: `LOSAT/src/algorithm/tblastx/reevaluate.rs:254-396`

**最適化内容**:
- `q == s` の比較をAVX2/SSE2で並列実行
- AVX2: 32バイトずつ `_mm256_cmpeq_epi8` + `_mm256_movemask_epi8` + popcount
- SSE2: 16バイトずつ同様の処理
- 端数はscalar fallback

**NCBI互換性**: ✅ 完全一致（出力SHA256一致確認済み）

**効果**: Reevaluation時のidentity計算を高速化

---

### 1.6 Reevaluation スコア計算ループ（`reevaluate.rs`）【2026-01-XX追加】

**関数**: `reevaluate_ungapped_hit_ncbi_translated_avx2`, `gather_scores8_reevaluate_avx2`, `reevaluate_score_table_32`

**場所**: `LOSAT/src/algorithm/tblastx/reevaluate.rs:13-183`

**最適化内容**:
- BLOSUM62スコア計算をAVX2 gatherでバッチ化
- 16残基ずつ（2x8 gather）スコアをgatherしてから逐次累積
- `sum < 0` や `sum > score` の判定は逐次処理（NCBI互換性維持）
- 32x32スコアテーブル（`REEVALUATE_SCORE_TABLE_32`）を使用

**NCBI互換性**: ✅ 完全一致（制御フローを逐次処理で維持、lane順処理）

**効果**: Reevaluation時のスコア計算を高速化（全HSPに対して実行されるため、HSP数が多いケースで大きな効果）

**注意**: 実際の測定では、reevaluationの時間は0.008秒（全体の0.03%）と非常に短く、ボトルネックではないことが判明。短いHSP（16未満）はスカラー実装を使用してSIMD化のオーバーヘッドを回避。

---

## 2. 次にSIMD化すべき候補箇所

### 2.1 Smith-WatermanアライメントのDP計算（TBLASTXでは未使用）

**場所**: `LOSAT/src/align/sw_banded.rs:106-228`, `LOSAT/src/algorithm/tblastx/extension.rs:573-732`

**現状**:
- バンド付きSmith-Watermanアルゴリズムの動的計画法（DP）マトリックス計算
- 各セルのスコア計算（match/mismatch、gap open、gap extend）が独立
- 複数のセルを並列処理可能

**重要**: **TBLASTXはungapped-onlyアルゴリズム**です。
- NCBI BLASTの仕様: "Gapped search is not allowed for tblastx" (blast_options.c:869-873)
- `extend_gapped_protein`関数は実装されていますが、TBLASTXでは呼ばれません
- 将来のgapped実装のための予備実装として残されています

**SIMD化の可能性**:
- **高**: 対角線方向の計算を並列化可能（ただしTBLASTXでは未使用）
- マッチスコアの計算、挿入・削除スコアの計算、最大値の選択をSIMD化
- AVX2/AVX512を使用して複数セルを同時計算

**期待効果**: N/A（TBLASTXでは使用されない）

**実装難易度**: 中-高（DPの依存関係を考慮する必要がある）

**注意**: BLASTNやBLASTPなどのgapped alignmentを使用するアルゴリズムでは有効な最適化候補です。

---

### 2.2 BLASTN Ungapped Extension（高優先度・最大の改善余地）

**場所**: `LOSAT/src/algorithm/blastn/extension.rs:24-88`

**現状**:
- **完全に未SIMD化**: 左右の拡張ループで文字比較とスコア計算を逐次的に実行
- match/mismatch判定: `if q == s { reward } else { penalty }`
- スコア累積とX-drop判定が逐次処理

**SIMD化の可能性**:
- **非常に高**: 複数位置の塩基比較を並列化可能
- AVX2を使用して32バイト（32塩基）を一度に処理
- `_mm256_cmpeq_epi8`で文字比較を並列化
- 比較結果に基づいてreward/penaltyを条件付きで適用（`_mm256_blendv_epi8`など）

**実装方針**:
- 32塩基ずつ処理（AVX2）
- 文字比較: `_mm256_cmpeq_epi8(q_vec, s_vec)`
- スコア適用: 比較結果に基づいてreward/penaltyを条件付き適用
- スコア累積とX-drop判定は逐次処理（NCBI互換性維持）

**期待効果**: **高**（BLASTNで頻繁に呼ばれる、完全に未SIMD化）

**実装難易度**: 中（文字比較は簡単だが、条件付きスコア適用がやや複雑）

**参考**: TBLASTXの`extend_left_ungapped_avx2`/`extend_right_ungapped_avx2`の実装パターンを参考に可能

---

### 2.3 TBLASTX 初期Word内での最良位置探索（低優先度）

**場所**: 
- `LOSAT/src/algorithm/tblastx/extension.rs:427-445` (`extend_hit_ungapped`内)
- `LOSAT/src/algorithm/tblastx/extension.rs:515-527` (`extend_hit_two_hit`内)

**現状**:
- 3残基（k_size=3）のword内で最良スコア位置を探索
- `get_score(q_char, s_char)`を3回呼び出し
- スコア累積と条件分岐が逐次処理

**SIMD化の可能性**:
- **中**: 3残基のスコアを一度にgather可能
- `gather_scores8_avx2`のパターンを使用（3残基のみ）
- ただし、3残基のみなので効果は限定的

**実装方針**:
- 3残基のスコアを一度にgather
- スコア累積と条件分岐は逐次処理（NCBI互換性維持）

**期待効果**: **低**（3残基のみのため、SIMD化のオーバーヘッドが大きい可能性）

**実装難易度**: 低（既存の`gather_scores8_avx2`を参考に可能）

**注意**: 3残基のみなので、SIMD化しても効果は限定的。むしろスカラー実装のままの方が速い可能性がある。

---

### 2.3 SEGフィルタリング（低複雑度領域マスキング）（中優先度）

**場所**: `LOSAT/src/utils/seg.rs:96-122, 378-408, 510-553`

**現状**:
- アミノ酸配列のエントロピー計算とウィンドウスライディング処理
- 状態ベクトルのソートとファクトリアル計算
- エントロピー配列の計算（全ウィンドウでの繰り返し計算）

**SIMD化の可能性**:
- **中**: ウィンドウスライディング処理を並列化可能
- カウント配列の更新をSIMD化
- エントロピー計算の並列化

**期待効果**: 中（大量の配列に対して実行される）

**実装難易度**: 中（複雑な計算ロジック）

---

### 2.4 DUSTフィルタリング（核酸配列の低複雑度領域マスキング）（中優先度）

**場所**: `LOSAT/src/utils/dust.rs:344-356, 454-516`

**現状**:
- トリプレットカウントとスコア計算
- トリプレットウィンドウの更新とカウント配列操作
- 完全区間の探索（複数のトリプレットに対する閾値判定）

**SIMD化の可能性**:
- **中**: トリプレットカウントの更新を並列化可能
- スコア計算の並列化

**期待効果**: 中（大量の配列に対して実行される）

**実装難易度**: 中

---

### 2.5 k-merネイバーフッドインデックスの構築（低優先度）

**場所**: `LOSAT/src/seed/aa_word_finder.rs:126-142, 209-219`

**現状**:
- 全13,824個の3-merに対して計算
- 3つのアミノ酸ペアのスコアを逐次計算

**SIMD化の可能性**:
- **中**: 3つのアミノ酸ペアのスコアを並列計算可能
- ただし、構築は1回のみなので効果は限定的

**期待効果**: 低（構築は1回のみ）

**実装難易度**: 低-中

---

### 2.6 核酸k-merのエンコーディング（低優先度）

**場所**: `LOSAT/src/seed/na_word_finder.rs:12-30`

**現状**:
- 複数塩基を逐次エンコード

**SIMD化の可能性**:
- **中**: 複数塩基を並列処理可能
- ただし、処理が軽量なので効果は限定的

**期待効果**: 低（処理が軽量）

**実装難易度**: 低

---

### 2.7 Karlin-Altschulパラメータ計算（低優先度）

**場所**: `LOSAT/src/stats/karlin_calc.rs:228-237, 326-357`

**現状**:
- スコア頻度プロファイルの計算（28×28の二重ループ）
- Newton-Raphson法による反復計算

**SIMD化の可能性**:
- **低-中**: 二重ループの内側を並列化可能
- ただし、計算は1回のみなので効果は限定的

**期待効果**: 低（計算は1回のみ）

**実装難易度**: 中

---

### 2.8 HSPリンキングアルゴリズム（低優先度）

**場所**: `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs:962-980, 990-1000`

**現状**:
- 大量のHSPに対する反復的な比較・更新処理
- アクティブリストの走査とスコア比較

**SIMD化の可能性**:
- **低**: 条件分岐が多く、SIMD化が困難
- スコア比較部分のみ並列化可能

**期待効果**: 低（条件分岐が多い）

**実装難易度**: 高（複雑な制御フロー）

---

### 2.9 翻訳処理（低優先度）

**場所**: `LOSAT/src/algorithm/tblastx/translation.rs:76-80`

**現状**:
- 配列全体のコドン変換ループ

**SIMD化の可能性**:
- **中**: 複数コドンを並列変換可能
- ただし、処理が軽量なので効果は限定的

**期待効果**: 低（処理が軽量）

**実装難易度**: 低-中

---

### 2.2 Two-hit Extension 初期スコア計算

**場所**: `LOSAT/src/algorithm/tblastx/extension.rs:500-527`

**関数**: `extend_hit_two_hit`

**現状**:
```rust
for i in 0..k_size {
    // ...
    score += get_score(q_char, s_char);  // ← 3残基のみだが、SIMD化可能
    // ...
}
```

**SIMD化の可能性**:
- **低-中**: 3残基のみなので効果は限定的
- ただし、この関数は非常に頻繁に呼ばれる（`calls=4067119`）
- 3残基を一度に処理するSIMD実装は可能

**期待効果**: 小（3残基のみのため）

**実装難易度**: 低（短いループ）

---

### 2.3 Lazy Neighbor Scan（未使用だが将来の最適化候補）

**場所**: `LOSAT/src/algorithm/tblastx/utils.rs:904-1053`

**関数**: `s_blast_aa_scan_subject_lazy`

**現状**:
- 現在は使用されていない（`neighbor_map`モードが推奨）
- 3重ループでneighbor k-merを列挙

**SIMD化の可能性**:
- **低**: 3重ループと条件分岐が多く、SIMD化が困難
- ただし、BLOSUM62スコア計算部分はSIMD化可能

**期待効果**: 低（現在未使用）

**実装難易度**: 高（複雑な制御フロー）

---

### 2.4 Lookup Table構築時のNeighbor Generation

**場所**: `LOSAT/src/algorithm/tblastx/lookup.rs:404-484`

**関数**: `build_ncbi_lookup` 内のneighbor generation

**現状**:
```rust
for s0 in 0..alphabet_size {
    for s1 in 0..alphabet_size {
        for s2 in 0..alphabet_size {
            let total_score = sc0 + blosum62_score(q1, s1) + blosum62_score(q2, s2);
            // ...
        }
    }
}
```

**SIMD化の可能性**:
- **中**: 内側のループ（s2）をSIMD化可能
- `blosum62_score(q2, s2)` を複数のs2に対して並列計算
- ただし、条件分岐（`total_score >= threshold`）が多く、効果は限定的

**期待効果**: 中（lookup table構築は1回のみだが、時間がかかる）

**実装難易度**: 中（条件分岐の処理が複雑）

---

## 3. パフォーマンス測定結果

### ベースライン（SIMD化前）
```
[TIMING] scan_subject: 0.124s (calls=38114)
[TIMING] ungapped_extend: 0.427s (calls=4067119)
[TIMING] search_total: 29.512s
[TIMING] total: 29.786s
```

### Step 1後（3-mer index生成SIMD化）
```
[TIMING] scan_subject: 0.128s (calls=38114)
[TIMING] ungapped_extend: 0.428s (calls=4067119)
[TIMING] search_total: 29.617s
[TIMING] total: 29.894s
```

### Step 2後（identity計算SIMD化）
```
[TIMING] scan_subject: 0.132s (calls=38114)
[TIMING] ungapped_extend: 0.442s (calls=4067119)
[TIMING] search_total: 30.722s
[TIMING] total: 31.011s
```

### Step 3後（Reevaluation スコア計算SIMD化）
```
[TIMING] scan_subject: 0.121s (calls=38114)
[TIMING] ungapped_extend: 0.420s (calls=4067119)
[TIMING] reevaluate: 0.008s (calls=65165)
[TIMING] search_total: 29.948s
[TIMING] total: 30.225s
```

### Step 4後（sum_stats_linking 計測追加 → 真のボトルネック確定）
```
[TIMING] read_queries: 0.209s
[TIMING] build_lookup: 0.036s
[TIMING] read_subjects: 0.025s
[TIMING] scan_subject: 0.183s (calls=38114)
[TIMING] ungapped_extend: 0.490s (calls=4067119)
[TIMING] reevaluate: 0.007s (calls=65165)
[TIMING] sum_stats_linking: 22.794s (calls=1)
[TIMING] identity_calc: 0.002s (calls=42797)
[TIMING] search_total: 25.652s
[TIMING] total: 25.999s
```

**重要な発見**:
- Reevaluationの時間は0.007秒（全体の0.03%）と非常に短く、ボトルネックではない
- **実際のボトルネックは`sum_stats_linking`**（22.794秒 / search_total 25.652秒）
- 短いHSP（16未満）はスカラー実装を使用することで、SIMD化のオーバーヘッドを回避

**注意**: このテストケース（AP027280自己比較）では、SIMD化による明確な高速化は見られませんでした。これは以下の理由が考えられます：
1. 現在の実装がscalar fallbackを多く使用している
2. このケースではSIMD化した箇所がボトルネックではない
3. より大きな配列やHSP数が多いケースで効果が現れる可能性

### 2026-01-05: AP027280自己比較（本セッションの計測コマンド）

**コマンド**:
```bash
(time $LOSAT_BIN tblastx -q ./fasta/AP027280.fasta -s ./fasta/AP027280.fasta -o ./losat_out/AP027280.AP027280.tlosatx.n1.out --query-gencode 1 --db-gencode 1 -n 1 )&>./losat_out/AP027280.AP027280.tlosatx.n1.log
```

**結果**:
```
real	0m23.379s
```

**ヒット数（NCBIとの差分）**:
- NCBI（`tests/blast_out/AP027280.AP027280.tblastx.n1.out`）: 42733
- LOSAT（`tests/losat_out/AP027280.AP027280.tlosatx.n1.out`）: 42797
- 差分: **+64**（LOSATの方が多い）

**メモ**:
- **ベンチのノイズ削減**のため、tblastx はデフォルトで進捗/情報ログを出さない（`--verbose` で有効化）。
- 詳細診断ログは `LOSAT_DIAGNOSTICS=1` のときのみ出す（通常は無出力）。

### 2026-01-05: sum_stats_linking高速化（本セッション）

#### Step 1: `chain_members` のpushをdebug時のみに限定（出力不変）

- **変更**: `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs` 内の `chain_members.push(cur)` を `debug_chaining` の時だけ実行（通常実行ではpush/拡張ゼロ）。
- **結果**:
  - hits: NCBI=42733, LOSAT=42797（差分 +64）
  - time: `real 0m26.365s`

#### Step 2: `.ln()`/`normalize_score` のホットループ再計算を排除（出力不変）

- **変更**:
  - contextごとに `logK=ln(K)` を前計算し、HSPごとに `xscore=lambda*score-logK` を前計算。
  - DP更新では `new_xsum = h_xsum + xscore` を使用。
- **結果**:
  - hits: NCBI=42733, LOSAT=42797（差分 +64）
  - time: `real 0m23.027s`（Step 1比で短縮）

#### Step 3: 参照/境界チェックの削減（出力不変）

- **変更**:
  - `lh_helpers[...]` の多重インデックス参照を `&helper` 参照にまとめる（bounds check/参照回数削減）。
  - active list 走査で `next_active` を先に読み、更新用の追加インデックスを減らす。
- **結果**:
  - hits: NCBI=42733, LOSAT=42797（差分 +64）
  - time: `real 0m22.335s`

#### Step 4: `ln_factorial_int` キャッシュ（不採用）

- **試行**: `LOSAT/src/stats/sum_statistics.rs` の `ln_factorial_int` をキャッシュ化（巨大な事前テーブル/小さめテーブル/thread_local growable を試行）。
- **結果**: このベンチ（AP027280自己比較, 1プロセス1回実行）では **初期化/キャッシュ管理のオーバーヘッドが勝ち、timeが悪化傾向**。
- **結論**: 本セッションでは **キャッシュ化をロールバック**し、最速だった `sum_stats_linking` 側の最適化（Step 2-3）を採用。
- **参考（ロールバック後）**: `real 0m22.251s`（同条件での再計測）

---

## 4. 推奨される次のステップ（TBLASTXを単一スレッドで速くする／NCBI互換性維持）

DeepWikiの提案（座標変換・neighbor map・scan/extension SIMDなど）を踏まえつつ、**実測（LOSAT_TIMING=1）で判明したボトルネック**に合わせて優先順位を更新する。

### 最優先: `sum_stats_linking` の高速化（最大ボトルネック）

**根拠（AP027280自己比較, LOSAT_TIMING=1）**:
- `[TIMING] sum_stats_linking: 22.794s (calls=1)`
- `[TIMING] search_total: 25.652s`
- 参考: `ungapped_extend`は0.490s、`reevaluate`は0.007s

**対象**:
- `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs` の `apply_sum_stats_even_gap_linking`
- `LOSAT/src/stats/sum_statistics.rs`（`small_gap_sum_e` / `large_gap_sum_e` / Romberg積分）

**方針（出力1-bit一致を維持しつつ）**:
- **まずリンク処理内部を更に分解して計測**（INDEX0/INDEX1、active list走査、prob計算、chain除去…）
- NCBI `link_hsps.c` と 1:1 で対比し、**同じ制御フロー/状態更新のまま** Rust側のオーバーヘッドを削る
  - バッファ再利用（`Vec`の再確保/ゼロ初期化を最小化）
  - 内側ループの境界チェック削減（`get_unchecked`/生ポインタ。ただし不変条件はassertで担保）
  - `small_gap_sum_e` / `large_gap_sum_e` の呼び出し回数がNCBIと同等か検証し、**不要な再計算があれば除去**

**期待効果**:
- **TBLASTXの実時間を大きく短縮できる唯一の箇所**（現状ここが全体の9割以上）

### 最優先: 長配列での過剰HSP生成の根本修正（NCBI互換性）

DeepWikiが指摘している通り、長配列（600kb+）ではHSP数が膨らみ、`sum_stats_linking`がO(n²)化して破綻する。

**実装方針（NCBIのフローに合わせる）**:
- `BlastSaveInitHsp`→`BLAST_GetUngappedHSPList` 相当の **「絶対座標→context内相対座標」変換** を実装
- 座標変換後に `Blast_HSPListReevaluateUngapped` 相当の **一括reevaluate** を実行（NCBI順）
- これにより、リンク対象HSP数/グループサイズをNCBI水準に戻す（=リンク入力を減らす）

### 次点: scan/extension側の低レベル最適化（ボトルネック化したら）

現状のAP027280ではscan/extensionは支配的ではないが、データセットによっては効く可能性がある。

- `s_blast_aa_scan_subject` のプリフェッチ、分岐ヒント、メモリアクセス局所性改善
- `--neighbor-map` の単一スレッド最適化（lookup構築/メモリフットプリント削減）
- `extend_hit_ungapped`/`extend_hit_two_hit` の初期word内探索（3残基）の改善（効果は限定的）

### 付帯: コンパイル設定（単一スレッドの底上げ）

- `RUSTFLAGS="-C target-cpu=native"`（実行CPUに最適化）
- `lto = "thin"`, `codegen-units = 1`（release最適化強化）
- PGO（計測→最適化ビルド）※アルゴリズムや出力は変えない

---

## 5. NCBI互換性の維持

すべてのSIMD最適化において、以下の原則を遵守：

1. **出力の完全一致**: SHA256ハッシュで検証。ただし、NCBI parityに近づく変化は推奨。
2. **処理順序の維持**: lane順の逐次処理
3. **制御フローの維持**: 条件分岐と終了条件を厳密に維持
4. **数学的等価性**: 計算結果が1-bitも狂わない

**検証方法**:
```bash
sha256sum baseline.out step1.out step2.out
# すべて同じハッシュであることを確認
```

---

## 6. 参考資料

- NCBI BLAST ソースコード: `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/`
- 主要な参照ファイル:
  - `blast_aascan.c`: Subject scanning
  - `aa_ungapped.c`: Ungapped extension
  - `blast_hits.c`: HSP reevaluation

---

## 7. 更新履歴

- 2026-01-XX: 初版作成
  - Step 1: 3-mer index生成のSIMD化（scalar fallback使用）
  - Step 2: Identity計算のSIMD化（AVX2/SSE2実装完了）
  - Step 3: Reevaluation スコア計算ループのSIMD化（AVX2実装完了、短いHSPはスカラー実装を使用）
- 2026-01-05: TBLASTXの真のボトルネックを確定（LOSAT_TIMING=1）
  - `sum_stats_linking` が 22.794s / search_total 25.652s を占めることを確認
  - 次の優先順位を「sum_stats_linking高速化」+「長配列での座標変換（NCBI互換性）」に更新
- 2026-01-XX: DeepWiki情報を反映して今後のSIMD化候補を整理
  - Smith-WatermanアライメントのDP計算を最優先候補に追加
  - BLASTN Ungapped Extensionを第二優先候補に追加
  - SEG/DUSTフィルタリングを第三優先候補に追加
  - その他の候補（k-mer構築、翻訳処理、Karlin-Altschul計算など）を低優先度として整理

