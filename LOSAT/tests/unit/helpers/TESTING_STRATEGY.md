# NCBI BLASTテスト戦略（実装上の落とし穴と対策）

## 概要

このドキュメントは、NCBI BLASTとの整合性テストを実装する際に遭遇する可能性のある「実装上の落とし穴」と、その対策をまとめたものです。

## 4つの主要な落とし穴

### 1. 実効長（Effective Length）計算の複雑さ ⚠️ 最重要

**問題**: NCBI BLASTの実効長計算は、単純な `(QueryLen - K) * (DbLen - K)` ではありません。

**詳細**:
- HSPの期待長（Expected HSP length）を統計パラメータ（K, H）から逆算して減算
- 短い配列の場合、期待長自体が変動するロジックが含まれる

**✅ 実装状況**: 
LOSATは既にNCBI BLASTの`BLAST_ComputeLengthAdjustment`関数をRustに移植しています。
`src/stats/length_adjustment.rs`の`compute_length_adjustment_ncbi`関数が該当します。

**対策**:
- ✅ **実装の検証**: NCBI BLASTの実装と1行ずつ比較し、完全一致を確認
- ✅ **小さな許容誤差**: 浮動小数点精度の違いを考慮し、**0.1%の相対誤差**を許容
- ⚠️ より大きな差（>1%）が見られる場合は、実装のバグの可能性があるため修正が必要

**実装例**:
```rust
// ❌ 悪い例: 完全一致を要求
assert_eq!(losat_space, ncbi_space);

// ✅ 良い例: 1%の許容誤差
let relative_diff = (losat_space - ncbi_space).abs() / ncbi_space;
assert!(relative_diff <= 0.01, "Effective space diff: {:.2}%", relative_diff * 100.0);
```

### 2. X-drop（拡張停止）のバタフライ効果

**問題**: アライメントの「終わり」を決めるX-dropアルゴリズムの挙動が、NCBIとLOSATで1塩基でもズレると、アライメント長が変わり、Raw Scoreが変わり、結果としてBit Scoreもズレます。

**詳細**:
- 特に「ギャップの直後でスコアが落ちた時」の判定順序などが影響
- 動的に計算したHSPを使用すると、X-dropの微細な違いがスコアに影響する

**対策**:
- ✅ **テストの分離**: 
  - 「アライメント拡張ロジック（どこまで伸びるか）」のテスト
  - 「スコア・E-value計算ロジック（伸びた結果をどう評価するか）」のテスト
  - これらを**完全に分離**
- ✅ **固定アライメントでのテスト**: 
  - E-valueのテストには、動的に計算したHSPではなく、**開始位置・終了位置・ギャップ位置をハードコードした固定のHSP**を入力として与える
  - 計算式の正しさだけを検証

**実装例**:
```rust
// ❌ 悪い例: 動的に計算したHSPを使用
let hsp = extend_hit_ungapped(...);  // X-dropの影響を受ける
let (bit_score, e_value) = calculate_evalue(hsp.score, ...);

// ✅ 良い例: 固定アライメントを使用
let fixed_alignment = FixedAlignment {
    q_start: 100,
    q_end: 200,
    s_start: 100,
    s_end: 200,
    raw_score: 100,  // ハードコード
};
let (bit_score, e_value) = calculate_evalue(fixed_alignment.raw_score, ...);
```

### 3. デフォルトのフィルタリング（DUST/SEG）の罠

**問題**: NCBI BLASTはデフォルトで低複雑度領域のフィルタリング（DUST/SEG）がONになっています。

**詳細**:
- 抽出スクリプトが拾う出力結果は「フィルタリング後の配列」に対するスコア
- LOSATのユニットテストで「生配列」をそのまま入力すると、フィルタリングされずにフルマッチしてしまい、スコアが爆発的に高くなる（あるいはその逆）

**対策**:
- ✅ **正解データ作成時にフィルタOFF**: 
  - BLASTN: `-dust no` を明示的に追加
  - TBLASTX: `-seg no` を明示的に追加
- ✅ またはMaskingのテスト: LOSAT側にもMaskingロジックが実装されているなら、入力配列を「小文字（soft masking）」にしてテストする

**実行例**:
```bash
# ❌ 悪い例: デフォルト（フィルタリングON）
blastn -query test.fasta -subject db.fasta -outfmt 7 > output.txt

# ✅ 良い例: フィルタリングOFF
blastn -query test.fasta -subject db.fasta -dust no -outfmt 7 > output.txt
```

### 4. NCBI BLASTのバージョン差異

**問題**: NCBI BLASTはバージョンによって、内部のデフォルトパラメータ（Gap Costのペナルティ等）や統計補正の微調整がサイレントに行われることがあります。

**詳細**:
- 手元の `blastn (v2.16.0)` で作った正解データが、サーバーの `blastn (v2.9.0)` のロジックと合わない
- デフォルトパラメータがバージョン間で変更される可能性がある

**対策**:
- ✅ **バージョン固定**: README.mdに「リファレンスデータの生成には NCBI BLAST+ 2.X.X を使用すること」と明記
- ✅ **パラメータの明示**: 正解データ生成時に、デフォルト任せにせず `-gapopen`, `-gapextend`, `-reward`, `-penalty` を全てコマンドライン引数で明示的に指定して固定

**実行例**:
```bash
# ❌ 悪い例: デフォルトパラメータに依存
blastn -query test.fasta -subject db.fasta -outfmt 7 > output.txt

# ✅ 良い例: 全てのパラメータを明示
blastn -query test.fasta -subject db.fasta \
       -dust no \
       -reward 1 -penalty -2 \
       -gapopen 0 -gapextend 0 \
       -outfmt 7 > output.txt
```

## 修正版アクションプラン

### 1. コマンドの厳格化

正解データ生成コマンドに以下を**必ず**追加：

**BLASTN**:
```bash
-dust no \                    # フィルタリングOFF
-reward 1 -penalty -2 \       # スコアリングパラメータ明示
-gapopen 0 -gapextend 0       # ギャップコスト明示
```

**TBLASTX**:
```bash
-comp_based_stats 0 \         # 組成補正OFF
-seg no \                     # フィルタリングOFF
-gapopen 11 -gapextend 1 \    # ギャップコスト明示
-matrix BLOSUM62              # スコアリング行列明示
```

### 2. アサーションの厳格化（実装が移植済みのため）

Effective Search Spaceのテストは、実装が移植済みであるため、非常に厳格に：

```rust
// 0.1%の相対誤差を許容（浮動小数点精度の違いのみ）
// より大きな差は実装のバグを示す
assert_relative_eq!(
    losat_effective_space,
    ncbi_effective_space,
    epsilon = 0.001  // 0.1%
);
```

**注意**: LOSATの実装はNCBI BLASTの直接移植であるため、完全一致または非常に近い値が期待されます。
大きな差（>1%）が見られる場合は、実装のバグを調査する必要があります。

### 3. テストの分離

以下の2つを**完全に分離**：

1. **アライメント拡張ロジックのテスト**:
   - X-dropの動作
   - アライメントの終了位置
   - ギャップ処理

2. **スコア・E-value計算ロジックのテスト**:
   - 固定アライメント座標を使用
   - 計算式の正確性のみを検証

## まとめ

このテスト戦略により、以下の問題を回避できます：

1. ✅ 実効長計算の複雑さ → 許容誤差で対応
2. ✅ X-dropのバタフライ効果 → テストの分離で対応
3. ✅ フィルタリングの影響 → コマンドの厳格化で対応
4. ✅ バージョン差異 → バージョン固定とパラメータ明示で対応

**重要な原則**: 
- 完全一致を目指さない（実装コストが高すぎる）
- 許容誤差を設ける（統計的に意味のある範囲内）
- テストを分離する（異なる関心事を混ぜない）
- パラメータを明示する（再現性を確保）

