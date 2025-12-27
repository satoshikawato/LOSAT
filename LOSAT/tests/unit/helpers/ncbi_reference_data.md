# NCBI BLAST参照データの準備方法（改訂版）

## 重要なテスト戦略の変更

単純にλとKを逆算するだけでは不十分です。以下の4つの問題に対処する必要があります：

### 1. 実効探索空間（Effective Search Space）の検証 ⚠️ 最重要

**問題**: NCBI BLASTは長さ調整（length adjustment）を行い、実際の文字数ではなく「実効長」を使用します。

**対策**: 
- NCBI BLASTの出力ヘッダーから `Effective search space` を抽出
- LOSATの計算結果と比較
- これが最も重要な検証項目です

**実装例**:
```rust
// NCBI BLASTヘッダーから抽出
# Effective search space used: 1234567890

// LOSATの計算結果と比較
let losat_space = SearchSpace::for_database_search(...);
assert_eq!(losat_space.effective_space, expected_space);
```

### 2. Composition-based Statistics の回避

**問題**: TBLASTX/BLASTPでは、デフォルトでアミノ酸組成に基づく補正がかかります。

**対策**: 
- NCBI BLAST実行時に `-comp_based_stats 0` を指定
- 純粋な統計値（固定のλ, K）と比較

**実行例**:
```bash
tblastx -query test.fasta -subject db.fasta \
        -comp_based_stats 0 \  # ← 重要！
        -outfmt 7 \
        -out ncbi_output.txt
```

### 3. 複数のパラメータセットでテスト

**問題**: λとKはスコアリングパラメータに依存します。

**対策**: 
- 複数のスコアリングパラメータセットでテスト
- デフォルト値だけでなく、様々な組み合わせを検証

**テストケース例**:
- Megablast: reward=1, penalty=-2
- Blastn task: reward=1, penalty=-3
- カスタム: reward=2, penalty=-5

### 4. 様々な長さの配列でテスト

**問題**: 短い配列では漸近理論が成り立たない可能性があります。

**対策**: 
- 短い配列（100bp）、中程度（1kbp）、長い配列（1Mbp）でテスト
- 長さ調整が正しく機能するか確認

## データ抽出の手順（改訂版）

### ステップ1: NCBI BLASTを実行（重要オプション付き）

**⚠️ 重要**: 以下のオプションを**必ず**指定してください：

```bash
# BLASTNの場合
blastn -query test_query.fasta \
       -subject test_subject.fasta \
       -dust no \  # ← 必須: フィルタリングをOFF（DUSTフィルタの影響を排除）
       -reward 1 -penalty -2 \  # ← 必須: パラメータを明示的に指定
       -gapopen 0 -gapextend 0 \  # ← 必須: ギャップコストを明示
       -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -out ncbi_output.txt

# TBLASTXの場合
tblastx -query test_query.fasta \
        -subject test_subject.fasta \
        -comp_based_stats 0 \  # ← 必須: 組成補正をOFF
        -seg no \  # ← 必須: フィルタリングをOFF（SEGフィルタの影響を排除）
        -gapopen 11 -gapextend 1 \  # ← 必須: ギャップコストを明示
        -matrix BLOSUM62 \  # ← 必須: スコアリング行列を明示
        -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -out ncbi_output.txt
```

**バージョン固定**: READMEに「NCBI BLAST+ 2.X.X を使用」と明記し、同じバージョンで統一してください。

### ステップ2: ヘッダーから実効探索空間を抽出

NCBI BLASTの出力ヘッダーを確認：
```
# BLASTN 2.17.0+  ← バージョン情報も記録
# Query: ...
# Database: ...
# Effective search space used: 1234567890  ← これを抽出
```

**✅ 実装状況：実効長計算は既に移植済み**

LOSATは既にNCBI BLASTの`BLAST_ComputeLengthAdjustment`関数をRustに移植しています：
- `src/stats/length_adjustment.rs`の`compute_length_adjustment_ncbi`関数
- NCBI BLASTの`blast_stat.c`の実装を直接移植

**対策**:
- **実装の検証**: NCBI BLASTの実装と1行ずつ比較し、完全一致を確認
- **小さな許容誤差**: 浮動小数点精度の違いを考慮し、**0.1%の相対誤差**を許容
- **バグ検出**: より大きな差（>1%）が見られる場合は、実装のバグの可能性があるため修正が必要

### ステップ3: テストケースの作成

```rust
NcbiEvalueTestCase {
    score: 100,
    q_len: 1000,
    db_len: 10000,
    db_num_seqs: 5,
    params: KarlinParams { ... },
    expected_bit_score: 190700.0,
    expected_evalue: 0.0,
    
    // ⚠️ 重要: 実効探索空間を指定
    expected_effective_space: Some(1234567890.0),  // NCBI BLASTヘッダーから
    expected_length_adjustment: Some(42),           // 計算またはNCBI BLASTから
    
    tolerance: 0.1,
    test_name: "megablast_default".to_string(),
}
```

## テストの役割分担

### ユニットテスト（このモジュール）
- **目的**: 計算ロジックの正確性を検証
- **検証項目**:
  - Bit score計算の正確性
  - E-value計算の正確性
  - **実効探索空間計算の正確性**（最重要）
  - パラメータセットごとの一貫性
  - 様々な長さでの動作

### 統合テスト（tests/run_comparison.sh）
- **目的**: エンドツーエンドのNCBI BLAST互換性を検証
- **検証項目**:
  - 実際の配列でのヒット数
  - E-value分布
  - Bit score分布
  - アライメント座標

## 推奨されるテストケース構成

### BLASTN用テストケース
1. **パラメータセット1**: Megablast (reward=1, penalty=-2)
   - 短い配列 (100bp)
   - 中程度 (1kbp)
   - 長い配列 (1Mbp)
2. **パラメータセット2**: Blastn task (reward=1, penalty=-3)
   - 同様の長さバリエーション
3. **データベース構成**:
   - 単一配列
   - 複数配列（5, 10, 100配列）

### TBLASTX用テストケース
1. **BLOSUM62** (gap_open=11, gap_extend=1)
   - 短いアライメント (50aa)
   - 中程度 (200aa)
   - 長いアライメント (1000aa)
2. **-comp_based_stats 0** で実行したデータのみ使用

## 注意事項

1. **E-value = 0.0 の扱い**
   - NCBI BLASTは `< 1e-180` を `0.0` と表示
   - テストでは `actual < 1e-180` で検証

2. **実効探索空間の取得**
   - NCBI BLASTのヘッダーから抽出（可能な場合）
   - または、LOSATの計算結果を手動で検証

3. **パラメータの一致**
   - LOSATとNCBI BLASTで同じスコアリングパラメータを使用
   - Karlin-Altschulパラメータも一致させる

4. **許容誤差**
   - E-value: 10-20%の相対誤差
   - Bit score: 0.01-0.1の絶対誤差
   - **実効探索空間: 1%の相対誤差**（完全一致は不可能なため、許容誤差を設ける）
   - 長さ調整値: 小さい値（<100）では1単位、大きい値では1%の相対誤差

5. **テストの分離**
   - **アライメント拡張ロジック**（どこまで伸びるか）のテスト
   - **スコア・E-value計算ロジック**（伸びた結果をどう評価するか）のテスト
   - これらを混ぜない（X-dropのバタフライ効果を避けるため）

6. **バージョン固定**
   - 正解データ生成には**NCBI BLAST+ 2.X.X**を使用（バージョンを明記）
   - パラメータを全てコマンドライン引数で明示的に指定して固定
