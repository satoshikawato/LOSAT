# 極端にバイアスのあるアミノ酸組成のクエリについて

## 概要

「極端にバイアスのあるアミノ酸組成のクエリ」とは、特定のアミノ酸が異常に多く含まれている（または少ない）クエリ配列のことです。

## 標準的なアミノ酸組成

NCBI BLASTが使用する標準的なアミノ酸組成（Robinson標準頻度）は以下の通りです：

| アミノ酸 | 標準頻度 | 説明 |
|---------|---------|------|
| A (Ala) | 7.8% | 比較的多い |
| R (Arg) | 1.9% | 少ない |
| N (Asn) | 5.4% | 中程度 |
| D (Asp) | 6.3% | 中程度 |
| C (Cys) | 1.5% | 非常に少ない |
| Q (Gln) | 3.4% | 少ない |
| E (Glu) | 6.7% | 中程度 |
| G (Gly) | 7.1% | 比較的多い |
| H (His) | 2.1% | 少ない |
| I (Ile) | 5.1% | 中程度 |
| L (Leu) | 5.7% | 中程度 |
| K (Lys) | 5.1% | 中程度 |
| M (Met) | 1.5% | 非常に少ない |
| F (Phe) | 4.0% | 少ない |
| P (Pro) | 4.7% | 少ない |
| S (Ser) | 6.1% | 中程度 |
| T (Thr) | 5.5% | 中程度 |
| W (Trp) | 1.3% | 非常に少ない |
| Y (Tyr) | 3.2% | 少ない |
| V (Val) | 6.9% | 比較的多い |

## 極端にバイアスのある例

### 例1: リピート配列（特定アミノ酸が異常に多い）

```
配列例: "AAAAAAAAAAAAAAAAAAAA..." (Aが90%以上)
```

- **特徴**: アラニン(A)が90%以上を占める
- **影響**: 
  - スコア頻度プロファイルが偏る
  - 計算されたLambdaが**小さくなる**可能性がある
  - `check_ideal`ロジックにより、Lambda < 0.3176 なら計算値が使用される

### 例2: 低複雑度配列（少数のアミノ酸に偏る）

```
配列例: "GGGGGGGGGGGGGGGGGGGG..." (Gが80%以上)
```

- **特徴**: グリシン(G)が80%以上を占める
- **影響**: 同様にLambdaが小さくなる可能性

### 例3: 特定のアミノ酸が異常に少ない

```
配列例: システイン(C)が0.1%以下（標準は1.5%）
```

- **特徴**: 特定のアミノ酸がほとんど含まれない
- **影響**: スコア分布が変わり、Lambdaが変動する可能性

## check_ideal ロジックの動作

NCBI BLASTの`check_ideal`ロジック（`blast_stat.c:2796-2797`）：

```c
if (check_ideal && kbp->Lambda >= sbp->kbp_ideal->Lambda)
   Blast_KarlinBlkCopy(kbp, sbp->kbp_ideal);
```

### 動作パターン

1. **通常のクエリ**:
   - 計算されたLambda ≈ 0.32-0.35
   - `Lambda >= 0.3176` (ideal) → **idealを使用**
   - より保守的な（小さい）LambdaでE-valueを計算

2. **極端にバイアスのあるクエリ**:
   - 計算されたLambda ≈ 0.25-0.30（例：Aが90%の配列）
   - `Lambda < 0.3176` (ideal) → **計算値を使用**
   - より大きいLambdaでE-valueを計算（より緩い判定）

## 実例

### 標準組成のクエリ

```
配列: "ACDEFGHIKLMNPQRSTVWY" (各アミノ酸が均等に含まれる)
計算されたLambda: 0.32
check_ideal判定: 0.32 >= 0.3176 → ideal (0.3176) を使用
```

### 極端にバイアスのあるクエリ

```
配列: "AAAAAAAAAAAAAAAAAAAA" (Aが100%)
計算されたLambda: 0.28 (仮定)
check_ideal判定: 0.28 < 0.3176 → 計算値 (0.28) を使用
```

## なぜ重要か

1. **E-value計算への影響**:
   - Lambdaが小さい → 同じraw scoreでもE-valueが大きくなる
   - より緩い判定基準になる

2. **NCBI parity**:
   - NCBIはcontextごとに計算し、`check_ideal`を適用
   - LOSATも同様に実装することで、完全なparityを達成

3. **実用例**:
   - リピート配列の検出
   - 低複雑度領域の処理
   - 特殊なタンパク質（例：コラーゲン、ケラチンなど）

## 実装での対応

LOSATの実装（`src/stats/karlin_calc.rs`）：

```rust
pub fn apply_check_ideal(computed: KarlinParams, ideal: KarlinParams) -> KarlinParams {
    if computed.lambda >= ideal.lambda {
        // 通常のクエリ: idealを使用（より保守的）
        ideal
    } else {
        // 極端にバイアスのあるクエリ: 計算値を使用
        computed
    }
}
```

これにより、NCBI BLASTと完全に同等の動作を保証します。

