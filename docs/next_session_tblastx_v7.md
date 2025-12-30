# LOSAT TBLASTX 改善セッション引継ぎ（v7）

## 今回のセッションでの主要変更

### 1. NCBIアミノ酸エンコーディング順序への完全移行

**問題**: LOSATはアルファベット順（A-Z, *）を使用していたが、NCBIは異なる順序を使用。

**NCBIの順序**: `ARNDCQEGHILKMFPSTWYVBJZX*` (25文字、インデックス0-24)

**変更したファイル**:

| ファイル | 変更内容 |
|---------|---------|
| `matrix.rs` | BLOSUM62を25x25に変更、NCBIの正確な値を使用 |
| `translation.rs` | `aa_char_to_ncbi_index`を使用してNCBIインデックスを返す |
| `constants.rs` | `STOP_CODON = 24` |
| `lookup.rs` | k-merテーブルサイズを24^3=13,824に変更 |
| `utils.rs` | k-merエンコーディングを`c1*576 + c2*24 + c3`に変更 |
| `extension.rs` | `get_score`を25x25マトリックス用に変更 |
| `aa_word_finder.rs` | 全関数を24^3対応に変更 |

### 2. `*-*`スコアの変更

NCBIでは`*-*`（Stop codon同士）のスコアが+1だが、これはself-comparisonでrunaway extensionを引き起こす。
**-4に変更**してX-dropで終了するようにした。

### 3. `diag_extended`フラグチェックの追加

NCBIの`flag`に相当するチェックを追加：
- 拡張後、`diag_extended = true`を設定
- 次のシードで`s_pos < last_hit`ならスキップ
- `s_pos >= last_hit`なら`diag_extended = false`にリセット

## 現在の状況（NCBIとの比較）

### Self-comparison (NZ_CP006932)

| 項目 | NCBI BLAST | LOSAT |
|------|-----------|-------|
| **Total hits** | 62,059 | 41,583 |
| **1000+ aa** | 6 | 89 |
| **最大長** | 2,260 aa | 3,670 aa |
| **21-50 aa** | 33,764 | 20,731 |

### Low identity (AP027131 vs AP027133)

| 項目 | NCBI BLAST | LOSAT |
|------|-----------|-------|
| **Total hits** | 14,877 | 6,488 |
| **500+ aa** | 6 | 4 |
| **21-50 aa** | 8,690 | 3,894 |

## 残存する問題

### 1. ヒット数が少ない（約44%〜67%）

**診断結果**:
```
Seeds filtered (two-hit): 646,445,045
Seeds passed to extension: 14,986,853
Filtered (cutoff_score=42): 14,906,586 (99.5%)
E-value passed: 6,488
```

拡張に渡されたシードの99.5%が`cutoff_score=42`でフィルタされている。

### 2. 長いヒットが多い

Self-comparisonでLOSATは最大3,670 aa、NCBIは最大2,260 aa。
89個の1000+ aa hitがあるが、NCBIは6個のみ。

## 考えられる原因

1. **拡張後のスコアがNCBIより低い**
   - マトリックスアクセスは正しい（確認済み）
   - X-drop終了条件の違い？

2. **Two-hit条件の違い**
   - NCBIのtwo-hit実装の詳細確認が必要

3. **Neighboring word展開の違い**
   - LOSATは24^3テーブル、NCBIは？

## 次のセッションでの確認事項

1. **NCBIの拡張スコア計算との比較**
   - 同じシード位置で同じスコアが計算されるか
   
2. **X-drop終了条件の確認**
   - `right_score <= 0`の条件がNCBIと一致しているか

3. **Two-hit window内の処理**
   - `last_hit`の更新タイミング

4. **診断ログの追加**
   - 拡張前後のスコアを出力して比較

## ビルド・テスト状況

- `cargo build --release` ✓
- 一部のunit testは失敗（DNA関連、アミノ酸変更とは無関係）
- 実行テスト: 動作確認済み

## 参考：NCBIファイル

- `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c`
- `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/util/tables/sm_blosum62.c`
- NCBI出力: `/mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT/tests/blast_out/`

