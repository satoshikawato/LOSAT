# 長い配列での過剰ヒット問題の調査手順

このドキュメントは、長い配列（600kb+）での過剰ヒット問題を調査するための手順を説明します。

## 問題の概要

- **現象**: AP027131 vs AP027133 (600kb+): LOSAT 29,766 vs NCBI 14,871 (約2倍)
- **根本原因**: HSP生成数が過多 (LOSAT 338,859 vs NCBI 推定 30,000-45,000、約8-11倍)
- **影響**: 低スコアヒット (bit < 30) が NCBI の 2.5 倍

## 調査手順

### ステップ1: テストケースの準備

長い配列のFASTAファイルを準備します：

```bash
# 例: AP027131 vs AP027133
QUERY="AP027131.fna"
SUBJECT="AP027133.fna"
GENCODE=4  # バクテリア用
```

### ステップ2: LOSATとNCBI BLAST+の実行

比較スクリプトを実行します：

```bash
cd tests
./compare_long_sequences_debug.sh "$QUERY" "$SUBJECT" "$GENCODE"
```

このスクリプトは以下を実行します：
1. LOSAT TBLASTXを実行（デバッグ出力付き）
2. NCBI BLAST+ TBLASTXを実行（可能な場合）
3. デバッグ情報を抽出
4. 統計情報を表示

### ステップ3: デバッグ出力の分析

詳細な分析を実行します：

```bash
python analyze_debug_output.py \
    debug_results_*/losat_debug.log \
    debug_results_*/ncbi_debug.log
```

### ステップ4: 比較項目の確認

以下の項目を比較します：

#### 1. Cutoff値の比較

**LOSATデバッグ出力から**:
- `[DEBUG CUTOFF_CALC]`: `eff_searchsp`, `cutoff_score_max`の計算過程
- `[DEBUG CUTOFF_UPDATE]`: 最終cutoff値の決定過程

**確認ポイント**:
- `eff_searchsp`がNCBIと一致しているか
- `cutoff_score_max`がNCBIと一致しているか
- 最終cutoff値（`MIN(update_cutoff, gap_trigger, cutoff_score_max)`）が正しいか

#### 2. Linking Cutoffの比較

**LOSATデバッグ出力から**:
- `[DEBUG LINKING_CUTOFF]`: `cutoff_small_gap`, `cutoff_big_gap`の計算
- `[DEBUG LINKING_FILTER]`: フィルタリング率

**確認ポイント**:
- `search_sp > 8 * window_sq`の条件が正しいか
- `cutoff_small_gap`と`cutoff_big_gap`がNCBIと一致しているか
- フィルタリング率が適切か

#### 3. HSP保存統計の確認

**LOSATデバッグ出力から**:
- `[DEBUG HSP_SAVING]`: 保存数、フィルタリング数、スコア分布

**確認ポイント**:
- 低スコアHSP (bit < 30) の割合
- cutoffによるフィルタリング率
- スコア分布の異常

#### 4. ヒット数の比較

**出力ファイルから**:
- LOSAT: `debug_results_*/losat_output.tsv`
- NCBI: `debug_results_*/ncbi_output.tsv`

**確認ポイント**:
- 総ヒット数の差
- スコア分布の差
- 低スコアヒットの数

## 期待される発見

以下のいずれかが原因である可能性があります：

1. **Cutoff計算の差異**: `eff_searchsp`や`cutoff_score_max`の計算がNCBIと異なる
2. **Linking Cutoffのスケーリング**: 長い配列での`cutoff_small_gap`/`cutoff_big_gap`の計算が不適切
3. **HSP保存条件**: `score >= cutoff`の条件が正しく適用されていない
4. **フィルタリング率**: Linking段階でのフィルタリングが不十分

## トラブルシューティング

### LOSATのデバッグ出力が表示されない場合

長い配列（600kb+）でのみデバッグ出力が有効になります。配列の長さを確認してください：

```bash
# 配列の長さを確認
grep -v "^>" "$QUERY" | tr -d '\n' | wc -c
grep -v "^>" "$SUBJECT" | tr -d '\n' | wc -c
```

### NCBI BLAST+がインストールされていない場合

NCBI BLAST+の出力は必須ではありませんが、比較には有用です。インストール方法：

```bash
# Ubuntu/Debian
sudo apt-get install ncbi-blast+

# macOS
brew install blast

# またはNCBIからダウンロード
# https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
```

## 次のステップ

比較結果を基に、以下のいずれかを実施します：

1. **Cutoff計算の修正**: `ncbi_cutoffs.rs`の修正
2. **Linking Cutoffの修正**: `sum_stats_linking.rs`の修正
3. **HSP保存条件の修正**: `utils.rs`の修正
4. **追加の調査**: より詳細なデバッグ出力の追加

## 関連ファイル

- `compare_long_sequences_debug.sh`: 比較スクリプト
- `analyze_debug_output.py`: デバッグ出力分析スクリプト
- `src/algorithm/tblastx/ncbi_cutoffs.rs`: Cutoff計算
- `src/algorithm/tblastx/sum_stats_linking.rs`: Linking Cutoff計算
- `src/algorithm/tblastx/utils.rs`: HSP保存ロジック

