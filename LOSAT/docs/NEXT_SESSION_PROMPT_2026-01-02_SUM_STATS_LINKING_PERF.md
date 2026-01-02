# Next Session Prompt: sum_stats_linking パフォーマンス改善

## 現状

### 完了した作業（2026-01-02 セッション）

1. **sum_stats_linking NCBI完全準拠修正**
   - `LhHelper`に`maxsum1`フィールド追加
   - センチネル初期化を`sum: [0, 0]`に修正（NCBI calloc動作）
   - `can_skip_ncbi`条件から余分な`linked_to >= 0`チェック削除
   - `next_larger`ループ条件を`prev > 0`に修正
   - `maxsum1`計算を3箇所に追加（構築時、INDEX 1後、lh_helpers[1]初期化）

2. **HSP重複除去関数追加**
   - `purge_hsps_with_common_endpoints`実装（NCBI `Blast_HSPListPurgeHSPsWithCommonEndpoints`相当）
   - 効果: 232,359 → 232,275 (84個削除) - self-comparisonでは限定的

### 現在の問題

| 指標 | 値 |
|------|-----|
| 生成HSP数 | 232,359 (self-comparison) |
| グループ数 | 4 (++, +-, -+, --) |
| グループあたりHSP | ~58,000 |
| INDEX 0計算量 | O(N²) worst case |
| 実行時間 | 3分38秒 |

**根本原因**: self-comparisonで対角線d=0上に巨大な完全一致領域があり、周辺対角線にも多数のヒットが集中。INDEX 0のinner loopでearly break (`qo > h_qe_gap + TRIM_SIZE`)が効かず、O(N²)が支配的。

**NCBIでも同様に遅くなる可能性が高い**:
- `hsp_num_max`はデフォルトINT4_MAX（無制限）
- INDEX 0には`next_larger`最適化がない
- early breakは同様に効かない

## 次セッションのタスク



### 優先度1: INDEX 0最適化の検討（NCBI非準拠オプション）

**目的**: INDEX 0のO(N²)をO(N log N)またはO(N)に改善

**検討事項**:
1. **INDEX 0に`next_larger`最適化を追加**
   - INDEX 1と同様のジャンプ最適化
   - NCBIにはないが、大幅な高速化が期待できる
   - オプション化（`--optimize-index0`）でNCBI準拠モードと切り替え可能に

2. **空間インデックス（KD-tree/R-tree）の導入**
   - `(q_off_trim, s_off_trim)`空間での近傍探索
   - 実装コストが高い

3. **早期終了条件の強化**
   - `h_sum`が一定値以上の場合、inner loopを早期終了
   - 出力に影響する可能性があるため慎重に検討

**制約**:
- **出力同等性を絶対維持**: 最適化は出力に影響を与えてはならない
- NCBI準拠モードと最適化モードを明確に分離
- 最適化モードは明示的なオプションで有効化

### 優先度2: upstream HSP数削減の検討

**目的**: extension段階でHSP数を削減

**検討事項**:
1. **cutoff計算の見直し**
   - 現在のcutoffが適切か確認
   - NCBIのcutoff計算と比較

2. **two-hit window制限の強化**
   - `window_size`パラメータの調整
   - ただし、NCBI準拠を維持する必要がある

3. **diag_arrayの早期クリア**
   - 一定数のHSPが生成されたら、低スコアの対角線をクリア
   - 出力に影響する可能性があるため慎重に

### 優先度3: 診断機能の強化

**目的**: パフォーマンスボトルネックの可視化

**実装内容**:
- INDEX 0/INDEX 1のinner loop反復回数カウント
- `next_larger`ジャンプ回数カウント
- `use_current_max`有効化回数カウント
- グループサイズ分布の出力
- タイミング情報（各フェーズの実行時間）

**環境変数**: `LOSAT_DIAGNOSTICS=1`で有効化

## 制約事項

### 絶対に守るべき原則

1. **NCBI実装が唯一の正解 (GROUND TRUTH)**
   - 出力を1ビットの狂いもなく一致させること
   - 出力に影響を与える"簡略化"は厳禁

2. **出力同等性の絶対維持**
   - 最適化は出力に影響を与えてはならない
   - 副作用・計算量を減らせる場合は変更を許容するが、出力同等性は絶対維持

3. **アルゴリズムの計算量 (Big O) と論理フローは維持**
   - 性能上の理由でRustらしく変形（借用チェッカーへの適合など）するのは許容
   - ただし、アルゴリズムの核心部分は維持

4. **推測を排除**
   - NCBIソースとLOSATソースを対比して差異を特定
   - 不一致を発見した場合は即座に修正

### 許容される変更

- **オプション化された最適化**: NCBI準拠モードと最適化モードを分離
- **診断機能の追加**: 出力に影響しない計測・ログ機能
- **HSP数制限オプション**: デフォルト無制限でNCBI準拠、明示的指定時のみ制限

## 参考ファイル

### NCBI参照コード

| ファイル | 関数/構造体 | 行番号 | 説明 |
|----------|-------------|--------|------|
| link_hsps.c | `s_BlastEvenGapLinkHSPs` | 400-1000 | メインリンキング関数 |
| link_hsps.c | `LinkHelpStruct` | 100-109 | ヘルパー構造体 |
| link_hsps.c | `LinkHSPStruct` | 76-94 | HSPリンク構造体 |
| link_hsps.c | INDEX 0ループ | 690-768 | 小ギャップリンキング |
| link_hsps.c | INDEX 1ループ | 770-896 | 大ギャップリンキング（next_larger最適化あり） |
| blast_hits.c | `Blast_HSPListPurgeHSPsWithCommonEndpoints` | 2455-2534 | HSP重複除去 |

### LOSAT実装ファイル

- `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs`: sum_stats_linking実装
- `LOSAT/src/algorithm/tblastx/utils.rs`: HSP生成と重複除去
- `LOSAT/src/algorithm/tblastx/args.rs`: コマンドライン引数定義

## 成功基準

1. **パフォーマンス改善（オプション）**
   - INDEX 0最適化が実装された場合、出力同等性を維持
   - 診断機能で改善効果を確認可能

2. **ドキュメント更新**
   - DEVLOGに実装内容を追記

## 注意事項

- **テスト実行は最後に**: 実装完了後にテストを実行
- **段階的実装**: 優先度1から順に実装
- **NCBI準拠を最優先**: 性能改善よりも正確性を優先

