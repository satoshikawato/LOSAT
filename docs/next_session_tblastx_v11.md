# TBLASTX実装乖離調査 - セッション11への引継ぎ

## 現在の状況サマリー

### ヒット数比較
| 実装 | ヒット数 | 処理時間 | 備考 |
|------|----------|----------|------|
| NCBI BLAST+ (SEG有り) | 62,053 | - | 目標値 |
| LOSAT frame-based (SEGなし) | 33,873 | ~1分25秒 | strand_group_v2.out |
| LOSAT strand-based (SEGなし) | 9,077 | 2分以上 | 最新実装、遅すぎる |

### 主要な問題

#### 1. sum_stats_linkingの速度問題
- **原因**: strand-basedグループ分け（4グループ）により、1グループあたり約11万HSP
- **影響**: O(N²)アルゴリズムで12兆回の比較 → 処理時間が爆発
- **NCBIの解決策**: `lh_helper`配列と`next_larger`ポインタによる最適化
- **現状**: next_larger最適化を実装したが、まだ遅い（whileループ内で毎回再計算している）

#### 2. SEGフィルタのリリースビルドクラッシュ
- **症状**: `--seg`オプション有効時にリリースビルドでSegfault
- **デバッグビルドでは動作する** → 最適化関連のバグ
- **暫定対応**: SEGのデフォルトをfalseに変更（args.rs）
- **未解決**

#### 3. ヒット数の乖離
- frame-basedで33,873ヒット（NCBIの54%）
- strand-basedで9,077ヒット（NCBIの15%）
- whileループ（複数チェーン抽出）を実装したが、効果未確認（速度問題で検証不能）

## 技術的詳細

### NCBIのsum_stats_linking最適化（link_hsps.c）

```c
// lh_helper構造体（100-109行）
typedef struct LinkHelpStruct {
  LinkHSPStruct* ptr;
  Int4 q_off_trim;
  Int4 s_off_trim;
  Int4 sum[2];        // small_gap, large_gap
  Int4 maxsum1;       // 早期終了用閾値
  Int4 next_larger;   // 最適化の鍵：より大きいsumを持つHSPへのポインタ
} LinkHelpStruct;

// next_largerの計算（675-684行）
Int4 cur_sum = lh_helper[H_index].sum[1];
Int4 prev = H_index - 1;
while ((cur_sum >= prev_sum) && (prev > 0)) {
    prev = lh_helper[prev].next_larger;
    prev_sum = lh_helper[prev].sum[1];
}
lh_helper[H_index].next_larger = prev;

// リンキングループでの使用（827-860行）
for (H2_index = H_index-1; H2_index > 1;) {
    sum = H2_helper->sum[index];
    next_larger = H2_helper->next_larger;
    
    b0 = sum <= H_hsp_sum;  // sumが小さすぎる
    
    H2_index--;
    if (b0) {  // sumが小さい場合、next_largerにジャンプ
        H2_index = next_larger;
    }
    
    if (!(b0|b1|b2)) {  // 全条件クリア
        // リンク候補として採用
    }
}
```

### 現在のLOSAT実装の問題点

1. **whileループ内でnext_largerを毎回再計算**
   - NCBIは初回のみ計算し、HSP削除時に更新
   - LOSATは毎イテレーションで全再計算 → O(N²)が解消されていない

2. **リンク情報のリセットが毎回発生**
   - `helpers[i].sum = helpers[i].score` を毎回実行
   - NCBIは変更されたHSPのみ再計算（`changed`フラグ使用）

3. **active配列のスキャンがO(N)**
   - best chain検索で全HSPをスキャン
   - NCBIはlinked listで効率的に管理

## ファイル構成

### 変更されたファイル
- `src/algorithm/tblastx/sum_stats_linking.rs` - メインのリンキング実装
- `src/algorithm/tblastx/args.rs` - SEGデフォルトをfalseに変更
- `src/algorithm/tblastx/utils.rs` - mask_keyにs_f_idx追加（以前の修正）

### 参照すべきNCBIファイル
- `ncbi-blast/c++/src/algo/blast/core/link_hsps.c` - リンキング実装
- `ncbi-blast/c++/src/algo/blast/core/blast_parameters.c` - cutoff計算

### テスト出力ファイル
- `tests/losat_out/strand_group_v2.out` - frame-based、33,873ヒット
- `tests/losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.out` - 最新、9,077ヒット
- `tests/blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out` - NCBI参照、62,053ヒット

## 今後の方針

### 優先度1: sum_stats_linkingの最適化
1. **NCBIの実装を完全に再現する**
   - `changed`フラグによる差分更新
   - next_largerの初回計算と差分更新
   - linked list的な構造でactive HSPを管理

### 優先度2: SEGクラッシュの修正
- デバッグビルドで動作確認済み
- リリースビルドの最適化オプションを調整
- または、SEG実装のバグを特定

### 優先度3: ヒット数の検証
- 速度問題解決後、ヒット数を再検証
- E-value分布をNCBIと比較

## 検証コマンド

```bash
# LOSAT実行
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT
cargo build --release
time ./target/release/losat tblastx \
  -q tests/fasta/NZ_CP006932.fasta \
  -s tests/fasta/NZ_CP006932.fasta \
  -o tests/losat_out/test.out

# 結果確認
wc -l tests/losat_out/test.out
wc -l tests/blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out

# E-value分布比較
awk -F'\t' '{print $11}' tests/losat_out/test.out | sort | uniq -c | sort -rn | head -10
```

## 重要な教訓

1. **LLMの弱点**: トークン列に引きずられて俯瞰的視点を失いやすい
   - コード修正に夢中になり、「そもそもの設計問題」を見落とす
   - 「ディスプレイがつかない」→電源確認という発想が必要

2. **グループサイズの影響**: O(N²)アルゴリズムではグループサイズが致命的
   - 36グループ→4グループで処理時間81倍以上

3. **NCBIの最適化は必須**: strand-basedグループ分けにはlh_helper最適化が不可欠

