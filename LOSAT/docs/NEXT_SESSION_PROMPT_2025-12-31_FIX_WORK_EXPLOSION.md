# 次セッションプロンプト: TBLASTX 仕事量爆発の修正

## 絶対命令

**Lookup entries の爆発（0.7M → 8.1M）を修正し、`NZ_CP006932 self` (full 6×6, single-thread) を 30秒未満で完了させる。**

## 前提

- NCBI BLAST+ (C/C++) が ground truth。出力を一致させる。
- 出力に影響を与える"簡略化"は禁止。
- 副作用・計算量を減らせる場合は変更を許容。
- 大規模な変更も躊躇なく行え。NCBI BLAST+ (C/C++) の ground truthに合わせろ。

## 現状の問題

1. **実行が10分以上かかる（ハング状態）**: 30秒目標に対して論外
2. **Lookup entries が 8.1M に爆発**: 以前の 0.7M から 10倍以上
3. **プログラムが "Searching..." で止まる**: scan loop が終わらないか、MPSC channel deadlock の可能性

## 最優先で確認・修正すべき箇所

### 1) Lookup entries 爆発の原因

`LOSAT/src/algorithm/tblastx/lookup.rs` の `build_ncbi_lookup()`:

- neighbor 生成で **同じ offset を重複追加**していないか？
- `exact_offsets` への追加と neighbor 展開の両方で同じ offset が入っていないか？
- `all_entries[nidx].extend_from_slice(offsets)` が複数回呼ばれて爆発していないか？

### 2) NCBI alphabet size との差異

```c
// NCBI: BLASTAA_SIZE = 28
// LOSAT: alphabet_size = 25
```

- NCBI は 28文字（0-27）: `-ABCDEFGHIKLMNPQRSTVWXYZUOJ*`
- LOSAT の `matrix.rs` は 25文字（0-24）: `ARNDCQEGHILKMFPSTWYVBJZX*`

**この差異がエンコード/デコードに影響していないか確認。**

### 3) SEG masking の適用タイミング

- SEG filter が seed 生成**前**に適用されているか？
- マスクされた領域から seed が生成されていないか？

### 4) MPSC channel の deadlock

`LOSAT/src/algorithm/tblastx/utils.rs` の `run()`:

- `drop(tx)` が `par_iter().for_each_init()` の**後**に正しく呼ばれているか？
- writer thread が `rx.recv()` で永遠に待機していないか？

## 参照ファイル

- `/mnt/c/Users/kawato/Documents/GitHub/LOSAT/src/algorithm/tblastx/lookup.rs` - lookup 生成
- `/mnt/c/Users/kawato/Documents/GitHub/LOSAT/src/algorithm/tblastx/utils.rs` - scan/extension/channel
- `/mnt/c/Users/kawato/Documents/GitHub/LOSAT/src/utils/matrix.rs` - BLOSUM62 matrix (25x25)
- `/mnt/c/Users/kawato/Documents/GitHub/LOSAT/docs/TBLASTX_NCBI_PARITY_DEVLOG.md` - 経緯

## NCBI 参照コード

- `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/include/algo/blast/core/blast_encoding.h` - `BLASTAA_SIZE = 28`
- `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c` - lookup 生成
- `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c` - scan loop

## 運用ルール
- **調査ではなく修正を最優先**。コードを読んだら即座に直す。
- **10分毎に `TBLASTX_NCBI_PARITY_DEVLOG.md` へ追記**して思考を保存
- パッチが当たらない場合は素直にファイルを読み直して追記


## ベンチマーク

```bash
cd /mnt/c/Users/kawato/Documents/GitHub/LOSAT/LOSAT/tests
cargo build --release
time LOSAT_DIAGNOSTICS=1 ../target/release/LOSAT tblastx \
  -q ./fasta/NZ_CP006932.fasta \
  -s ./fasta/NZ_CP006932.fasta \
  -o ./losat_out/test.out \
  --query-gencode 4 --db-gencode 4 -n 1
```

目標: **30秒未満**

