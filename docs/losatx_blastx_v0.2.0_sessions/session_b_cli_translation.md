# Session B — CLI・翻訳・query 準備

## INSTRUCTION PROMPT

LOSATX v0.2.0 の段階 B を **`feature/losatx-blastx-v0.2.0`** で実行する。
ブランチ・作業ツリー、`AGENTS.md`、`verify-ncbi-parity-and-speed`、
[総合計画書](../losatx_blastx_v0.2.0_plan.md)の X01〜X05・X09・X14・
第 4〜8 節、および段階 A の source map・defaults・fixture 証拠を確認する。
A の必須差分が残るなら先に解消する。

1. 固定 NCBI の `blastx_args.cpp`、`blastx_options.cpp`、`blast_query_info.c`、
   翻訳・SEG/mask・入力バッチ処理と caller を追い、option 抽出順と各 context の
   状態を記録する。Rust の移植箇所には直上の NCBI 断片・パス・行番号を記す。
2. `LOSAT blastx`、`-task blastx`、help、query/subject FASTA、stdout/`-out`、
   既定値と有効値、未対応設定の明示拒否を実装する。parser 受理と検索対応を
   混同しない。公開検索は C〜E 完了まで未実装経路を明示的に拒否する。
3. 六つの翻訳 context、26 query genetic codes、strand、NCBISTDAA/sentinel、
   stop・曖昧・端数 codon、元の核酸長と context offset を NCBI に合わせる。
   code 32 は BLASTX では拒否し、TBLASTN 限定の製品判断を広げない。
4. SEG、小文字、soft/hard mask の核酸区間から frame への変換と、
   検索用・非マスク配列の保持を実装する。単一/複数 record と batch 境界で
   S01〜S05/S11、26×64 codon の内部状態を比較する。

成果物は解決済み options、context/mask の比較、無効入力の拒否例、
source 行参照を持つ unit/integration fixture と再実行可能な証拠。
option、配列 bytes、context、mask に未説明差がないことを B gate とする。

## 終了・引き継ぎ

合格・未実行・初回差分を分け、残件があれば修正してから
[Session C](session_c_preliminary_search.md) へ渡す。
セッションで許可された commit/push の範囲に従い同ブランチへ反映する。
