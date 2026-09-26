# Session C — seed・preliminary 探索

## INSTRUCTION PROMPT

LOSATX v0.2.0 の段階 C を **`feature/losatx-blastx-v0.2.0`** で実行する。
ブランチ・作業ツリー、`AGENTS.md`、`verify-ncbi-parity-and-speed`、
[総合計画書](../losatx_blastx_v0.2.0_plan.md)の X06・X08・X14 と第 4〜8 節、
段階 A/B の証拠を確認する。前段の必須差分を解決してから着手する。

1. 固定 NCBI の `blast_aalookup.c`、`blast_aascan.c`、`aa_ungapped.c`、
   `blast_engine.c`、`blast_gapalign.c`、`blast_parameters.c`、split-query の
   caller と状態を追う。word 3/5、閾値、one/two-hit window、圧縮 lookup、
   context-aware 対角管理、cutoff の時点と比較器を記録する。
2. protein subject の入力 OID 順、seed/ungapped、preliminary gapped HSP、
   common-endpoint と containment、保存・pruning の順序を NCBI に合わせる。
   `-ungapped -comp_based_stats 0` は独立の有効分岐として扱う。
3. 長 query の chunk/overlap と六 frame の再統合を移植する。
   S06/S07/S10/S12、短縮実データと NCBI unit fixture で、候補集合、
   HSP の score・内部座標・順序、split 前後を段階別に照合する。
4. Rust の各変更直上に対応 NCBI C/C++ 断片、ファイル・行番号を付ける。
   NCBI 実行ファイルを runtime/build/fallback に入れない。

成果物は seed→preliminary→split/rejoin の比較 trace と最初の差分、
fixture/command/hash を持つ証拠。選んだ必須 C fixture の内部 HSP が
score・座標・順序まで一致し、NCBI 記録値を Rust 計算の代わりに注入して
合格にしていない場合に C を完了とする。

## 終了・引き継ぎ

C の相違を残したまま統計や formatter の検証へ進めない。
結果と証拠を同ブランチに整理し、セッションで許可された commit/push の
範囲に従って反映して [Session D](session_d_statistics_linking.md) へ渡す。
