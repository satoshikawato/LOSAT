# Session D — 統計・linking・Kappa・filter

## INSTRUCTION PROMPT

LOSATX v0.2.0 の段階 D を **`feature/losatx-blastx-v0.2.0`** で実行する。
ブランチ・作業ツリー、`AGENTS.md`、`verify-ncbi-parity-and-speed`、
[総合計画書](../losatx_blastx_v0.2.0_plan.md)の X07・X10〜X12 と第 4〜8 節、
段階 A〜C の証拠を確認する。C の候補/HSP 差分が残るなら先に解消する。

1. 固定 NCBI の `blast_setup.c`、`blast_stat.c`、`blast_parameters.c`、
   `link_hsps.c`、`blast_kappa.c`、`blast_traceback.c`、`blast_hits.c`、
   culling/subject-besthit の caller を追い、effective length、search space、
   sum statistics、query 側 gap/intron、link/reap の順序を記録する。
2. six-context Kappa、composition 2 の redo、行列選択・scale・丸め、
   composition 0 の通常 traceback、score/bit/E-value、identity/positive、
   post-traceback pruning を移植する。六つの独立 BLASTP の結合で代替しない。
3. max targets/HSP、culling、subject_besthit、heap、early termination を
   NCBI のタイミングと比較器で行う。TBLASTN の translated-subject 統計を
   protein subject に適用しない。Rust の変更直上に NCBI 断片・行番号を付ける。
4. S08/S09/S14、全 code の縮小 fixture と実データで最終 HSP membership、
   raw/bit score、full-precision E-value、削除理由、同点順序を照合する。
   不一致は最初の段階まで戻って原因を解決する。

成果物は実行順を含む source/Rust 対応、統計・filter trace、再実行コマンド、
最初の差分と解消記録。最終 HSP 集合・順序・数値に未説明差がなければ D 完了。
表示は未完成でも、公開 CLI は未実装境界を維持する。

## 終了・引き継ぎ

D の差分を E の表示差として移さない。証拠と残件を明記し、
セッションで許可された commit/push の範囲に従い同ブランチへ反映して
[Session E](session_e_reporting_native.md) へ渡す。
