# Session C — 探索本体

## INSTRUCTION PROMPT

`feature/tlosan-tblastn-v0.2.0` で TLOSAN v0.2.0 の段階 C を実行する。
開始時にブランチと作業ツリー、総合計画書 `docs/tlosan_tblastn_v0.2.0_plan.md`、
`AGENTS.md`、`verify-ncbi-parity-and-speed`、段階 A・B の証拠を確認する。
段階 B の遺伝暗号検証が未完なら先に解消する。

1. NCBI `blast_engine.c:747-888` を起点に、TBLASTN のタンパク質 query 符号化、
   lookup、subject 六フレーム翻訳、各フレーム探索、HSP 統合の実際の呼び出し順を
   現行 checkout で追う。query/subject の座標系、センチネル、マスク、
   フレーム順、入力状態を記録する。
2. 対応する LOSAT の BLASTP/TBLASTX 経路を調べ、NCBI の引数と状態が等しい
   演算だけを共有する。異なるタイミングや入力状態の関数をそのまま流用しない。
   六フレームの初期 HSP、統合、再評価、削除を NCBI と同じ順序で実装する。
3. +1/+2/+3/−1/−2/−3、末端の部分コドン、曖昧塩基、終止、
   低複雑度、no-hit、同点を含む小型ケースを固定する。
   新規生成した NCBI 出力と段階別トレースを比較し、最初の相違を
   候補、raw score、内部座標、順序へ分解して解消する。
4. Rust の移植箇所には NCBI ファイル・行番号と該当 C/C++ 断片を直上に記す。
   後段の統計・表示が未完成なら、公開 CLI は明示的な未実装エラーを維持する。

成果物は NCBI/LOSAT 呼び出し経路対応表、六フレームの固定 fixture、
段階別差分、変更箇所のソースコメント、候補・HSP 段階の一致である。
段階 D が扱う統計・連結の未解決を区別して引き継ぐ。

## 終了・引き継ぎ

段階の完了条件を満たすまで必要な実装・検証・修正を続ける。
作業終了時に変更と再実行可能な証拠を `feature/tlosan-tblastn-v0.2.0` へ
コミットし、`origin/feature/tlosan-tblastn-v0.2.0` にプッシュする。
最終回答にコミット SHA、プッシュ結果、検証結果、残件を記し、
次セッション（D — 統計・連結）でそのまま使える Codex 用
INSTRUCTION PROMPT を全文で提示する。
段階 C の完了条件が残る場合は、その解消を次プロンプトの最初の作業とし、
完了や認証を宣言しない。
