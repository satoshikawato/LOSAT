# Session D — 統計・連結

## INSTRUCTION PROMPT

`feature/tlosan-tblastn-v0.2.0` で TLOSAN v0.2.0 の段階 D を実行する。
開始時にブランチと作業ツリー、[総合計画書](../tlosan_tblastn_v0.2.0_plan.md)、
`AGENTS.md`、`verify-ncbi-parity-and-speed`、段階 A〜C の証拠を確認する。
探索段階に既知の NCBI 相違が残るなら、その境界を解決してから統計を検証する。

1. 現行 NCBI ソースの `blast_parameters.c`、`link_hsps.c`、
   `blast_kappa.c`、`blast_hits.c` と呼び出し元を追い、
   有効長、翻訳 subject の長さ換算、sum statistics、
   intron/linking、組成補正モード 2、Kappa redo、E-value、
   bit score、フィルター、hitlist 順序を実際の実行順に記録する。
   `-max_intron_length 0` を linking 無効として扱わない。
2. BLASTP の `do_link_hsps=false` を TBLASTN へ流用せず、
   NCBI の比較器・同点処理・浮動小数点精度を移植する。
   Rust の移植箇所の直上に NCBI ファイル・行番号と C/C++ 断片を記す。
3. 段階 C の fixture で最初の差分を追い、コード 1 の HSP と統計値を
   NCBI ローカルと比較する。非標準コードは段階 A の限定契約に従い、
   選択した subject 翻訳だけに起因する差かをコード 1 対照と
   NCBI `-db` 参照で説明する。HSP 数だけで合格にしない。
4. 0/6/7 の表示が未完成なら公開 CLI の未実装境界を維持する。
   未実装の統計・補正モードを暗黙に既定値へ戻さない。

成果物は関数・実行順対応表、固定 fixture の raw/bit score・E-value・
順位・削除順の比較、全残差の分類である。段階 E に表示だけの差を引き継ぐ。
