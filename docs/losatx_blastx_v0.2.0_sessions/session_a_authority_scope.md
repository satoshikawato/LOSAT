# Session A — 権威・範囲・fixture

## INSTRUCTION PROMPT

LOSATX v0.2.0 の段階 A を **`feature/losatx-blastx-v0.2.0`** で実行する。
開始時にブランチと作業ツリーを確認し、`AGENTS.md`、
`verify-ncbi-parity-and-speed`、[総合計画書](../losatx_blastx_v0.2.0_plan.md)の
第 1〜8・12 節を読む。今回の段階は権威と再実行可能な基準の固定である。
検索エンジンの実装完了を主張しない。

1. TLOSAN 実装後の基点 SHA、NCBI ソース commit と該当ファイル checksum、
   BLAST+ `blastx` の版・バイナリ hash・build provenance、target と環境を固定する。
   計画中の行番号と既定値を現行の固定 C/C++ source と CLI で再確認する。
2. NCBI `blastx -help` の全キーを抽出し、X01〜X18 の成功条件、明示拒否、
   help/version、非該当に漏れなく分類する。値域、相互排他、暗黙 default と
   明示 default の差を `scope.tsv` と `defaults.json` に記す。
3. CLI から formatter まで caller/callee の実行順を追い、NCBI 関数・行範囲、
   Rust owner、入力単位・座標・精度・比較器を `source_map.tsv` に記す。
   対応ソースが見つからない能力は推測で実装せず要求根拠を再検討する。
4. 総合計画書の S01〜S15、NCBI unit fixture、R01〜R08 の取得・生成手順、
   accession/version、元 bytes と SHA、期待出力の oracle を固定する。
   初期の必須 subset と、後段で取得する大型 fixture を区別する。
5. 基点にある BLASTN/megablast、BLASTP、TBLASTX、TBLASTN の
   適用可能な現行 gate を確認し、共有部変更前の baseline を記録する。

成果物は `docs/evidence/losatx_stage_a/` の scope、source map、oracle、
fixtures、defaults、coverage、NCBI raw expected、再実行コマンドと baseline。
全必須要求に source owner と test ID があり、help に未分類キーがなく、
入力・期待値を再生成できる場合に段階 A を完了とする。

## 終了・引き継ぎ

未解決の基準や fixture を残したまま B へ進めない。結果に初回差分、
通過・未実行 gate、残件、次の具体的作業、証拠 path を記録する。
セッションで許可された commit/push の範囲に従い同ブランチへ反映し、
次は [Session B](session_b_cli_translation.md) を渡す。
