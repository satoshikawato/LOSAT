# Session F — 並列・WASI・Web/reactor

## INSTRUCTION PROMPT

LOSATX v0.2.0 の段階 F を **`feature/losatx-blastx-v0.2.0`** で実行する。
ブランチ・作業ツリー、`AGENTS.md`、`verify-ncbi-parity-and-speed`、
[総合計画書](../losatx_blastx_v0.2.0_plan.md)の X15〜X17・第 8/10 節、
段階 A〜E の証拠を確認する。Native serial の必須 raw 差を先に解消する。

1. 既存 search-scoped pool で subject 計算を分配し、入力 OID 順に回収する。
   collector、heap、Kappa redo と同点比較を NCBI 順で処理する。
   並列化する order-sensitive 処理には byte 同値の独立証拠を要する。
2. Native threads 1/2/4/8 と worker 数の前後の subject 数、同点、
   多数 query、26 codes、0/6/7、反復実行で M8 を照合し、実 worker 稼働を示す。
3. serial `wasm32-wasip1` command と、`wasm32-wasip1-threads` +
   `wasm-threads` command をそれぞれ実行し、Native と同じ bytes を確認する。
   serial target を threaded と表示しない。
4. `losat_web_run_pair` と `_handles` に同じ Rust BLASTX engine を接続する。
   serial/threaded reactor、in-memory/handle、0/6/7/custom、失敗後回復、
   invalid handle、pool failure、result clearing、worker/memory 再利用を検証する。
   Node reactor と実ブラウザ smoke を区別する。
5. NCBI owner に対応する共有処理の変更にはソース断片・行番号を付け、
   変更した BLASTN/BLASTP/TBLASTX/TBLASTN caller の該当回帰を実行する。
   offline bundle/依存関係を変えた場合のみ browser-offline-qa を適用する。

成果物は M8/M9 の全 command/target/binary hash、出力比較、worker lifecycle
証拠。serial/threaded と各入口の bytes 一致、実並列、エラー後回復を
確認した場合に F 完了。build 成功だけでは完了しない。

## 終了・引き継ぎ

未説明の byte 差や lifecycle failure を解決し、通過・未実行を区分する。
セッションで許可された commit/push の範囲に従い同ブランチへ反映して
[Session G](session_g_certification.md) へ渡す。
