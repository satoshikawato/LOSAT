# Session G — 総合認証・性能

## INSTRUCTION PROMPT

LOSATX v0.2.0 の段階 G を **`feature/losatx-blastx-v0.2.0`** で実行する。
ブランチ・作業ツリー、`AGENTS.md`、`verify-ncbi-parity-and-speed`、
[総合計画書](../losatx_blastx_v0.2.0_plan.md)の第 7〜14 節、段階 A〜F の
証拠を確認する。未解決の必須差分があれば原因を解決してから認証する。

1. coverage ledger を X01〜X18→NCBI owner→Rust owner→fixture→evidence で
   検査する。M1 の 26 codes×0/6/7×threads {1,4}=156 unique cases、
   M2〜M11 の必須境界を新しい実行場所で再現する。missing、SKIP、
   重複計上、未知 fingerprint、raw byte 差を PASS にしない。
2. Native 配布 platform、serial/threaded command WASI、serial/threaded
   reactor と実ブラウザの指定入口を同じ options で確認する。
   build/test/fmt/clippy と関連回帰を実施する。基点の BLASTN/megablast、
   BLASTP、TBLASTX、TBLASTN を同じ候補 SHA で回帰検証する。
3. パリティ確立後、第 11 節の固定 workload で各 mode を untimed warmup 1 回、
   timed 3 回ちょうど測る。中央値、三標本の min–max、全 raw samples、
   output hash を保存する。NCBI の時間は事前 `makeblastdb -dbtype prot` の
   `-db` 検索で測り、DB 作成の版・command・入力 hash・経過時間を別に残す。
   NCBI local `-subject` の分布 oracle と DB timing output を混同しない。
4. source/input/binary identity、全コマンド、初回差分と修正、例外範囲、
   未対応機能を evidence template に記す。発売向けパリティ・性能主張の前に
   `ncbi_parity_auditor` の独立した読み取り専用監査を依頼し、指摘を解消する。
5. README、CHANGELOG、scope、manifest を証明済みの候補 SHA と能力に揃える。
   明示された製品範囲外を実装済みと表示しない。

X01〜X18 と第 14 節の全項目を満たした場合だけ v0.2.0 gate を合格とする。
有限 fixture の一致を全入力の証明と表現しない。

## 終了・引き継ぎ

認証記録に candidate SHA、全 gate、監査結果、残件と適用範囲を残す。
セッションで許可された commit/push の範囲に従い同ブランチへ反映する。
タグ、release 公開、package publication は別の明示的な許可境界とする。
