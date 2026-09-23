# Session G — 認証

## INSTRUCTION PROMPT

`feature/tlosan-tblastn-v0.2.0` で TLOSAN v0.2.0 の段階 G を実行する。
開始時にブランチと作業ツリー、[総合計画書](../tlosan_tblastn_v0.2.0_plan.md)、
`AGENTS.md`、`verify-ncbi-parity-and-speed`、段階 A〜F の証拠を確認する。
未解決の NCBI 相違があれば認証済みという主張をせず、原因を解決して再検証する。

1. マニフェストの全 27 ID × 0/6/7 × threads {1,4} の 162 条件を
   新しい一時ディレクトリで実行する。コード 1 は NCBI ローカル
   `-subject` と生バイト一致を求める。非標準コードは明示された
   subject genetic code の限定契約で差を説明し、それ以外の差を拒否する。
   コード 32 は比較専用 NCBI API オラクルとの HSP・統計・表示照合を必須とする。
2. 代表的なコード 1/4/11 と実データで threads 2/8、複数 subject、
   同点、no-hit を追加する。Native と認証対象 Wasm の出力を比較し、
   既存 BLASTN/BLASTP/TBLASTX 回帰、build、test、fmt、clippy、
   適用可能な Wasm ゲートを実施する。
3. 性能は固定 release build、入力、出力先、thread 数を使い、
   untimed warmup 1 回と timed 3 回の中央値・最小最大を報告する。
   NCBI の速度測定は事前 `makeblastdb` の `-db` 検索を使い、
   DB 作成コマンド・時間・版・入力 SHA-256 を検索時間と別に保持する。
   速度測定出力と分布比較出力を混同しない。
4. `verify-ncbi-parity-and-speed/references/evidence.md` の各欄を埋め、
   入力・出力 checksum、ソース参照、最初の差分、例外、全コマンド、
   生サンプル、未対応範囲を残す。発売向けのパリティ・性能主張は
   `ncbi_parity_auditor` に独立した読み取り専用監査を依頼する。

全条件と監査が揃った場合だけ v0.2.0 の認証結果を記す。
満たせない条件は対象外または未解決として根拠とともに記録し、
部分的な一致を完全パリティと表現しない。
