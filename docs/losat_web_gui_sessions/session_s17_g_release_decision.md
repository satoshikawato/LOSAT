# Session S17 — G：公開判定

## INSTRUCTION PROMPT

LOSAT Web の公開判定の材料をそろえる。公開するかどうかは保守者が決める（README の規則 7）。先に [セッション README](README.md) の共通規則を読み、それに従う。完了条件の正本は、総合計画書 §7 の S17 の行である。

1. 前提を確かめる：`PD-LOSAT-WEB-APP-BOUNDARY` が Accepted であること。BLASTX の v0.2.0 の認証が LOSATX 計画（`docs/losatx_blastx_v0.2.0_plan.md`）で完了し、SX が終わっていること。欠けているものがあれば、判定の記録の最初に書く。
2. 最終の commit で、5 program の認証のゲート（Gate A、TLOSAN の Stage G、LOSATX の比較、BLASTP・BLASTN・TBLASTX の既存の比較、このブランチで足した outfmt 0/7 と範囲指定の比較）と、V-ABI・V-BR を再実行する。
3. 受入表：`docs/web/requirements_trace.tsv` の initial の要求ごとに、満たしたことを示す試験と証拠のファイルを対応付ける。満たしていない要求があれば、公開できない理由として書く（設計書 §2.2「初期必須を省いたものを正式公開しない」）。
4. 対応環境の表：OS × ブラウザ（Chromium 系、Firefox、Safari / WebKit）× 経路（threaded / serial）ごとに、V-BR の結果と実機の記録（V-MOB：iOS Safari、Android Chrome）を書く。
5. 実測した上限値（入力の大きさ、メモリ、端末ごと）と性能（V-PERF：1 回の暖機と 3 回の計測、中央値と範囲）を記録し、公開する文面の案を作る。上限値を推測で書かない。
6. 既知の例外（承認済みの遺伝暗号の例外、認証済みプロファイルの外の設定、native-equivalent の升目）の一覧を、検証バッジの判定表と照らして作る。
7. 独立監査（エンジン側の主張）と画面レビュー（公開する画面）を受ける。
8. 版の名前の案と、ReleaseManifest の内容を記録する。

成果物は `docs/evidence/losat_web_g/README.md`（判定の材料と、保守者の判断を記録する欄）。

## 終了・引き継ぎ

README の規則 8 に従う。公開・タグ付け・デプロイは、保守者の判断の後に別の作業として行う。
