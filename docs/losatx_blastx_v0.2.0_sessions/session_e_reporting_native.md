# Session E — 表示・Native CLI

## INSTRUCTION PROMPT

LOSATX v0.2.0 の段階 E を **`feature/losatx-blastx-v0.2.0`** で実行する。
ブランチ・作業ツリー、`AGENTS.md`、`verify-ncbi-parity-and-speed`、
[総合計画書](../losatx_blastx_v0.2.0_plan.md)の X01〜X03・X13・第 3.3/4/7/8 節、
段階 A〜D の証拠を確認する。内部 HSP と統計の差は E 着手前に解消する。

1. NCBI `blast_seqalign.cpp` と formatter の実際の BLASTX caller を追い、
   query 核酸/subject アミノ酸座標、正負 frame、alignment・BTOP、
   ID/title、丸め、no-hit、prolog/footer を固定する。Rust 変更直上に
   NCBI C/C++ 断片・パス・行番号を付ける。
2. 既定 outfmt 0、標準 6/7、総合計画書の全 30 custom fields と `std`、
   列順・重複・未知列拒否を実装する。省略/明示 default、複数 query、
   batch reset、stdout/`-out`、stderr/exit を同じ解決済み options で照合する。
3. Native serial の M1〜M7/M10 を固定 NCBI local `-subject` の raw bytes と
   比較する。26 code×0/6/7、six frames、mask、limits/filter、長 query、
   no-hit と無効入力を省かない。format 6 の数値だけで完了としない。
4. 全必須検索経路が通る時点で公開 `LOSAT blastx` を接続し、help と文書を
   実測した能力に揃える。未対応 option は暗黙の default に変換しない。

成果物は fixture ごとの入力/期待/実出力 hash、最初の byte diff、
全 native serial coverage と CLI/error 記録。必須 M1〜M7/M10 の
未説明 raw byte 差がゼロなら E 完了。

## 終了・引き継ぎ

差があれば対応する NCBI owner に戻って修正する。結果と未実行 gate を記し、
セッションで許可された commit/push の範囲に従い同ブランチへ反映して
[Session F](session_f_parallel_wasi_web.md) へ渡す。
