# Session E — 表示

## INSTRUCTION PROMPT

`feature/tlosan-tblastn-v0.2.0` で TLOSAN v0.2.0 の段階 E を実行する。
開始時にブランチと作業ツリー、[総合計画書](../tlosan_tblastn_v0.2.0_plan.md)、
`AGENTS.md`、`verify-ncbi-parity-and-speed`、段階 A〜D の証拠を確認する。
探索・統計の既知差を表示の加工で隠さない。

1. NCBI `tabular.cpp`、`showalign.cpp`、`format_flags.cpp` と
   TBLASTN の formatter 呼び出し元を現行 checkout で追う。
   タンパク質 query 座標と subject 塩基座標、正負六フレーム、ID、
   アラインメント文字列、数値丸め、ヘッダー・フッターの変換順を記録する。
2. `-outfmt 0` のアラインメント本文・統計部、`6` の各列、
   `7` の各クエリ見出し・no-hit・フッターを実装する。
   既存の簡易 pairwise writer が異なる場合は NCBI 表示経路を優先する。
   Rust の移植箇所の直上に NCBI ファイル・行番号と C/C++ 断片を記す。
3. コード 1 は同じローカル `-subject` 条件で、単一・複数 query、
   六フレーム、同点、複数 HSP、no-hit、0/6/7 の生バイトを比較する。
   非標準コードは段階 A の限定契約で差分を分類し、書式・順序の
   無関係な差を許さない。コード 32 は比較専用 API オラクルも使う。
4. 差があれば最初の相違を query/subject ID、フレーム、座標、score、
   E-value、HSP 順序、表示文字列へ分解し、全該当ケースを一括修正する。

成果物は固定コマンド・入力 SHA-256・生バイト比較、非標準コードの
限定差分の説明、0/6/7 の宣言済み範囲と未対応オプションの拒否表である。
段階 F に serial の基準出力 SHA-256 を渡す。
