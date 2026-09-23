# Session B — CLI・遺伝暗号

## INSTRUCTION PROMPT

`feature/tlosan-tblastn-v0.2.0` で TLOSAN v0.2.0 の段階 B を実行する。
開始時にブランチと作業ツリー、[総合計画書](../tlosan_tblastn_v0.2.0_plan.md)、
`AGENTS.md`、`verify-ncbi-parity-and-speed`、段階 A の証拠を確認する。
段階 A の境界記録が未完成なら、その不足を先に解消する。

1. NCBI の `tblastn_args.cpp`、`blast_args.cpp`、
   `tblastn_options.cpp` を現行 checkout で追い、
   `tblastn` の引数、既定値、禁止組合せとエラーを Rust CLI に移植する。
   `tblastn-fast`、PSI-TBLASTN、DB 検索など未実装経路を明示的に拒否する。
2. NCBI `gc.prt` の 27 ID を明示的に検証する。TBLASTN CLI は
   段階 A の製品判断に従い 32 も受理し、不正 ID をコード 1 に置換しない。
   `utils/genetic_code.rs` と `core/gencode_singleton.rs` の重複表・
   暗黙フォールバックを調べ、NCBI と一致する単一の動作に揃える。
   既存 TBLASTX の契約を壊さない。
3. 各 ID の 64 コドン、曖昧塩基、終止と不正 ID を単位検証する。
   コード 32 は段階 A の比較専用 API オラクルを使う。
   Rust の移植箇所には `AGENTS.md` 所定の NCBI ファイル・行番号と
   C/C++ 断片を直上に記す。
4. 既存の TBLASTX 遺伝暗号回帰を実行する。探索・表示が未完成の間は
   TBLASTN の部分的成功を一般利用可能な機能として公開しない。

成果物は CLI 受理・拒否表、27 × 64 の照合結果、不正 ID のエラー、
既存回帰の結果、変更箇所の NCBI 対応表である。段階 C へ残す未対応を明記する。
