# TLOSAN（TBLASTN）v0.2.0 セッション用 INSTRUCTION PROMPTS

[総合計画書](../tlosan_tblastn_v0.2.0_plan.md)から実行指示を分離した。
以下の各ファイルは、対応する一つの作業セッションに渡す独立した指示文である。
A から G の順に進め、各セッションの引き継ぎ記録と未解決事項を次へ渡す。

作業ブランチは `feature/tlosan-tblastn-v0.2.0`。
このブランチは、依頼者が `origin/dev` の代わりに指定した
`origin/main` の 2026-09-23 時点の最新 commit
`15ea61c101ca401a6742896f2cb356ed60406350` から作成された。
各セッションの冒頭でブランチと作業ツリーを確認する。

| 順序 | セッション | 主な完了条件 |
| --- | --- | --- |
| A | [権威・基準値・製品判断](session_a_authority.md) | 再実行可能な基準、NCBI ソース対応表、遺伝暗号の境界記録 |
| B | [CLI・遺伝暗号](session_b_cli_gencode.md) | 27 ID と各 64 コドン、無効 ID の拒否 |
| C | [探索本体](session_c_search_engine.md) | 六フレームの候補・HSP 段階の一致 |
| D | [統計・連結](session_d_statistics_linking.md) | スコア、統計、連結、削除順の一致 |
| E | [表示](session_e_reporting.md) | outfmt 0/6/7 の宣言済み出力契約 |
| F | [並列・Wasm](session_f_parallel_wasm.md) | thread 数と適用可能な target 間の同一出力 |
| G | [認証](session_g_certification.md) | 全マトリクス、回帰、性能、独立監査 |

実装中はリポジトリの `AGENTS.md` と
`verify-ncbi-parity-and-speed` を適用する。NCBI ソースを唯一の動作上の権威とし、
NCBI 実行ファイルは比較オラクルにだけ使う。未実装の経路は明示的に拒否する。
