# LOSATX（BLASTX）v0.2.0 セッション用 INSTRUCTION PROMPTS

[総合計画書](../losatx_blastx_v0.2.0_plan.md)にある契約と合格条件を、
セッション A〜G で実行するための指示文に分けた。各ファイルの
「INSTRUCTION PROMPT」以降を次の作業セッションへ渡す。

実装ブランチは **`feature/losatx-blastx-v0.2.0`**。
TLOSAN 実装後の `origin/feature/tlosan-tblastn-v0.2.0` の
`2976bd5f427cc5448315b6605ab85ee0787545a5` から作成された。
各セッションでは開始時にブランチと作業ツリーを確認し、このブランチ上で作業する。
別のブランチへ黙って切り替えない。

| 順序 | 指示文 | 完了条件の要点 |
| --- | --- | --- |
| A | [権威・範囲・fixture](session_a_authority_scope.md) | NCBI 分岐と全公開キー、入力・oracle・既存基準を固定 |
| B | [CLI・翻訳・query 準備](session_b_cli_translation.md) | option、六 context、26 code、mask と入力境界を照合 |
| C | [seed・preliminary 探索](session_c_preliminary_search.md) | word 3/5、gapped/ungapped、split 前後の HSP を照合 |
| D | [統計・linking・Kappa・filter](session_d_statistics_linking.md) | 最終 HSP、数値、順序、削除理由を照合 |
| E | [表示・Native CLI](session_e_reporting_native.md) | outfmt 0/6/7、30 列と Native serial の raw parity |
| F | [並列・WASI・Web](session_f_parallel_wasi_web.md) | threads と入口をまたいで同一 bytes、worker lifecycle |
| G | [総合認証・性能](session_g_certification.md) | 全 coverage、回帰、性能標本、独立監査 |

各セッションで [AGENTS.md](../../AGENTS.md) と
[verify-ncbi-parity-and-speed](../../.agents/skills/verify-ncbi-parity-and-speed/SKILL.md)
を適用する。NCBI C/C++ は唯一の挙動権威、実行ファイルは比較専用 oracle。
Rust の変更には直上に該当 NCBI 断片・ファイル・行番号を付ける。
未完成の公開機能は明示的に拒否し、raw byte 差を許容差で合格にしない。

各段階の証拠には入力・実行ファイル・ソースの identity、全コマンド、
最初の差分、合格・未実行・未解決の区分、次段階への引継ぎを残す。
後続段階の着手前に前段の未解決必須差分を解消する。
