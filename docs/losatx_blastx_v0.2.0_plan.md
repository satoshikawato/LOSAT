# LOSATX（BLASTX）v0.2.0 総合実装・互換性検証計画

作成日: 2026-09-26。状態: **計画。実装・認証は未実施**。

**実装ブランチ: `feature/losatx-blastx-v0.2.0`。** 2026-09-26 に取得した
TLOSAN 実装後の `origin/feature/tlosan-tblastn-v0.2.0` の `2976bd5f427cc5448315b6605ab85ee0787545a5` を基点とする。
各実装セッションは
このブランチを使用し、開始時に `git branch --show-current` と `git status` で確認する。
段階別の実行指示は[セッション用 INSTRUCTION PROMPTS](losatx_blastx_v0.2.0_sessions/README.md)に分離した。

LOSATX は LOSAT による NCBI BLASTX の純 Rust 実装を指す。公開検索コマンドは
`LOSAT blastx` とする。本書はユーザーの「v0.2.0 では LOSATN、LOSATP、
TLOSATX と同等の NCBI 互換性を必須とする」という要件を実装・検証可能な
条件に分解したものである。既定設定だけの試作版を v0.2.0 完了とは扱わない。

## 1. 目的と完了の定義

核酸 query とタンパク質 subject のローカル FASTA 検索を実装する。
query の六フレーム翻訳、アミノ酸検索、ギャップ伸長、組成補正、HSP 連結、
統計、枝刈り、順位付け、出力を NCBI と同じタイミング・順序・入力状態で行う。

「同等の互換性」は次のすべてを満たすことと定義する。

1. **出力精度**: 宣言した対応範囲では、NCBI BLAST+ の出力を行順・書式を含め
   生バイトで再現する。ヒット数、座標、スコアだけの一致では完了しない。
2. **機能範囲**: 既存三プログラムで実装されている機能のうち、NCBI BLASTX
   に対応する機能を第 3 節の必須表に取り込む。LOSATP のカスタム列や
   TLOSATX の culling を、既定 outfmt 6 の認証だけで省略しない。
3. **実行面**: Native、serial WASI、threaded WASI、既存 Web/reactor の
   BLASTX 入口で同じ検索意味論を提供する。threads 数や繰り返しで結果を変えない。
4. **入力・失敗面**: 複数 query/subject、no-hit、短い入力、警告、無効入力、
   未対応設定の拒否を仕様化し、成功例だけを認証しない。
5. **既存機能の維持**: BLASTN/megablast、BLASTP、TBLASTX、および同じツリーの
   TBLASTN の既存契約を維持する。共通部品の変更は該当回帰で検証する。
6. **実装独立性**: NCBI バイナリ・ライブラリ・FFI・subprocess を LOSAT の
   runtime、build、fallback、未対応機能の実装に一切使わない。

互換性の対象は既存 LOSAT と同じローカル比較製品である。NCBI BLAST+ 全製品の
全タスク、DB 形式、remote 検索、全出力形式への対応を意味しない。
一方、第 3 節で必須とした機能を「未対応として拒否できる」だけでは合格しない。
必須機能に未実装・不一致が残れば **LOSATX の v0.2.0 gate は不合格**とする。

網羅性は、要求・NCBI 分岐・入力境界・相互作用・ターゲットを追跡することで
示す。有限の fixture 合格を「すべての入力で証明済み」とは表現しない。

## 2. 調査基点と権威

### 2.1 初期調査時の状態と実装基点

| 項目 | 確認した状態 |
| --- | --- |
| 初期調査時の LOSAT HEAD | `5b9f6b9b3002390261aa8a8cafd2a6e9f6866200` |
| 初期調査時の作業ブランチ | `feature/tlosan-tblastn-v0.2.0` |
| 作業ツリー | 既存の未コミット変更あり。HEAD だけを調査したソース全体の識別子にしない |
| NCBI ソース | `/mnt/c/Users/genom/GitHub/ncbi-blast/` |
| NCBI ソース commit | `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4` |
| 比較実行ファイル候補 | `/home/kawato/micromamba/bin/blastx` |
| その版 | `blastx: 2.17.0+`、build `Aug 11 2025 09:46:06` |
| その SHA-256 | `5c29dd472f7b18db37b688077ea4a48eea8f46ccfbeab64e11ca5efddc142c05` |

この表は元の計画を作成した際の調査記録であり、リリース候補 SHA の指定ではない。
実装基点は TLOSAN 実装後の `2976bd5f427cc5448315b6605ab85ee0787545a5`、作業先は
`feature/losatx-blastx-v0.2.0` である。段階 A で関連差分、ソース・入力・
実行ファイルの checksum を再固定する。TLOSAN 側の過去の合格記録を
BLASTX の合格記録へ転用しない。

### 2.2 優先する資料

- [AGENTS.md](../AGENTS.md): 現行の最上位リポジトリ規則。
- [verify-ncbi-parity-and-speed](../.agents/skills/verify-ncbi-parity-and-speed/SKILL.md):
  差分凍結、NCBI owner 追跡、比較、性能、証拠記録。
- [CLI v2](cli_v2_migration.md): single-dash/snake_case、公開能力と内部設定の区別。
- 現行 `LOSAT/src/algorithm/{blastn,blastp,tblastx}/args.rs` と検索入口の validation。
- [TLOSAN 計画](tlosan_tblastn_v0.2.0_plan.md)と
  [段階 G 記録](evidence/tlosan_stage_g/STAGE_G_CERTIFICATION_20260926.md):
  再利用候補、バッチ・統計・出力の検証方法。BLASTX の証拠への転用は不可。
- [NCBI CLI マニュアル](https://www.ncbi.nlm.nih.gov/books/NBK279684/): 引数表の補助。
  挙動の最終権威は固定 C/C++ ソースと対応する executable oracle。

古い計画や README の概要と現行 validation が異なる場合、公開引数、検索の
実対応、認証済み範囲を分けて記録する。過去の未解決件数や速度を現状と断定しない。

### 2.3 既存三プログラムから継承する下限

| 比較対象 | 現行コードで確認した境界 | LOSATX への対応 |
| --- | --- | --- |
| LOSATN | local FASTA、megablast/blastn、標準 6/7、HSP/target 制限、小文字マスク、subject_besthit | local 入力、制限、マスク、該当するフィルター、標準 tabular の精度を必須化 |
| LOSATP | ordinary blastp、実行時は BLOSUM62・gap 11/1・word 3/5・composition 2、0/6/7 とカスタム列 | 同等以上の蛋白検索構成、word 5 の圧縮 lookup、0/6/7、全 30 カスタム列を必須化 |
| TLOSATX | query genetic code 26 ID、word 3、SEG、one/two-hit window、culling、標準 6 | query 六フレーム、同じコード受付、mask/window/culling を必須化 |
| 共通基盤 | single-dash CLI、Native/WASI、Web/reactor、thread pool、比較 manifest | 同じ入口・実行モードへの統合と決定的な出力を必須化 |

BLASTP の parser が受理しても検索時に拒否する別 matrix/gap、composition 1/3、
word 2/4/6/7、Smith-Waterman は「既存実対応」と数えない。
逆に古い v0.1.0 認証が outfmt 6 のみであっても、現行 LOSATP の 0/7 と
カスタム列を LOSATX の計画から削らない。
BLASTN 専用 reward/penalty、DUST、核酸 lookup 等は BLASTX へ移植しない。

現行 Web の `parse_blastp_args` は native CLI の task 制限を通らず、内部の
`blastp-fast` 構成へ到達できる。しかし [CLI v2 の記録](cli_v2_migration.md)は
その historical fast-default に retained-row 差があることを明記しており、
公開 CLI では task を `blastp` に制限している。この入口間の不整合を記録し、
未認証の内部到達性を本計画の「対応済み」下限に数えない。LOSATX の Web 入口は
本書の対応表と同じ validation を通し、未対応 task を迂回受理しない。

## 3. v0.2.0 必須対応表

### 3.1 公開機能

| ID | 必須能力 | 具体的範囲・成功条件 |
| --- | --- | --- |
| X01 | CLI | `LOSAT blastx`、`-task blastx`、CLI v2 の canonical 表記、help、既存 root version 規約 |
| X02 | ローカル入力 | 必須 `-query` 核酸 FASTA、`-subject` 蛋白 FASTA。単一・複数 record、順序、ID・説明文、LF/CRLF、終端改行なし |
| X03 | 出力先 | stdout と `-out`。同じ report bytes、I/O エラーを成功扱いしない |
| X04 | 翻訳 | `-query_gencode` 26 ID、`-strand both/plus/minus`、正負六フレーム、曖昧コドン、stop、不完全末端コドン |
| X05 | スコア | BLOSUM62、gapped gapopen 11/gapextend 1、省略・明示同値。異なる有効コストへの暗黙丸め禁止 |
| X06 | lookup/seed | `-word_size 3/5`、`-threshold`、`-window_size`。0 の one-hit と正値の two-hit。word 5 の NCBI 圧縮 lookup |
| X07 | 組成補正 | 既定 2 と別名 D/d/T/t、対照 0/F/f。six-context Kappa、再アラインメント、正規化、削除を検索経路として実装 |
| X08 | ungapped | `-ungapped -comp_based_stats 0`。composition 有効との不正組合せを NCBI に従い拒否。既定 gapped へのすり替え禁止 |
| X09 | マスク | `-seg no/yes/"window locut hicut"`、`-soft_masking true/false`、`-lcase_masking`。核酸区間→翻訳 frame の変換と非マスク配列保持 |
| X10 | 統計・連結 | `-evalue`、`-sum_stats true/false`、`-max_intron_length` の NCBI CLI 有効範囲。有効長、search space、gap 条件、連結・再評価の時点 |
| X11 | 結果制限 | `-max_target_seqs`、`-max_hsps`。同点・上限境界、複数 frame の統合後の制限、query ごとの hitlist |
| X12 | フィルター | `-culling_limit`、`-subject_besthit`。NCBI BLASTX の context/座標/タイミングで適用 |
| X13 | 表示 | 既定 `-outfmt 0`、標準 6/7、6/7 の第 3.3 節の全カスタム列 |
| X14 | バッチ・長配列 | 複数 query のバッチ、長い単一 query の分割・再統合、フレーム・重複領域・出力順の維持 |
| X15 | 並列 | Native threads 1/2/4/8、subject 数が worker 数より少ない・多い場合、同一プロセスでの繰り返し |
| X16 | Wasm/Web | serial/threaded command WASI、既存 serial/threaded reactor と in-memory/handle 入力の BLASTX 統合 |
| X17 | エラー・診断 | no-hit と error の区別、警告、未対応・無効設定の明示拒否、古い結果や半端な成功結果を返さない |
| X18 | 回帰・証拠 | 既存四プログラムの適用可能な gate、機械可読 coverage、入力・出力 SHA、独立監査 |

X07 の mode 0、X08、BLASTX 固有の strand/sum_stats/intron は、本計画で
NCBI BLASTX の動作を直接検証するための必須追加範囲とする。
これらを既存 BLASTP で実装済みと説明しない。
値域の厳密な上限、相互排他、明示指定時の派生 default は段階 A/B で
固定ソースから登録する。NCBI が無効とする組合せは成功 matrix に含めない。

### 3.2 既定値と省略時動作

| 項目 | 基準 |
| --- | --- |
| task / strand / genetic code | `blastx` / `both` / `1` |
| matrix / gaps | BLOSUM62 / 11,1 |
| word size / threshold / window | 通常 word 3 / 12 / 40。word 5 明示時の派生値を別登録 |
| E-value / max targets / threads | 10 / 500 / 1 |
| composition / gapped | 2 / gapped |
| outfmt | 0。省略で pairwise 出力を実行可能にする |
| soft masking / lowercase | false / 明示フラグなし |
| SEG | CLI help は `12 2.2 2.5`。API constructor の `SetSegFiltering(false)` と区別し、CLI 抽出後の状態・無指定出力を段階 A で固定 |
| sum statistics | `blastx_options.cpp:86-90` の設定と CLI 抽出後の状態を固定 |
| max intron | CLI default 0。0 の意味は `blast_parameters.c` と linking の実行から確定し、無効化と推測しない |
| max HSPs / culling / subject best hit | 未指定上限なし / 0 / false |

NCBI の API constructor、CLI help、CLI の option 抽出は別の層である。
既存 BLASTP/TBLASTN の default を丸ごとコピーしない。
少なくとも「すべて省略」「各値を明示」「全基準値を明示」の三者で、
同値になるべきケースと NCBI 自身が差を生むケースを登録する。
pairwise の description/alignment 件数は formatter の default を確認し、
`max_target_seqs` の一律な後段切り詰めで代替しない。

### 3.3 カスタム tabular 出力

現行 LOSATP の parser にある次の **30 フィールド**と `std` を必須とする。

```text
qseqid qacc qaccver qlen sseqid sacc saccver slen
qstart qend sstart send qseq sseq
evalue bitscore score length pident nident mismatch positive gapopen gaps ppos
qframe sframe frames btop stitle
```

標準 12 列の指定は固定 NCBI の `std`、すなわち
`qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore`
を参照する。単純な FASTA で同じ文字列になることを理由に、qseqid と qaccver
の区別を消さない。指定順、`std` 展開、重複指定、未知列、空 specification、
outfmt 0 に列を付けた場合を検証する。

BLASTX の qstart/qend/qlen は核酸側、subject 座標・長さはアミノ酸側として
NCBI formatter の単位を追う。alignment length、qseq/sseq、BTOP の単位を
核酸座標から推定しない。qframe は六フレーム、蛋白 subject の sframe と
frames の表記は NCBI が生成した期待値を使う。正鎖用の換算式を負鎖へ流用しない。

### 3.4 遺伝暗号の契約

必須受付 ID は、既存 TLOSATX の query-gencode と固定 BLASTX CLI に揃える。

```text
1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,33
```

- 全 26 ID で翻訳・探索・組成補正・統計・座標・0/6/7 まで通す。
- 無効 ID を 1 に置換しない。コードごとの差が現れるコドンを入力に入れる。
- `db_gencode` は蛋白 subject を持つ BLASTX の公開オプションにしない。
- TBLASTX/TBLASTN のローカル subject genetic-code 例外は BLASTX に適用しない。
- **コード 32 は現契約では拒否する。** 全 27 ID への拡張は別の明示的な製品判断が
  必要であり、TBLASTN 限定 `PD-TLOSAN-LOCAL-GENCODE-32` から許可を推論しない。
  拡張時は query code 32 を通した比較専用 NCBI C++ API の検索・format oracle と、
  code 1 の同入力 CLI 校正を用意する。翻訳テーブル単体の一致だけでは認証しない。

コード 32 の追加判断がなくても、本書の 26 ID 契約を完了できる。
共有 genetic-code 表は再利用するが、program ごとの CLI 許可集合は分離する。

### 3.5 今回の製品範囲外と明示拒否

次は既存三プログラムの実対応を下回らない範囲外項目として登録する。
無視、近い値への置換、外部 BLAST 呼出しを禁止する。

| 分類 | 対象 |
| --- | --- |
| 別 task | `blastx-fast`、NCBI BLASTX に存在しない独自 task |
| DB / remote | `-db`、`-remote`、GI/seqid/taxid/IPG 制限、Entrez、DB mask、taxonomy DB 依存 |
| 未対応 scoring | BLOSUM62 以外、gapped 11/1 以外、word 2/4/6/7 等、composition 1/3、`-use_sw_tback` |
| 別 alignment | CLI で無効化されている out-of-frame/frameshift 経路。通常の別 frame HSP と混同しない |
| 未対応追加調整 | `query_loc`、`subject_loc`、`dbsize`、`searchsp`、xdrop 三種の明示上書き、`qcov_hsp_perc`、best_hit_overhang/score_edge、`mt_mode` 明示指定 |
| 未対応表示 | 0/6/7 以外、未登録 custom field、`delim`、HTML、show_gis、num_descriptions/num_alignments/line_length/sorthits/sorthsps 明示変更 |
| 入出力拡張 | stdin query/subject、`parse_deflines` 明示指定、search strategy import/export |

未公開オプションの **既定動作** は対応範囲に必要な限り移植する。
例えば xdrop の明示上書きが範囲外でも、NCBI 既定 xdrop 計算は必須である。
BLASTP 特有 unified-P の `u` suffix を BLASTX に持ち込まない。
既存 CLI v2 の完全 token 検証と NCBI BLASTX の語彙を切り分けて説明する。

段階 A で固定 `blastx -help` の全キーを機械的に列挙し、X01–X18 の成功テスト、
この表の拒否テスト、help/version の既存製品規約のいずれかへ必ず割り当てる。
本書から漏れたキーを自動的に受理しない。将来、比較元三プログラムの適用可能な
実対応が広がった場合は scope ledger を更新し、同等性の下限を再評価する。

## 4. NCBI ソース対応表

以下の相対パスは固定 NCBI checkout の `c++/src/` 以下を指す。
行番号は調査版の入口・確認箇所であり、移植時には caller と callee の全範囲を
読んでソース対応記録に保存する。単一関数の一致をパイプライン一致の証明にしない。

| 処理 | NCBI owner / 確認箇所 | 移植時に固定する内容 |
| --- | --- | --- |
| CLI と task | `algo/blast/blastinput/blastx_args.cpp:44-149` | 引数群の登録・抽出順、task、query code、batch size |
| default | `algo/blast/api/blastx_options.cpp:47-90`、`blast_prot_options.cpp`、`blastinput/blast_args.cpp:257-420,825-899` | 閾値、matrix/gap、SEG、composition、sum stats、省略状態 |
| app 順序 | `app/blast/blastx_app.cpp:207-295` | FASTA→batch→CLocalBlast→formatter、empty input、query 順 |
| 入力バッチ | `algo/blast/blastinput/blast_input_aux.cpp:119-136`、`blast_input.cpp` | 通常 blastx の 10002、実際の batch 切替と例外条件 |
| query 分割 | `algo/blast/api/split_query_cxx.cpp`、`split_query_blk.cpp`、`core/split_query.c` | 一配列の chunk/overlap、context offset、再統合 |
| query context | `algo/blast/core/blast_query_info.c:59-101,141-220`、`api/blast_setup_cxx.cpp` | 六 context、無効 frame、元長、effective search space |
| 翻訳 | `algo/blast/core/blast_util.c`（長さ関数 `:923`）、`api/blast_aux.cpp:588-603`、`objects/seqfeat/gc.prt` | NCBISTDAA、NULLB、曖昧コドン、末端・負鎖、遺伝暗号 |
| マスク | `algo/blast/core/blast_setup.c:615-652`、`blast_filter.c`、`blast_seg.c` | translated AA mask、DNA mask、lookup と extension、nomask |
| lookup/seed | `algo/blast/core/blast_aalookup.c`、`blast_aascan.c`、`aa_ungapped.c` | word 3/5、近傍語、overflow、one/two-hit、対角管理 |
| preliminary search | `algo/blast/core/blast_engine.c`、`blast_gapalign.c`、`blast_parameters.c` | context 選択、cutoff、gapped HSP、containment、初期 pruning |
| 統計 | `algo/blast/core/blast_setup.c:699-847`、`blast_stat.c`、`blast_parameters.c` | query AA 長、protein subject 長、有効長、丸めと精度、cutoff 更新 |
| Kappa | `algo/blast/core/blast_kappa.c:2309,2424-2427,2981-3060,3244-3250`、`composition_adjustment/*` | translated query 情報、frame 数、redo 入力、行列、heap、early termination |
| linking | `algo/blast/core/link_hsps.c:1444-1462,1613-1810` | BLASTX の gap_q/gap_s、same/uneven gap、連結と E-value |
| traceback/pruning | `algo/blast/core/blast_traceback.c`、`blast_hits.c:1330,2268,2455`、`blast_itree.c` | edit script、同端点、包含、同点、スコア順・E-value 順の時点 |
| 座標 | `algo/blast/api/blast_seqalign.cpp:659-745,1324-1380`、`objtools/align_format/tabular.cpp` | gapped の `s_BlastHSP2SeqAlign`、ungapped の `x_UngappedHSPToStdSeg` から formatter へ。query nt / subject aa、両鎖 |
| filter | `algo/blast/api/setup_factory.cpp:296-400`、`core/hspfilter_culling.c`、`core/blast_hits.c:2537-2606`、`core/blast_engine.c:587-591`、`core/blast_traceback.c:859` | writer/pipe 登録と subject_besthit の preliminary/traceback caller。範囲外の best-hit pipe と混同しない |
| formatter | `algo/blast/format/blast_format.cpp`、`objtools/align_format/tabular.cpp`、`showalign.cpp`、`format_flags.cpp` | 0/6/7、IDs、frames、BTOP、no-hit、統計部、headers |
| unit tests | `algo/blast/unit_tests/api/{bl2seq,blastfilter,linkhsp,split_query}_unit_test.cpp` | 既存の BLASTX 境界条件と独立した期待値 |

表に挙げた入口から呼出し先を追い、段階 A で全分岐の owner を補完する。
対応 NCBI 実装が見つからない機能は新規推測実装せず、その要求の根拠を再検討する。
Rust を変更する際は各変更の直上に **正確なパス・行番号・C/C++ 断片**を付ける。
`Blast_HSPGetAdjustedOffsets` という名前だけで translated report の変換を決めず、
BLASTX CLI が実際に通る Seq-align 生成・formatter の caller を権威にする。

## 5. Rust 構成と再利用方針

### 5.1 追加・変更 owner

| 責務 | 予定 owner | 方針 |
| --- | --- | --- |
| BLASTX パイプライン | 新規 `LOSAT/src/algorithm/blastx/` | args、query setup、search、statistics/redo、report の責務を分ける。細分化は NCBI owner に合わせる |
| 公開 CLI | `src/cli.rs`、`src/main.rs`、`src/algorithm/mod.rs` | `blastx` dispatch と対応表に基づく validation |
| 翻訳・遺伝暗号 | `algorithm/tblastx/translation.rs`、`utils/genetic_code.rs`、既存 core | 重複表を追加しない。BLASTX context 構築は固有 owner に置く |
| AA lookup/extension | `algorithm/tblastx/lookup/`、`algorithm/blastp/{extension,gapalign}.rs` | プリミティブと scratch を再利用。BLASTP 検索全体を frame ごとに呼ばない |
| 組成補正・統計 | `core/composition_adjustment/`、`stats/`、BLASTP/TBLASTN の対応箇所 | NCBI で共有される計算だけ共有。translated query 分岐を明示する |
| HSP/filter/report | `post/`、`report/`、必要最小限の共通型 | 内部座標を保持し、format 前に一度だけ変換 |
| library / Web | `api/local_blast.rs`、`web_api.rs`、既存 wasm host | 既存所有権・ABI に BLASTX dispatch を追加。別の検索実装を作らない |
| 検証 | `LOSAT/tests/`、`docs/evidence/losatx_stage_*/` | source/run/fixture identity と raw diff を記録 |

上記は所有責任であり、同名ファイルを無条件にすべて作る指示ではない。
既存 TBLASTN の private helper を公開するために巨大な一般化を行わない。
必要最小限の共通化には旧 caller と新 caller の両方の NCBI 根拠を記す。

### 5.2 処理の順序と単位

```text
CLI / Web 入力検証
  → NCBI-compatible query batch
  → 六 context の翻訳・mask・score block・lookup
  → protein subject ごとの seed / ungapped / preliminary gapped HSP
  → preliminary pruning / statistics / linking / collection
  → NCBI が選ぶ composition redo または通常 traceback
  → 再評価・統計・包含除去・hitlist・post-filter
  → query 順 / subject 順 / HSP 順の確定
  → query nucleotide 座標と subject amino-acid 座標へ変換
  → 0 / 6 / 7 の formatter
```

この図は工程の目次である。詳細な link/reap/redo/filter の反復と呼出し順は
段階 A の call graph、段階 C/D の trace で確定し、図の単純化を実装に持ち込まない。
sum_stats=false、composition=0、ungapped で通る分岐も独立して記録する。

### 5.3 保存すべき状態

- 元 query の長さ、query index、context/frame、連結 offset、chunk の元座標。
- sentinel を含む検索用配列と未マスク配列。identity を検索用 X 配列から数えない。
- subject の入力 OID と表示 ID。ID が同じでも入力 record を勝手に統合しない。
- HSP の score、context、canonical offsets、gapped starts、edit script、link 情報。
- context ごとの Karlin/length adjustment/cutoff、各 caller が使った search space。
- redo 前後の matrix、scale factor、rounding、heap insertion/replacement と削除理由。

全コードで同じ浮動小数点型・計算順・NCBI 比較器を保つ。
TBLASTX の翻訳 subject 用 length adjustment や TBLASTN の DB 長÷3 を
BLASTX の protein subject に適用しない。

## 6. テストデータと取得・固定方法

### 6.1 データ層

| 層 / ID | 入力 | 検証目的 | 実行区分 |
| --- | --- | --- | --- |
| S: 合成小型 | 下記 S01–S15 | 最初の相違と境界条件を小さい入力で特定 | PR gate |
| N: NCBI 由来 | BLASTX、mask、link、split-query の公式 unit fixture | NCBI が既に固定している分岐・不具合防止 | PR / focused |
| R01 | `AvCLPV.fasta → AvCLPV.faa` | 標準コード self、416069 nt / 120 proteins | 短縮版 PR、全長 release |
| R02 | `AvCLPV.fasta → PsCLPV.faa` と逆方向の核酸/蛋白組 | 近縁 cross、複数 HSP、composition、PsCLPV 179 proteins | release |
| R03 | `AP027131.fasta → AP027131.faa`、query code 4 | 非標準コード self、662108 nt / 595 proteins | 短縮版 PR、全長 release |
| R04 | `AP027131.fasta → AP027133.faa` と逆方向、query code 4 | 非標準コード cross、AP027133 543 proteins | release |
| R05 | `MelaMJNV.fasta → PemoMJNVA.faa` | 別の配列群、287061 nt / 115 proteins | release |
| R06 | 上記核酸の固定 fragment 集合→対応・近縁 FAA | 150/300/1000 nt 等の短い query、部分 CDS、逆鎖、no-hit、batch | PR / release / benchmark |
| R07 | `EDL933.fna` / `Sakai.fna` と対応 GB から抽出した FAA | bacterial code 11、長 query、多 subject、メモリと並列 | extended / release |
| R08 | HBB `NM_000518.5 → NP_000509.1` と `NG_000007.3` の HBB 注釈領域→同蛋白 | transcript と intron を含む genomic query、code 1、複数 HSP/linking | 段階 A で取得・領域・bytes を固定、release |

既存実データの場所は [LOSAT/tests/fasta](../LOSAT/tests/fasta)。
R01/R02 の核酸 accession.version はそれぞれ `LC738883.1`、`AP027154.1`。
R03/R04 は `AP027131.1`、`AP027133.1`。ファイル名を checksum の代わりにしない。
表のサイズは計画作成時の集計であり、実行時に record 数・長さ・checksum を再記録する。

R07 の FAA 抽出は GenBank `/translation` と protein_id を優先し、location、
codon_start、transl_table、pseudo/partial の扱いを manifest に残す。
LOSAT の新規翻訳コードで生成した蛋白だけを正解側に置かない。
R08 の transcript/protein 対は
[NCBI の NP_000509.1 レコード](https://www.ncbi.nlm.nih.gov/protein/NP_000509.1)で確認した。
genomic 入力候補は [NG_000007.3](https://www.ncbi.nlm.nih.gov/nuccore/NG_000007.3)。
今回配列そのものは取得していないため、段階 A で GenBank の HBB annotation、
切出しの元座標、向き、元 record 全体と派生 FASTA の checksum を検証・固定する。
transcript と genomic query の出力が同じになるとは仮定しない。
ネット接続がない CI は保存済み fixture で完結させる。
大きい NR/Swiss-Prot 全体は初期 gate に不要。
R08 の生物例だけで intron linking の分岐を保証せず、S08/N fixture と併用する。

### 6.2 合成 fixture の必須一覧

| ID | 生成・選定する条件 | 必須観測 |
| --- | --- | --- |
| S01 | 同じ非低複雑性蛋白を +1/+2/+3/-1/-2/-3 の各核酸 frame に配置 | frame、座標、qseq、score、query/subject 長 |
| S02 | 長さ 0/1/2/3 と word が作れる前後、3n/3n+1/3n+2、末端一致 | 無効 context、sentinel、部分コドン、警告/no-hit |
| S03 | 全 64 codon、全 26 code、IUPAC ambiguity を含む識別 probe | 翻訳表、曖昧コドン、stop、コード別 end-to-end |
| S04 | 内部/末端 stop、連続 stop、全 N、曖昧蛋白 B/Z/J/X/U/O/* | encoding、extension 停止・許容、再評価、identity |
| S05 | SEG 境界、小文字 mask、plus/minus、soft/hard、全 mask | lookup/extension/nomask の区別、mask の座標変換 |
| S06 | one-hit/two-hit、window 境界、threshold 前後、word 3/5 | seed 集合、対角更新、圧縮 lookup overflow |
| S07 | insertion/deletion、X-drop 境界、同端点と包含、同 score | preliminary/traceback/edit script、survivor と tie order |
| S08 | 同 frame・別 frame の複数 HSP、間隔/重なり/intron 境界 | link 可否、sum_stats、連結集合、query 側 gap 条件 |
| S09 | biased composition、短長の差、adjust/no-adjust、複数 subject | Kappa matrix rule、redo window、scaled score、early termination |
| S10 | 同点 subject 8 個以上、入れ替え、max targets/HSP 直前直後 | 入力 OID、heap、削除順、並列集約 |
| S11 | hit/no-hit/invalid-context query の混在、10002 前後の総長 | batch 境界、出力 footer、異常 query の影響範囲 |
| S12 | 一配列の NCBI split 閾値の前後、overlap を跨ぐ HSP | split/rejoin、重複除去、負鎖 frame、元座標 |
| S13 | local/accession/pipe ID、説明文、重複 ID、複数 record | qseqid/qacc/qaccver、stitle、順序、ヘッダー |
| S14 | target/HSP/culling/subject_besthit と E-value の交差 | 全 filter の適用順・範囲・閾値 |
| S15 | 不正 CLI、I/O failure、pool failure、再利用後の失敗→成功 | 明示エラー、結果 clearing、worker 回収、状態漏れなし |

合成 generator は seed、生成規則、元配列、座標、reverse-complement 手順を固定する。
「stop を越えないはず」「この HSP は一つのはず」のような手作りの検索期待値は禁止。
期待 report は固定 NCBI CLI、内部段階は固定 NCBI source oracle から生成する。
翻訳の unit expected は NCBI のテーブル/翻訳実装に基づき、テスト対象 Rust 関数と
同じコードを使って再計算しない。

### 6.3 NCBI unit fixture の移植

- `bl2seq_unit_test.cpp:1638` の両鎖、`:1876` の plus、`:1890` の minus。
- `blastfilter_unit_test.cpp:949` 以降の BLASTX lowercase mask。
- `linkhsp_unit_test.cpp:542` の uneven-gap BLASTX と cutoff/境界ケース。
- `split_query_unit_test.cpp` の translated-query 分割と context のケース。

古い GI や object-manager のオンライン取得へ CI を依存させない。
必要な配列を NCBI 同梱 ASN.1 または明示した取得元から FASTA に materialize し、
元の testcase・入力 identity・変換手順・利用条件を記録する。
NCBI の API 期待値と CLI の formatter/batch default の違いを明示し、API 出力を
未校正のまま CLI の期待 report に使わない。

## 7. 比較 matrix と合格条件

### 7.1 必須 matrix

| Matrix | 必須構成 | 合格条件 |
| --- | --- | --- |
| M1: 遺伝暗号 | 26 codes × 3 formats（0/6/7）× native threads {1,4} = **156 unique cases** | 各行で同じ入力・query code の NCBI CLI 生バイト一致 |
| M2: strand/frame | 6 frame × strand 指定、長さ剰余 0/1/2、コード 1/4/11 | 各出力と無効 context が NCBI に一致 |
| M3: 検索 profile | word 3/5 × composition 0/2 × SEG no/yes/explicit × window 0/default/boundary | 有効組合せの分岐到達証拠と raw parity |
| M4: 統計 | gapped、ungapped+mode0、sum_stats true/false、intron/長さ/E-value 境界、S08/N/R08 | HSP 集合、連結、score、full-precision 診断、report が一致 |
| M5: 制限/filter | max targets/HSP の省略/1/境界、culling、subject_besthit、同点/入力順 | 除去結果と report 順が一致 |
| M6: formatter | 0/6/7、全 30 custom fields、std、列順/重複、IDs、no-hit/multiquery | 本文・見出し・脚注・空白・丸めを含む一致 |
| M7: batch/long | S11/S12、R01–R07 の全長または固定 multiquery | 分割による欠落・重複・frame shift・統計変化なし |
| M8: Native parallel | threads 1/2/4/8、worker 未満/同数/超過 subject、同点、再実行 | 全出力同一、実並列の稼働証拠 |
| M9: WASI/Web | 第 10 節の各入口・各 mode と同一入力/解決済み options | Native と同じ report、失敗時 lifecycle 合格 |
| M10: negative | 全範囲外キー・値・不正組合せ・入力と I/O/pool failure | 所定のエラー契約、黙示 fallback なし |
| M11: existing regression | BLASTN/megablast、BLASTP、TBLASTX、TBLASTN | 既存の各 authority/exception を維持 |

M1 の入力には各 code の結果に意味を持つ probe と各 strand の hit を含める。
全 64 codon の unit matrix と M1 の検索 matrix は別であり、片方で代用しない。
156 は M1 の最低行数であり、全体の認証件数ではない。
同じ行の再実行を unique coverage として重複計上しない。

既知に相互作用する `code×frame×mask`、`composition×linking×length`、
`ties×limits×threads`、`batch×frame×report` は専用ケースを必須にする。
任意の全直積を無計画に増やすのではなく、残る独立設定は pairwise coverage と
境界値で補完する。選択規則を generator/manifest に残す。

### 7.2 比較の手順

1. program/task、実行 target、版/hash、query/subject bytes、全 argv、環境を固定。
2. 現在の release build と固定 NCBI を新しい作業ディレクトリで実行。
3. stdout または `-out` の report を加工せず byte compare。通常 stderr と exit code
   も保存し、成功・警告・失敗の契約を照合する。
4. 不一致なら最初の byte、query/subject、HSP、field を特定する。
5. seed→preliminary→link/reap→redo→post-filter→formatter の初回相違を trace。
6. 同じ NCBI owner の全 offender を分類・修正してから focused gate を再実行。
7. 関連 program/target gate、最後に必要な回帰・認証を実行する。

行の sort、E-value tolerance、空白・ヘッダー除去、同点 HSP の独自 canonicalization
で合格にしない。診断用の normalized table は raw output と別保存する。
raw score/full-precision E-value など通常 outfmt 6 にない情報は、追加 custom 出力
または比較専用 trace で調べる。trace あり/なしの通常 stdout/stderr も校正する。
作業 directory、入力の lexical path、locale、改行、浮動小数点・SIMD build 条件も
固定する。`BATCH_SIZE`、`OLD_FSC` 等の NCBI 実験用環境変数と LOSAT 診断変数の
継承を防ぎ、有効な environment を記録する。バッチ実装の検証を環境変数による
小型化だけで代用せず、通常環境の境界ケースを必ず残す。

### 7.3 エラー契約の区別

NCBI 自身の有効検索で生じる no-hit、全無効 context、配列警告は、NCBI の
exit code と report/警告のタイミングを再現する。
CLI v2 の未知引数、未対応機能、help/root version と Web の status/error buffer は
既存 LOSAT のインターフェース規約に従う。NCBI の usage 全文コピーとは分ける。
その違いを機械可読 ledger に明記し、検索出力の差をエラー表現差で隠さない。

## 8. 段階別実装計画

各段階は prerequisites、変更 owner、成果物、exit gate を持つ。
前段の未解決必須差分を「後で直す」として認証段階へ進めない。
内部 unit 用入口と未完成の公開検索を区別し、公開検索は実装済み能力だけを受理する。

### A. 権威・scope・fixture の固定

作業:

- 作業ツリーと基点を記録し、既存変更を保護する。
- 固定 NCBI の source commit、binary hash、build provenance、help、実行環境を保存。
- 第 3 節のキー、値域、default、派生 default、相互排他、拒否を ledger 化する。
- BLASTX call graph と source→Rust の owner map を作る。SEG/default と batch を校正。
- S01–S15 と R/N 層の最小セットを固定し、NCBI raw expected を生成する。
- 既存 program の適用可能な gate を登録し、共有部変更前の baseline を固定する。

成果物（予定）: `scope.tsv`、`source_map.tsv`、`oracle.json`、`fixtures.tsv`、
`defaults.json`、`coverage.tsv`、生成スクリプト、NCBI expected、baseline 記録。

完了条件: 必須要求に owner と test ID がある。NCBI help の未分類キーがない。
入力と期待値を再現でき、API/CLI default の未説明差を残さない。

### B. CLI・query preparation・翻訳

作業:

- `BlastxArgs`、解決済み options、CLI dispatch と明示拒否を実装。
- FASTA query/subject、batch、六 context、genetic code、strand、mask を接続。
- 共有 sentinel/translation と protein encoding の単位・offset を合わせる。
- 全 26 code×64 codon、ambiguity、S01–S05、S11 の内部状態を照合。
- public engine は C–E が完成するまで検索未実装を明示する。

完了条件: option 解決、配列 bytes、context、mask が一致。
必須能力の単なる parser 受理を検索対応済みと表示しない。

### C. seed・preliminary search・gapped/ungapped

作業:

- word 3/5、one/two-hit、context-aware 対角管理と cutoff 計算を接続。
- protein subject の順序、preliminary gapped extension、common endpoints と
  containment、HSP 保存タイミングを移植。
- ungapped+mode0 の分岐を NCBI の独立経路として実装。
- 長 query の chunk/overlap を移植し、再統合前後を照合。
- S06/S07/S10/S12、実データ縮小版で候補・raw HSP を照合する。

完了条件: 選んだ全 C fixture の seed/preliminary HSP、score、座標、順序が一致。
NCBI の記録済み cutoff を Rust の実際の計算の代わりに渡す診断は、C 完了に数えない。

### D. 統計・linking・Kappa・filter

作業:

- translated query と protein subject の effective length/search space を実装。
- sum statistics、BLASTX の query 側 gap/intron、link/reap の呼出し順を移植。
- six-context query info、composition mode 2、redo、行列・scale、score 正規化を接続。
- mode 0 の通常 traceback、identity/positive、post-traceback pruning を移植。
- max targets/HSP、culling、subject_besthit、heap、early termination を接続。
- S08/S09/S14 と全コード・実データ縮小版で HSP membership と統計を照合。

完了条件: 最終 HSP 集合、順序、score、full-precision 診断が一致。
内部で 6 回の独立 BLASTP を実行して出力を結合する代替経路を残さない。
TBLASTN の translated subject 分岐や BLASTP 固有の統計選択を流用しない。

### E. report・public Native CLI

作業:

- query nt / subject aa 座標、frame、ID、alignment、BTOP、全 30 field を実装。
- outfmt 0/6/7 の prolog、query ごとの no-hit、statistics、footer を実装。
- batch の reset/format 時点、stdout/`-out`、通常 stderr/exit を照合。
- M1–M7/M10 を Native serial で通し、完成した公開 CLI へ接続する。
- 省略既定値・明示同値・custom format の全能力を help と文書へ反映する。

完了条件: 必須 Native serial matrix に raw byte 差がない。
内部だけ一致、outfmt 6 だけ一致、code 1 だけ一致では E 完了にしない。

### F. Native parallel・WASI・Web/reactor

作業:

- 独立 subject 計算を既存 search-scoped pool に配分し、入力 OID 順に回収。
- collector/heap/redo の order-sensitive 部分は NCBI 順で処理し、並列化するなら
  その同値性を別証拠で示す。
- command WASI serial/threaded と Web/reactor の両入口へ同じ engine を接続。
- thread 数変化、再利用、pool failure、invalid handle、result clearing を検証。
- M8/M9 と同点・多数 query/subject、codes、formatter の横断ケースを実行。

完了条件: serial/threaded の bytes 一致、実並列稼働、lifecycle 合格。
thread target の build 成功だけで並列対応を認証しない。

### G. 総合 gate・性能・独立監査・release handoff

作業:

- 全必須 matrix の coverage と結果を集計し、欠落・skip・重複計上を拒否。
- 基点にある BLASTN/megablast、BLASTP、TBLASTX、TBLASTN と、
  build/test/fmt/clippy、Native platform の必要 gate を同じ候補で実施。
- 第 11 節の性能測定を実施し、全 raw samples と output hash を保存。
- `ncbi_parity_auditor` に source/output/order/benchmark の独立 read-only 監査を依頼。
- README、CHANGELOG、release scope、manifest、使用例を実証した能力に揃える。
- 正確な candidate SHA と build/config/input identity で認証記録を固定する。

完了条件: X01–X18 に未解決項目がなく、第 14 節の checklist がすべて満たされる。
タグ、release 公開、package publication は別途の許可境界を維持する。
この計画書とセッションプロンプトの commit/push は今回の依頼範囲に含まれる。

## 9. Native platform と既存契約の維持

既存配布の Linux x86_64、Windows x86_64、macOS arm64/x86_64 を認証対象に含める。
追加 architecture は段階 A で現在の release workflow と対応表を照合する。
各 platform で source identity、toolchain、build flags、binary hash を残す。

[PD-NCBI-PLATFORM-VARIANCE](product_decisions/PD-NCBI-PLATFORM-VARIANCE.md)の
既存 PR 6 契約を変更しない。

- 既存 Gate A: LOSAT は凍結 PR 5 raw bytes に一致する。
- 既存 Gate B: 六つの official NCBI search は登録済み platform fingerprint に一致する。
- platform-local NCBI の bytes を既存 LOSAT expected に置き換えない。
- 既存 BLASTN 同点分類、TBLASTX/TBLASTN の genetic-code 例外を BLASTX へ拡張しない。

BLASTX は既存六検索の登録対象ではない。BLASTX 用の source authority と基準
raw output を段階 A/E で固定し、各配布 platform の oracle を特性化する。
不一致が出た場合は新しい authority version と review が必要であり、既存の
登録表にない fingerprint を黙認しない。platform 別 LOSAT アルゴリズムで吸収しない。
既存 PD は BLASTX への拡張を承認していない。BLASTX 用の bounded platform
authority を新設・拡張して受け入れるには、特性化と review に加え、明示的な
製品・governance 承認を必要とする。それまでは BLASTX の未説明 raw variance は
hard-fail のままとし、既存 PD を一般的な許容差として使わない。
対象 platform の検証が未実施なら、その対象を認証済みとした v0.2.0 handoff は不可。

## 10. Wasm と Web/reactor の検証

| 入口 | 必須確認 |
| --- | --- |
| Native CLI | threads 1/2/4/8、同一 query/subject/options、stdout/`-out` |
| `wasm32-wasip1` command | 実際は serial。Native と同一 bytes。スレッド数指定の扱いは既存 serial-target 契約に固定 |
| `wasm32-wasip1-threads` command + `wasm-threads` | threads 1/2/4/8、共有 memory、実 worker spawn/join、同一 bytes |
| serial reactor | 文字列入力・保存済み FASTA handle、0/6/7/custom、連続呼出し |
| threaded reactor | 上記に加え pool サイズ変更、spawn failure、回収、回復後の成功 |

`web_api.rs` の `losat_web_run_pair` と `losat_web_run_pair_handles` に `blastx` を追加し、
同じ Rust engine を使う。既存 ABI、alloc/dealloc、store/release/clear、result/error
buffer の所有権を維持する。新しい公開 SDK の安定性保証や gbdraw 側の UI 改修は
この計画の前提にしない。

Web の既定 outfmt 6 と CLI の既定 0 は既存入口規約として明示し、NCBI との比較は
**同じ解決済み options** で行う。Web の ID/title 処理も manifest に記録し、
必要な共通 FASTA identity を用意して raw report を比較する。
エラー時に前回の結果を返さず、失敗後の再実行を成功させる。
異なる gencode/mask/options を同じ handle で交互に実行して cache 汚染を調べる。

既存 `check_wasi_reactor.js`、`check_wasm_threading.py` 等へ BLASTX を追加する。
同一入力の反復で worker と memory が際限なく増えないことを確認する。
Node の reactor 検証とブラウザ実行を区別し、ブラウザ向け配布入口では実ブラウザで
必要な worker/shared-memory 設定と成功・失敗の smoke を行う。
offline bundle/依存関係を変更する場合に限り browser-offline-qa の専用 gate を適用する。

## 11. 性能評価

パリティ確立後に性能を測る。候補集合・pruning・統計を変える速度改善は禁止する。
初回 LOSATX に旧認証版がない場合、既存別 program を baseline と偽らない。
最初の正しい serial BLASTX を baseline とし、その後の並列化・最適化を比較する。

| workload | 対象 |
| --- | --- |
| B1 小型 | 単一短 query / 少数 proteins。startup を含む利用時 latency |
| B2 fragment batch | R06 の固定多数 query / protein set。batch と throughput |
| B3 全長 viral | R01/R02/R05。長 query、composition、複数 HSP |
| B4 code 4 | R03/R04。非標準 query translation と linking |
| B5 bacterial | R07。多数 subjects、並列 scaling、peak memory |

各 case/mode で非計測 warmup **1 回**、保持する計測 **3 回ちょうど**。
中央値と全三標本の min–max、wall time、取得可能な CPU time/peak RSS、
出力 hash、threads、入力/実行ファイル/toolchain identity を報告する。
追加標本はユーザーの明示要求または三標本が実際に判定不能な場合だけ行い、理由を記録。
最速値の選択、遅い標本の無断除外はしない。

NCBI の threaded search 時間は `makeblastdb -dbtype prot` を事前実行した `-db` で測る。
DB 作成の command/version/input checksum/経過時間は別記録として残す。
LOSAT local-subject と NCBI DB の統計・I/O 契約の違いを明示し、同条件の speedup と
短絡しない。BLASTX の parity/hit-distribution oracle は同じ query_gencode を
指定した NCBI **local `-subject`** とし、timing DB output とは別名・別 metadata にする。

Wasm は cold command と warm reactor を分ける。ロード、コンパイル、worker startup
をどこまで含むかを固定する。スレッドを増やして遅い workload も結果に含める。
一定倍率の NCBI 超えを根拠なく release 条件にしない。
共有部変更の既存 program 性能悪化は同条件 baseline と三標本で評価し、明確な悪化は
原因を解消するか、具体的な tradeoff の製品判断を得るまで未解決にする。

長時間 benchmark の status poll は原則 10 分間隔。
ユーザーへの作業状況更新と、実行中プロセスへの頻繁な status poll を区別する。

## 12. 検証・証拠ファイルの構成

以下は今後作る予定の配置であり、本書の保存時点で存在・合格を主張しない。

```text
LOSAT/tests/
  blastx_parity_manifest.tsv
  compare_blastx_parity.py
  blastx_cli.rs
  fixtures/blastx/
    synthetic/ ncbi/ real_subsets/
    generate_fixtures.py
    provenance.json
docs/evidence/losatx_stage_a/ ... losatx_stage_g/
  scope.tsv
  source_map.tsv
  coverage.tsv
  environment.json
  comparison.jsonl
  first_differences/
  expected/ actual/ traces/
  evidence.sha256
  STAGE_<letter>_GATE.md
```

manifest の最低項目:

```text
case_id, requirement_ids, fixture_family, input_paths, input_sha256,
accession_version_or_generator, program, task, full_argv,
resolved_options, genetic_code, strand, format_spec,
target, features, threads, entrypoint, repeat_id,
losat_source_identity, losat_binary_sha256,
ncbi_source_identity, ncbi_binary_sha256, oracle_kind,
expected_exit, expected_stdout_sha256, expected_stderr_sha256,
actual_exit, actual_stdout_sha256, actual_stderr_sha256,
first_difference, result, authority_version
```

scope 行の状態は `REQUIRED` / `UNSUPPORTED` / `NOT_APPLICABLE`、実装・検証状態は
別列にする。`SKIP`、未実行、既知差分を `PASS` に集計しない。
機械可読 coverage は requirement→source owner→production owner→test case→evidence
を双方向に辿れるようにする。必須 code/format/thread の欠落や重複も gate failure。

既存 evidence は code/input/environment/acceptance が同じ場合だけ再利用する。
新しい BLASTX 結果に TLOSAN/BLASTP の過去の合格件数を加算しない。
巨大な再生成可能 output を無条件に Git に入れず、再現に必要な小型 fixture と
恒久的な provenance・hash・取得手順を保存する。外部一時ディレクトリだけを指す
失われるリンクで release 証拠を完結させない。

## 13. 実行コマンドと確認の順序

以下の BLASTX コマンドは実装後の例であり、現在の CLI で実行済みではない。

```bash
# 同一 local-subject contract、code 4、標準 tabular
LOSAT blastx -query LOSAT/tests/fasta/AP027131.fasta \
  -subject LOSAT/tests/fasta/AP027133.faa \
  -task blastx -query_gencode 4 -outfmt 6 -num_threads 1 \
  -out /tmp/losatx-candidate.out

blastx -query LOSAT/tests/fasta/AP027131.fasta \
  -subject LOSAT/tests/fasta/AP027133.faa \
  -task blastx -query_gencode 4 -outfmt 6 -num_threads 1 \
  -out /tmp/losatx-ncbi.out

cmp /tmp/losatx-ncbi.out /tmp/losatx-candidate.out
```

runner は上記の固定例名をそのまま共有せず、実行ごとに新規 directory を作り、
NCBI の絶対 path/hash を検証して stdout/stderr/exit を保存する。

関連 parity sweep 完了後の基本 gate:

```bash
cargo build --manifest-path LOSAT/Cargo.toml --release --locked
cargo test --manifest-path LOSAT/Cargo.toml --locked
cargo fmt --manifest-path LOSAT/Cargo.toml -- --check
cargo clippy --manifest-path LOSAT/Cargo.toml --all-targets --locked -- -D warnings
cargo build --manifest-path LOSAT/Cargo.toml --release --locked --bin LOSAT \
  --target wasm32-wasip1 --no-default-features
cargo build --manifest-path LOSAT/Cargo.toml --release --locked --bin LOSAT \
  --target wasm32-wasip1-threads --features wasm-threads
```

reactor は既存 workflow の CRT 初期化・link 設定で別 artifact として build する。
command wasm を reactor として代用しない。
比較 runner と既存 regression script は cleanup/input/oracle を読んでから実行する。
既知の NCBI 分岐欠落を残したまま「役に立つかもしれない」広域テストを繰り返さない。
production、tests、docs、generated evidence の diff を別々に確認し、変更点が増えた
部分だけ再確認する。

## 14. v0.2.0 完了 checklist

- [ ] X01–X18 の必須能力が Rust で実装され、未対応で代替していない。
- [ ] NCBI source map、default、全 help key の成功/拒否分類が固定されている。
- [ ] 全 26 query genetic codes、six frames、mask、末端・曖昧コドンが一致する。
- [ ] word 3/5、gapped/ungapped、composition 0/2、sum stats/link/filter が一致する。
- [ ] 0/6/7 と全 30 custom fields、no-hit、multiquery、batch/long query が raw parity。
- [ ] 必須 matrix に missing、skip、未知 fingerprint、未説明の byte 差がない。
- [ ] Native 1/2/4/8、serial/threaded WASI、Web/reactor の指定入口で同じ結果になる。
- [ ] 実並列、worker 回収、反復 memory、エラー後回復を確認している。
- [ ] 配布 Native platform の gate と既存 Gate A/B の契約を維持している。
- [ ] BLASTN/megablast、BLASTP、TBLASTX、TBLASTN に退行がない。
- [ ] build/test/fmt/clippy と適用可能な package/artifact smoke が合格している。
- [ ] NCBI が runtime/build/fallback 依存になっていない。
- [ ] 性能は warmup 1+計測3、全標本・中央値・範囲と output hash を記録している。
- [ ] 独立 read-only audit の指摘を解消し、適用範囲を明記している。
- [ ] 文書・help・manifest・release 記録が同じ source candidate と能力を示す。

## 15. 主なリスクと解消条件

| リスク | 検出・解消 |
| --- | --- |
| 六つの BLASTP として扱って統計が変わる | context/Kappa/link の NCBI trace。BLASTX 固有 pipeline を維持 |
| API default と CLI default の混同 | 無指定/明示設定の固定実行、option 抽出後の状態照合 |
| TBLASTN の query/subject 単位を誤転用 | nucleotide/AA/context 単位を型・記録で明示し両鎖 fixture を通す |
| custom output が検索精度から漏れる | 全 30 field と 0/7 本文を独立 gate にする |
| バッチ・chunk で frame/統計/順序が変わる | S11/S12 と全長実データを mandatory にする |
| 同点・thread scheduling が順序を変える | input OID、同点 HSP、limits、反復を横断検証 |
| NCBI platform variance を一般的許容差にする | authority version と未知 fingerprint hard-fail |
| 対応引数だけ増えて engine が未実装 | scope/implementation/evidence の状態を分離、public integration gate |
| 既存 TLOSAN 等を共通化で壊す | 変更した owner の全 caller と必要回帰を確認 |
| fixture が目的の分岐を通っていない | branch/stage trace で到達確認、no-hit の大量一致を過大評価しない |
| 時間不足で必須を任意化する | 未完了を明記し gate failure とする。下限変更は別の製品判断 |

## 16. セッション分割と引継ぎ

実装は原則 A→B→C→D→E→F→G。C/D が大きい場合は同じ exit 条件を維持して
小さな checkpoint に分ける。checkpoint は次段階の gate pass ではない。
セッションごとに `feature/losatx-blastx-v0.2.0` を使用し、commit/source identity、
変更 owner、現在の初回相違、通過 gate、未実行 gate、次の具体的作業、証拠 path を残す。

特に D の内部一致と E の formatter 一致、F の build 成功と runtime 一致を混同しない。
独立監査は release-facing parity/performance の受け入れ前に行う。
今回の計画書作成時点では、検索比較・benchmark・新規実装は行っていない。
