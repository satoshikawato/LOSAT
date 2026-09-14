# megablast 出力差と Wasm BLASTN・BLASTP 高速化の修正計画

作成日: 2026-09-13  
状態: ユーザー指示で現状を暫定採用・追加計測停止（2026-09-14）。全ゲート合格とは扱わない。証拠は [run-20260913-01](evidence/megablast_wasm_np_remediation_20260913/run-20260913-01/STATUS.md) に保存。  
計画時の HEAD: `521b0659f2b92ed8a71b6911ef685bfec539abaa`。多数の既存未コミット変更があるため、実行時の基準は HEAD だけで識別しない。

## 1. 目的と今回の到達点

次の三つを独立した成果目標として追跡する。

| ID | 成果目標 | 完了の条件 |
|---|---|---|
| G-M | 既存 megablast 出力差を解消する | NZ 自己比較と Sakai/MG1655 の現在の差を再現し、NCBI に根拠のある修正で raw bytes 一致を確認する。原因分類や既存差の再確認だけでは完了しない |
| G-N | Wasm BLASTN の検索時間を短縮する | 正しさを確認した基準版に対し、指定した二つの主入力の threaded n8 で各 5%以上短縮し、出力・他条件・メモリの採用基準を満たす |
| G-P | Wasm BLASTP の検索時間を短縮する | 正しさを確認した基準版に対し、指定した二つの主入力の threaded n8 で各 5%以上短縮し、出力・他条件・メモリの採用基準を満たす |

5% は実装開始後のユーザー指示（「5%以上でいいです」）による工学的な採用基準。前計画の 10% から変更する。NCBI の仕様ではない。二つの主入力の平均や最良入力だけで達成判定しない。serial n1 と threaded n1 は必須の診断・非回帰対象とし、単一スレッドの高速化を主張する場合は、その条件でも各 5%以上の短縮を別に証明する。

後続のユーザー指示「3.84%でもいいじゃん、ほかのやつはドカンと速くなったんだから。もともと短いわけでしょ」により、C-N4固定HSP方式について、短いBLASTN主入力 LC738874/LC738870 の5%下限を免除する。初回5反復の3.84%短縮を受け入れ、追加測定の変動も開示する。大きいBLASTN主入力とBLASTP両主入力の5%基準、全入力の出力・時間非回帰・メモリ・再利用検証は維持する。この限定変更を他入力の基準緩和に転用しない。

追加指示により、標準Wasm配布・通常検証はthreadedへ統一し、単一スレッド比較にも同じartifactのn1を使う。serialは明示指定の互換ビルド・検証とする。現在のgbdrawはAUTOやfallbackでserialを使用しているため、完全廃止はconsumer移行と別途扱う。既存serial測定、凍結された認証・release契約は保存する。この後の通常候補測定はnative n1/n8とthreaded n1/n8を主軸とし、serial性能測定は互換オプションに移す。

2026-09-14の追加承認により、基準版にもある反復時のWasm linear memory増加は共通runtimeの別件として残す。今回の採用では絶対的なplateau条件を外し、既定のRSS非回帰・メモリ予算、実保持量の比較、scratch再利用の検証を維持する。元の失敗結果は保存し、すべての観測memory byte数が同一とは主張しない。承認と適用境界は [runtime-memory-scope-decision.json](evidence/megablast_wasm_np_remediation_20260913/run-20260913-01/runtime-memory-scope-decision.json) に記録する。時間・出力・速度下限の条件は変更しない。

さらにユーザーからTLOSATXにも同じデータ移動削減を依頼された。BLASTPの判定を優先し、その後TBLASTXのHSP sort/link/replayを実測して、NCBIと同じ本体保持・参照操作へ改善する。TBLASTXの新しい短縮率下限は指定されていない。現在の出力・順序・遺伝暗号の契約と性能非回帰を維持し、既存の解決済みparity問題は回帰対象として扱う。

主実行環境は現在の Node/WASI command とする。ブラウザや別 runtime の高速化、native と同速、n8 の常時優位は保証しない。compiled module 再利用と同一 reactor 反復は別の測定軸として検証する。

計画作成時の依頼は計画書の保存まで。現在は後続の実装依頼に基づき進める。未実施の実装・テスト・リリース作業を実行済みとは扱わない。

## 2. 出発点と変更する進め方

以下は保存済みの観測であり、新しい基準測定には流用しない。

| 対象 | 既存記録 | 今回の判断 |
|---|---|---|
| BLASTN batch 16→32 | 一方で 8.03%短縮、もう一方で 18.72%悪化 | 棄却を維持する。別の値を根拠なく総当たりせず、待ち時間・計算量・負荷の偏りの未測定部分を先に調べる |
| BLASTP DP 行書き込み | 4.39%悪化 / 2.46%短縮 | 棄却を維持する。重い関数の一部分を簡略化したことと、全体時間を減らせることを区別する |
| BLASTP の追加三候補 | SEG cache 再利用、DP 関数境界、行列特殊化はいずれも不採用。最大 5.45%短縮 | 新たな原因証拠なしに復活・積み重ねしない |
| NZ 自己比較 megablast | 454 行中 1 行の identity/length/mismatch/gapopen 差。最初の因果的な不一致は未特定 | 最優先で段階ごとの状態を追う。task blastn の既存分類を適用しない |
| Sakai/MG1655 megablast | 6,476 行中 5 行。既存監査・製品判断に common-endpoint の source-underdetermined 分類がある | NZ と同じ原因とは扱わない。既存分類と現在の状態列を照合し、移植漏れとソースが一意に規定しない選択を分ける |
| TBLASTX | cold Node/WASI 比較に限る TurboFan 設定を採用済み | 維持する。BLASTN・BLASTP の実績に数えず、専用メモリ許容も転用しない |

根拠: [BLASTN 実験](evidence/wasm_performance_20260913/run-20260913-01/P2.md)、[BLASTP 実験](evidence/wasm_performance_20260913/run-20260913-01/P3.md)、[追加報告](evidence/wasm_performance_20260913/run-20260913-02/REPORT.md)、[megablast 監査](evidence/wasm_performance_20260913/run-20260913-01/audit-p4.md)。

実験には `ACCEPTED / REJECTED / INCONCLUSIVE`、成果目標には `NOT_STARTED / IN_PROGRESS / UNMET / BLOCKED / ACHIEVED` を使う。候補を棄却しても G-N/G-P は終了しない。次の調査対象・不足する証拠・再開条件を記録する。旧 P2/P3 の `COMPLETE` は実験終了の歴史的記録として保持し、本計画の目標達成へ引き継がない。

## 3. 守る契約

- [AGENTS.md](../AGENTS.md) と [verify-ncbi-parity-and-speed](../.agents/skills/verify-ncbi-parity-and-speed/SKILL.md) を適用する。NCBI C/C++ を唯一の挙動根拠とし、変更した Rust と移植関数のテストには対応するファイル・行番号・snippet を直上に記載する。
- NCBI executable は比較 oracle と診断専用。runtime、build、FFI、未実装機能の fallback に使わない。参照元にない機能や tie-break を追加しない。
- seed/candidate 集合、diagonal 状態、X-drop、HSP 構築・剪定・順序、統計・精度・丸め・書式を維持する。CBS/SEG の省略、近似 DP、独自の候補削減で高速化しない。
- raw bytes の比較を先に行う。構造化 diff、座標による対応付け、行数、再計算スコアは診断専用であり、並べ替え・header 除去・数値正規化を合格条件にしない。
- plain `wasm32-wasip1` は serial。並列は `wasm32-wasip1-threads` と `wasm-threads`。明示されたスレッド数を黙って縮小しない。search-scoped pool、終了待ち、エラー伝播、元の順序への結果復元を維持する。
- parity 修正と性能変更は別の候補・差分として扱う。原因に関係する既知差が残る間は、その追跡に必要な実行に限定し、無関係な性能測定や広いテストを増やさない。
- 既存の凍結出力、platform authority、旧 evidence、無関係な作業ツリー変更を上書きしない。公開・push・tag・deploy は本計画に含めない。

### Sakai の分類と凍結基準の扱い

[PD-BLASTN-HSP-CANONICALIZATION](product_decisions/PD-BLASTN-HSP-CANONICALIZATION.md) と [既存 footprint](../LOSAT/tests/blastn_v010_source_exceptions.tsv) は、NCBI comparator が等価と扱う異なる edit script の選択について、限定された既存分類を記録している。本計画はそれを削除・拡張せず、G-M の raw 一致を達成したことにも読み替えない。

NCBI ソースが一意な survivor を定めておらず、既存契約のまま特定バイナリへ一致させる根拠が見つからない場合、G-M の該当部分を `BLOCKED` とする。必要なのは同じ入力状態、comparator の全比較結果、sort 前後の HSP 列、選ばれた script を結び付けた証拠である。単に原因不明、調査が長い、候補が不採用という理由では `BLOCKED` にしない。独自の gap/identity tie-break、プラットフォーム分岐、出力書き換えで解決しない。

[PD-NCBI-PLATFORM-VARIANCE](product_decisions/PD-NCBI-PLATFORM-VARIANCE.md) の Gate A（凍結 PR5 の LOSAT raw bytes）と Gate B（登録済み native NCBI fingerprint）は別々に確認する。正当な修正が凍結 LOSAT 出力を変える場合も、旧 Gate A の失敗を残す。影響する case、旧・新 bytes、NCBI 根拠をまとめ、既存方針に従った新しい認証基準のレビューを別途必要とする。golden の自動更新や未知 fingerprint の自動登録は行わない。契約判断が必要になっても、影響しない G-N/G-P の作業は継続できる。

## 4. 実行順序と成果物

通常の依存順は `R0 → R1 → R2 → R3 → R4 → R5 → R6`。測定・build は競合させず、実装候補は一度に一つとする。

| 段階 | 作業 | 段階を終了するための証拠 |
|---|---|---|
| R0 | 現行入力・source・環境と基準を固定 | 新しい manifest、source snapshot、独立 build、再現 argv |
| R1 | megablast の最初の不一致を追跡 | 各差分の最初の不一致と直前の一致状態、NCBI/Rust owner |
| R2 | 根拠のある parity 修正と対象範囲の sweep | 修正前後の raw diff、境界テスト、全既知差の disposition、独立監査 |
| R3 | BLASTN の費用を測り、原因に対応する候補を実装 | 原因証拠、NCBI 同値性、採用 A/B・対照。未達なら R3 内で調査を続ける |
| R4 | BLASTP の費用を測り、原因に対応する候補を実装 | 原因証拠、NCBI 同値性、採用 A/B・対照。未達なら R4 内で調査を続ける |
| R5 | 採用変更の統合・必要な回帰・実測 | 最終 source での正しさ、各性能目標、thread/lifecycle/memory の結果 |
| R6 | 独立確認と引き渡し | G-M/G-N/G-P 別の達成状態、残件、再現手順、差分レビュー |

### R0 — 現在の状態を固定する

1. `git status` と対象ファイルを確認し、実装・テスト・文書・生成物を区別して初期差分を保存する。旧記録の artifact と現在の source が同じだと仮定しない。
2. 証拠を `docs/evidence/megablast_wasm_np_remediation_20260913/<new-run-id>/` に新規保存する。入力、source/build config、runner、Cargo.lock、toolchain、Node/V8、NCBI executable/source の版と SHA-256 を記録する。既存 run ID を再利用しない。
3. 同一 source から native command、serial WASI command、threaded WASI command を別 target directory で release build する。reactor の変更・検証が必要な段階では reactor も別に build する。NCBI 診断用 build は公式 oracle から分離する。
4. 比較器は既存 runner を読み、現在の argv/schema を確認して使用する。欠落や誤計測が確認された場合だけ最小限修正し、汎用基盤を作り直さない。
5. 基準を `B0=開始時実装`、`B1=parity 修正後の実装`、`C=性能候補` と命名する。B0→B1 は正しさの差分、B1→C は出力不変の性能差分として記録する。既知の誤出力を速度の合格基準にしない。限定された契約 blocker が残る場合は、B1 の合格 fixture と未合格 fixture を明示する。

### R1 — megablast の最初の不一致を特定する

対象は三入力。`task=megablast`、遺伝暗号既定、初期比較は同一条件の outfmt 6 とする。reward/penalty、gap、word size、dust/soft masking、strand、evalue、hitlist 制限を含む全 argv を旧 run と照合して固定する。outfmt 7 の manifest を無言で outfmt 6 の再現条件へ代用しない。

| 入力 | 役割 | 追跡の出発点 |
|---|---|---|
| NZ_CP006932.fasta self | 主な未解決差 | 旧記録の query 320597–321685 / subject 321850–322940。現在も同じ差なら追跡キーに使う |
| Sakai.fna / MG1655.fna | 既存 source-underdetermined footprint | 旧 5 座標キーを現在の結果と対応付ける。行番号を HSP 識別子にしない |
| EDL933.fna / Sakai.fna | 合格対照 | 現行 oracle との raw 一致を確認する |

1. まず native n1 と固定公式 oracle で各入力を新しい出力先へ再実行する。`cmp` による raw 比較後、既存 `compare_blastn_parity.py` で全フィールドと行順の差を保存する。差が変わった場合は旧説明を当てはめず、source・argv・入力の境界を調べる。
2. 同じ artifact/source・入力で native n8、serial n1、threaded n1/n8 を比較する。スレッド差があればその最初の不一致を先に切り分ける。Wasm 固有でない差に対して JIT を原因と決めない。
3. 差がある HSP について、seed、context、内部 query/subject offset、raw score、gapped start、候補の生成順、preliminary HSP 配列を記録する。出力座標だけでなく内部状態で対応付ける。
4. forward/reverse greedy extension、traceback の predecessor 選択、左右 script の結合、`s_ReduceGaps` 相当、再評価、common-start/common-end purge、最終 hitlist、identity/length 計算、formatter の順に比較する。各段階の入力と出力を残し、「初めて違う段階」と「最後に一致した段階」を一組で示す。
5. 同点候補は raw score と全 comparator key、edit script、sort 前後の列、survivor を記録する。identity/gap の追加 tie-break は提案しない。NCBI source build の診断結果は、同じ入力の未計装結果・公式 oracle との関係を確認し、公式版を置き換えない。
6. 既存の `LOSAT_TRACE_BLASTN_*` と greedy/purge trace を優先する。不足する状態だけを隔離診断 build で採取し、ログによる時刻・順序への影響を記録する。NCBI の計装が必要なら writable な隔離 source copy を使い、参照 tree を変更しない。
7. fixture の縮小は、同じ最初の不一致を保つと確認できた場合だけ行う。部分配列化で seed、mask、統計、候補順が変わって差が消えた場合は、その縮小例を修正根拠にしない。

成果物は case ごとの `first-divergence.md` と raw/structured diff、状態 snapshot、正確な再実行コマンド。NZ の原因が分かっても Sakai へ一般化せず、全差分を原因別に対応付ける。

### R2 — parity 修正を完了させる

1. 最初の不一致を起こす NCBI の caller、入力、呼び出し順、比較演算、更新タイミングを確認する。対応する Rust owner の最小箇所を修正する。
2. 同じ欠陥の全 call site を確認し、一つの座標だけに効く分岐を作らない。NZ/Sakai の全差分を説明し、同じ原因の修正を一巡させてからテストへ進む。
3. NCBI 移植関数の境界テストを追加する。実際の原因に応じて、等価 HSP、gap tie、負向き、空/短配列、sentinel、左右結合、gap reduction、再評価、purge の境界を選ぶ。可能なら NCBI unit test を参照する。
4. 三入力の raw 比較、serial/threaded/native の一致、該当 BLASTN task suite を順に実行する。修正が formatter を含む場合は 0/6/7/custom を同一契約で比較する。限定 source-underdetermined blocker が残る場合は、必要な非回帰とその未合格結果を別々に保存する。
5. `ncbi_parity_auditor` に実装・NCBI ソース・保存した実出力の独立 read-only 確認を依頼し、原因と修正の因果関係を確認する。raw 差が残る目標を解決済みにしない。
6. Gate A/B への影響を早期に確認する。契約変更を要する場合は §3 の手続きへ分岐し、旧基準を維持したまま修正候補をレビュー可能にする。

### R3 — Wasm BLASTN の全体時間を減らす

主対象は LC738874/LC738870 と AP027202/LC738875 の `task=blastn`。小入力 AP027152/AP027202、高密度 NZ self、fallback を通った LC738873/LC738871、および修正後 megablast 三入力を対照とする。

1. B1 の native n1/n8、serial n1、threaded n1/n2/n4/n8 を測定する。cold 全体時間と、起動・検索各段階・出力・終了の費用を分ける。計測区間が重なるものは足し合わせない。
2. 既存のバッチ数・計算数・消費数に加え、バッチ内の仕事時間分布、最後の slot 完了までの待ち、順序付き replay、preliminary/seed、reserve/realloc、実際の copy/初期化量を必要な範囲だけ計測する。scratch slot を worker ID と呼ばず、容量を allocation 回数や RSS に読み替えない。
3. 旧診断の 0/1-job batch と scratch constructor だけでは大きな改善を説明できなかった事実を使って、候補を優先順位付けする。捨てられた DP の件数比を時間比と扱わない。
4. 候補ごとに「全体のどの直列区間または並列完了待ちを何秒減らせるか」「追加費用」「NCBI と同じ計算を維持する根拠」を先に記述する。5%に届く見込みを根拠付きで評価できなければ、未測定部分の調査へ戻る。
5. 原因が負荷の偏りなら既存 pool/scratch 内の同じ仕事集合の割当、同期なら限定したバッチ設計、初期化・コピーなら実測された重複処理、単一スレッド kernel なら該当 hot path を一つずつ検討する。バッチ値や fixture 名に依存した選択を先に実装しない。
6. 先行 DP は純粋な計算に限定し、prelim index 順の containment、start adjustment、traceback 結果採否、identity 判定、tree 挿入を維持する。megablast の greedy 経路に DP 調整の効果を主張しない。
7. §6 の採用手順を実施する。不採用なら棄却理由で仮説を更新し、次の最も大きい未説明費用へ戻る。G-N は `IN_PROGRESS` または `UNMET` のままにする。

### R4 — Wasm BLASTP の全体時間を減らす

主対象は AP027078/AP027131 と AP027132/NZ_CP006932 の `.faa`。AP027131/NZ_CP006932 を長入力・正常終了対照、WSSV/PajaWSV と SicyWSV/CoBV を小入力対照とする。

1. B1 による native n1/n8、serial n1、threaded n1/n2/n4/n8 の検索と段階時間を再取得する。通常 Wasm build で `LOSAT_TIMING=1` だけでは必要な段階 timer が動かないため、実際の計装有効性を確認する。
2. main と全 worker のプロファイルを取得する。code load/move/delete と関数対応を確認し、未対応サンプルを残す。完全な対応付けが得られなければ、粗い段階 timer と call/cell/byte counter で不足部分を測る。未対応サンプルを既知関数へ比例配分せず、kernel の tick 比を全体時間比にしない。
3. query/subject 準備、scan/ungapped、preliminary DP、Kappa query redo・match redo、score-only DP、traceback、merge/format を区別する。cell 数、band 幅、再計算、row 確保・初期化・copy、matrix 参照、cache lookup/clone/insertion を必要に応じて測る。
4. SEG cache はヒット数増加だけで選ばない。検索全体の時間を減らせる費用が示されない限り、棄却済み候補を再実装しない。`map_init` / `for_each_init` の初期化数を worker 数と仮定しない。
5. 全体時間に十分な寄与がある一つの処理について、既存 row/scratch 配置、余分なコピー・初期化・dispatch、query/subject 処理の分担を検討する。既存 SIMD、BLOSUM62 特殊化、Kappa 並列化を未実装とみなして作り直さない。
6. 同じ cell、処理順、X-drop、gap tie、script、CBS/SEG、adjusted matrix、浮動小数点精度・丸め、restricted→exact retry を維持する。未対応の CLI オプションを対照のためだけに実装せず、既存 fail-fast を維持する。
7. 標準/adjusted matrix、forward/reverse、band/fence、scratch 再利用、near-identical/SEG、restricted→exact retry、single-query/many-subject、many-query のうち変更分岐を実際に通る対照を用意する。retry counter が 0 の run を retry 検証と呼ばない。
8. §6 の採用手順を実施する。不採用時は原因の調査へ戻り、G-P を完了にしない。

## 5. NCBI と LOSAT の調査箇所

NCBI root は `/mnt/c/Users/genom/GitHub/ncbi-blast/`、代替は AGENTS.md の指定を参照する。以下の行番号は既存監査と計画時のソース確認による入口であり、実装時に現在の行番号・caller を再確認する。コメントの説明だけでなく実際の比較式を読む。

| 境界 | NCBI の入口 | LOSAT の入口 |
|---|---|---|
| greedy predecessor 選択 | `c++/src/algo/blast/core/greedy_align.c:278–294`、`s_GetNextNonAffineTback` | [greedy.rs](../LOSAT/src/algorithm/blastn/alignment/greedy.rs) の `get_next_non_affine_tback` |
| 左右 extension、script 結合、gap reduction | `blast_gapalign.c:2669–2936`、`s_ReduceGaps` / `BLAST_GreedyGappedAlignment` | 同 `greedy_gapped_alignment_internal` / `reduce_gaps` |
| common-endpoint comparator と survivor | `blast_hits.c:2268–2387,2455–2535` | [purge_endpoints.rs](../LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs)、[hsp.rs](../LOSAT/src/algorithm/blastn/hsp.rs)、preliminary 選択 |
| containment、traceback、identity、tree 更新 | `blast_traceback.c:403–405,436–472,509–512,583–609`、`blast_itree.c` | [BLASTN run.rs](../LOSAT/src/algorithm/blastn/blast_engine/run.rs)、[gapped.rs](../LOSAT/src/algorithm/blastn/alignment/gapped.rs) |
| protein DP / traceback | `blast_gapalign.c:374` 以降の `ALIGN_EX`、特に `531–726` | [BLASTP gapalign.rs](../LOSAT/src/algorithm/blastp/gapalign.rs) の score-only / traceback 実装 |
| Kappa redo、SEG 条件、局所状態 | `blast_kappa.c:1414–1454,1626–1643,2942` 以降、`3493–3503`、`redo_alignment.c` | [blast_engine.rs](../LOSAT/src/algorithm/blastp/blast_engine.rs)、[kappa.rs](../LOSAT/src/algorithm/blastp/kappa.rs) |

## 6. 性能候補の採用・棄却・継続条件

### 固定する測定条件

- B1 と候補は同じ入力、検索 argv、出力先の種類、Node/V8、build profile/features、thread 数で比較する。NCBI は速度目標ではなく出力 oracle とする。native 対 Wasm の時間比は診断値として別に示す。
- BLASTN/BLASTP の Node 追加 flags は baseline/candidate とも同一に固定する。TBLASTX 専用 flags と追加メモリ許容は両プログラムに適用しない。
- 探索は条件ごとに 3 反復。採用候補は warmup 1 回を除外し、A/B・B/A を交互にした 5 反復以上の未計装 release を測る。全サンプル・範囲・中央値を保存する。曖昧なら追加 5 回、それでも不明なら `INCONCLUSIVE` とする。
- cold は起動・検査・compile・worker・検索・出力・終了を含む。module 再利用＋新 command instance、同一 reactor instance の初回/以降は別表にする。再利用の準備費用を除外して cold の改善と呼ばない。
- CLOCK_MONOTONIC を同じ境界で採取し、利用可能なら BOOTTIME と照合する。realtime の不一致は記録する。診断ログによる大幅な時間増加を採用用試料へ混ぜない。
- CPU/OS、CPU 制約、filesystem、CPU user/system、process peak RSS、Wasm linear memory を記録する。CPU/wall 比を worker 稼働率とみなさず、RSS と linear memory と scratch 容量を区別する。
- timeout は短中入力 300 秒、長入力 3,600 秒を出発点に、候補の結果を見る前に十分な共通上限を固定する。失敗・timeout・欠落出力は成功時間に変換せず、原因調査対象として全件残す。

### 採用ゲート

| 項目 | 条件 |
|---|---|
| 正しさ | 修正後基準と候補の raw bytes が一致し、該当する公式 oracle 契約も合格。既存未合格差を隠さない |
| 主対象 | G-N/G-P それぞれの二つの主入力で、threaded n8 の cold 中央値を各 5%以上短縮 |
| 非回帰 | 小入力、native n1/n8、serial n1、threaded n1、影響する他 task で、時間増加が `max(基準中央値の5%, 50 ms)` 以下 |
| メモリ | 同条件の測定 run の最大 peak RSS の増加が `max(基準の10%, 16 MiB)` 以下。Wasm memory/scratch 容量も記録し、反復での継続増加・予算超過を認めない |
| 並列性 | 同一 threaded artifact の n1/n2/n4/n8 で raw・要求 pool 数を確認。診断 run で各 worker の spawn→ready→exit と終了待ちを検証 |
| 適用範囲 | 変更分岐を通る境界・対照で非回帰。unrun の行列/retry/形式/反復を PASS としない |

候補が失敗したら、その候補の採用を止める。性能目標は継続し、次の反証可能な仮説を立てる。局所的な小改善を合わせる場合も、各要素が出力・非回帰を満たす独立した根拠を持ち、組合せを新しい候補として全体で検証する。未検証・悪化した候補を積み上げない。

主入力の追加や主採用軸の変更は理由を先に記録し、元の入力と目標の結果も残す。測定後に主対象や 5% 基準を下げて達成にしない。データや実行環境が不足する場合は必要な入力を特定する。性能上限は直列区間や並列完了待ちを含む実測に基づいて説明し、数件の棄却だけで「Wasm の限界」と結論しない。

## 7. R5/R6 — 統合、回帰、独立確認

1. 採用済みの変更だけを順に統合する。B0→parity 修正、B1→各候補、B1→最終版を区別し、最終 source から対象 artifact を再 build する。
2. parity sweep 後、変更した関数の focused tests、該当 task suite、共有箇所に必要な回帰を実施する。Rust を変更した場合は `cargo test --all-features`、`cargo clippy --all-features --all-targets -- -D warnings`、`cargo fmt --check` と両 WASI target の release build を実施する。
3. outfmt 0/6/7/custom、threaded n1/n2/n4/n8、native、serial について変更の影響範囲を確認する。pool/host/API/scratch 寿命を変更した場合は command/reactor を含む反復 `1→2→4→8→2→1`、engine 切替、部分 spawn 失敗、trap、次の呼び出しへの影響まで既存 runner で検証する。
4. TBLASTX の LC738874/LC738875 evalue 10/100/10000 と長入力 AP027131/AP027133 gencode 4 を保護する。共有検索処理や runtime に影響する変更では回帰を実行する。非既定 local-subject `db_gencode` は既存の DB oracle 等の契約を守る。
5. `ncbi_parity_auditor` が最終ソース、NCBI の挙動根拠、実出力、全サンプル、source/artifact/argv の対応、採用・除外条件を独立 read-only 確認する。指摘を直した source の測定を旧 artifact の証拠で代用しない。
6. 過去の固定回帰 99/145、未完了 46 件は当時の記録として維持する。以前のユーザー指示による延期を勝手に解除せず、今回必要な対象回帰と正式認証の残件を区別する。新しい runtime/source の成果を過去 99 件で代用しない。正式な全 native platform 認証は別の到達点である。
7. production、test/計測、文書、生成物の差分を分けてレビューする。実施した範囲、未達目標、契約 blocker、unsupported scope を最終報告する。実装時の commit handoff は最終差分に合う英語タイトルと短い説明を作成する。

## 8. 既存の実行入口と証拠の保存形式

新しい runner を増やす前に次を使う。引数は実行時に source を読み直し、全 argv を manifest に保存する。

| 用途 | 入口と注意 |
|---|---|
| 入力と検索引数 | [comparison_cases.tsv](../LOSAT/tests/comparison_cases.tsv)、[blastn_parity_manifest.tsv](../LOSAT/tests/blastn_parity_manifest.tsv)、[blastp_v010_parity_manifest.tsv](../LOSAT/tests/blastp_v010_parity_manifest.tsv)。case ID と `losat_stem` は区別する |
| BLASTN fresh paired / diff | [compare_blastn_parity.py](../LOSAT/tests/compare_blastn_parity.py)。`--fresh-paired --paired-output-dir`、直接 diff は `--ncbi --losat --fail-on-byte-diff`。構造化判定だけで raw PASS にしない |
| NZ megablast の登録 | comparison catalog には存在するが現行 BLASTN parity manifest には task blastn の NZ 行だけがある。R0 で新 run 専用の manifest 行を作り、固定 release manifest を変更せず megablast の argv を明示する |
| Wasm build | [build_wasi_artifacts.py](../LOSAT/tests/build_wasi_artifacts.py) と独立 `CARGO_TARGET_DIR`。command/reactor、serial/threaded の種類を検証する |
| A/B benchmark | [benchmark_wasm_threading.py](../LOSAT/tests/benchmark_wasm_threading.py)。candidate/baseline artifact と runner、oracle、fresh output directory を明示し、`--case` は catalog の exact `losat_stem` を使う |
| thread / lifecycle | [check_wasm_threading.py](../LOSAT/tests/check_wasm_threading.py)、[check_wasm_threading_regressions.py](../LOSAT/tests/check_wasm_threading_regressions.py)。既存の正式認証と focused 回帰を区別する |
| BLASTP 書式・threads | [run_blastp_threads_comparison.sh](../LOSAT/tests/run_blastp_threads_comparison.sh)、[run_blastp_pairwise_threads_comparison.sh](../LOSAT/tests/run_blastp_pairwise_threads_comparison.sh)、[run_blastp_outfmt7_comparison.sh](../LOSAT/tests/run_blastp_outfmt7_comparison.sh)、[run_blastp_custom_fields_comparison.sh](../LOSAT/tests/run_blastp_custom_fields_comparison.sh)。実行前に既存出力の上書き先を確認する |

主対象の `--case` は BLASTN が `MelaMJNV.PemoMJNVA.losatn.blastn` と `MjPMNV.MlPMNV.losatn.blastn`、BLASTP が `AP027078.AP027131.losatp` と `AP027132.NZ_CP006932.losatp`。表示名から FASTA 名や task を推測しない。

各 run は [evidence template](../.agents/skills/verify-ncbi-parity-and-speed/references/evidence.md) を満たし、最低限、次を対応付けて保存する。

- `manifest.json`: program/task、全 argv、入力/source/artifact/runner hash、実行環境、Node 設定、固定した測定方針。
- `parity/`: raw 出力、raw diff、構造化 diff、状態 snapshot、各差分の NCBI/Rust owner、最初の不一致。
- `candidate-<id>/`: 独立 patch/build log、仮説、同値性の説明、全 measured/diagnostic sample、失敗・除外理由、採否。
- `STATUS.md` と `REPORT.md`: 成果目標と実験状態を別記し、次の具体的な調査と不足証拠を記載。
- `SHA256SUMS`: 保存した実ファイルの一覧。大きな raw/profile は archive の場所と hash を記録し、Git にある範囲とローカル保存範囲を明示。

## 9. 開始時の目標台帳

| 目標 | 状態 | 次の作業 | 達成を妨げ得る未確認事項 |
|---|---|---|---|
| G-M | NOT_STARTED | R0/R1: 現行 megablast 三入力を再現し、NZ の最初の不一致を追跡 | NZ の原因、Sakai の現在の survivor 状態、凍結 Gate A への影響 |
| G-N | NOT_STARTED | R2 の修正後基準または限定 blocker の切り分け後に R3 | 検索全体の未説明費用、負荷の偏り、copy/初期化の実費用 |
| G-P | NOT_STARTED | R4: Wasm の段階時間と kernel の費用を対応付ける | 未対応 profile、DP/Kappa の改善可能量、主入力間の差 |

三つとも `ACHIEVED` になった場合にのみ、この修正計画の成果目標を完了とする。一部が `UNMET/BLOCKED` なら、完了した実装・測定の範囲を示しつつ、計画全体は未完了として引き渡す。

## 10. ユーザー指示による暫定採用（2026-09-14）

「計測はこれ以上繰り返さないで。まとめよう」「とりあえず採用しようよ」の指示により、現状の実装と検証済みnative設定を暫定採用し、追加計測を停止した。全589対象sourceの一致を確認。G-Mの現行比較は達成、G-N/G-Pの主入力速度は短いBLASTNの既承認免除を含め条件を満たす。一方、TBLASTX反復時間9/12未達、後続6群と最終native megablast再確認は未実行のまま保存する。G-N/G-Pの最終統合ゲートおよび計画全体を全条件達成と読み替えない。詳細は上記STATUSとuser-provisional-adoption-20260914.jsonを参照。
