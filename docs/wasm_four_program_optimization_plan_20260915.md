# megablast / blastn / tblastx / blastp Wasm最適化実装計画

作成日: 2026-09-15

状態: **実装済み・統合版を採用（2026-09-15、ユーザー指示「じゃあいいじゃん。採用」）。測定は最初の3組で終了。**

採用記録: [ADOPTION.md](evidence/wasm_four_program_20260915/run-01/ADOPTION.md)。統合版の時間基準は39/45条件、RSS基準は45/45条件で合格。速度低下を含む結果を確認したうえで採用し、未実施の測定や既存課題は記録を維持する。

最終変更: 追加測定を中止し、body/reuse/TurboFan等は未実施として記録する。以下の5組の記述は当初計画であり、今回の最終集計には3組を用いる。

実装記録: [2026-09-15 run-01](evidence/wasm_four_program_20260915/run-01/REPORT.md)。以下の計画作成時の記述と、今回の検証結果は区別する。

## 1. 目的と優先順位

NCBI BLASTの検索結果を保ち、LOSATの重複計算、コピー、割り当て、同じ状態の再処理を減らす。WasmとNativeの速度差は、同じ入力・実スレッド数・計測区間で評価する。Native同等の達成は現時点では保証しない。

| 対象 | 最初の実装候補 | 次の候補 | 難易度・主要リスク |
|---|---|---|---|
| `blastn -task megablast` | 非アフィンgreedy行の往復コピー除去 | 最初の不一致を探す整数比較の高速化 | 中〜高。scratchの残存内容、再試行、境界 |
| `blastn -task blastn` | 並列traceback開始位置の再計算除去 | バッチ用配列再利用、queryの4塩基変換再利用 | 低〜中。包含判定の時点とquery encoding |
| `tblastx` | 最大値管理木の変化しない更新伝播を止める | 実測で選ぶリンク判定ループのメモリアクセス改善 | 中〜高。同点順位、リンク状態、浮動小数点順序 |
| `blastp` | 不変の初期行列生成を検索内で再利用 | 単一query準備再利用、diagonal初期化削減、先行redoの有界化 | 低〜高。経路によって適用範囲が異なる |

推奨順序は、共通baseline確定 → N1 → P1 → M1 → X1 → 各対象の次候補。X0の診断で木の費用が小さければX1を見送り、X2の調査へ進む。番号順の全候補採用を目的にしない。

## 2. 現在地と適用範囲

ソース読解時のLOSAT HEADは`1de348ad`。作業ツリーには多数の既存変更があり、HEADだけでは実装を識別できない。NCBI参照ツリーは`/mnt/c/Users/genom/GitHub/ncbi-blast/`、HEADは`598d8ae6a72b923127ba2fbfaffd48e4c83bfb4`。以下の行番号はこの読解時点の目印で、実装時に再照合する。

すでに入っている変更を新規候補として数えない。

- BLASTN endpoint pruningの固定HSP格納と小さなハンドルの並べ替え。
- BLASTP traceback行の予約済み領域への書き込みと行単位の長さ更新。
- TBLASTX最終リンクソートの参照利用、最大値管理木、既存SIMD。
- Wasmの`simd128`ビルド設定、compiled moduleのworkerへの受け渡し。
- `num_threads=N`に呼び出し元を含め、子スレッドをN−1本作る現行契約。
- Node比較用のTBLASTX専用コンパイルフラグ。これはブラウザの通常設定とは別の実行条件。

過去の不採用候補をそのまま復活させない。BLASTNの固定DP band変更・batch32、BLASTPのSEG cache・単純な行push変更・行列dispatch・関数境界変更、TBLASTXのループ抽出・不要field除去には不採用または不確定の記録がある。再調査には、新しいprofileや異なる変更点を明記する。

根拠: [直近の統合報告](evidence/megablast_wasm_np_remediation_20260913/run-20260913-01/REPORT.md)、[調査履歴](evidence/megablast_wasm_np_remediation_20260913/run-20260913-01/performance-investigation.md)、[現在のthread契約](wasm_total_threads_20260914.md)。履歴内の途中経過より、統合報告と現行ソースを優先する。

## 3. 共通の開始条件・計測・採用条件

### S0: 現在のbaselineを固定する

1. 既存変更を保持し、tracked/untrackedの必要なbuild入力・fixture・runnerをhash付きで保存する。baseline/candidateは別の出力先へrelease buildする。実装、計測コード、文書、生成物を別々にレビューする。
2. 各fixtureのprogram/task、FASTA hash、全argv、遺伝暗号、出力形式、thread数、NCBI binary/version/hashを記録する。古いmanifestの分類や保存出力を現在の結果と見なさない。
3. 実際の利用可能な各targetでbaselineの現在の比較を実行し、対象経路を通ることを確認する。新しい差があればNCBIの最初の相違箇所を特定し、性能変更と分けて扱う。該当経路の差が未解決なら、その性能候補の採用へ進まない。
4. Native、serial Wasm、threaded Wasmを区別する。threadedは同一artifactのn1/n2/n4/n8で比較し、n8が呼び出し元1＋子7であることを確認する。plain `wasm32-wasip1`に並列速度向上を期待しない。
5. Node/WASIとブラウザの結果を分ける。双方を対象とするが、主となるruntime・入力・cold/反復の評価区間は候補の測定前に固定する。ブラウザ未実行なら、その範囲は未確認と報告する。

### 予定fixture

以下は開始時に実在・経路・raw比較を確認する候補。経路を通らないfixtureは実装前に理由付きで差し替え、測定後に都合のよい入力へ変更しない。

| 対象 | 主fixture候補 | 保護する対照 |
|---|---|---|
| megablast | EDL933/Sakai、NZ_CP006932 selfを明示的に`-task megablast`で実行 | Sakai/MG1655、短配列、no-hit、非ゼロgap設定のgreedy別経路 |
| task blastn | LC738874/LC738870、AP027202/LC738875 | compact multi-query、NZ_CP006932 self、word size 7/11、負strand、no-hit |
| tblastx | MjeNMV/MelaMJNV、AP027280 self | MelaMJNV/PemoMJNVA、LC738874/LC738875のE-value 10/100/10000、AP027131/AP027133のcode 4、no-hit |
| blastp | AP027078/AP027131、AP027132/NZ_CP006932 | SicyWSV/CoBV、AP027131/NZ_CP006932、単一query＋複数subject、no-hit、対応する行列・組成補正設定 |

既存の`comparison_cases.tsv`にある名称と実際のFASTA名の対応を保存する。特にBLASTPの多query fixtureだけで、単一queryの最適化を評価しない。単一query候補には別の主fixtureを実装前に固定する。

### 計測方式

- 別の診断buildでcall数、処理cell数、コピーbyte数、allocation数、各段階の時間を取る。cellごとの時計読み出しや共有atomic更新を入れず、slot内で集計して段階終了時に回収する。
- 採用時間は診断なしのrelease buildで測る。baseline/candidateの実行順を交互にし、同じCPU条件、同じrunner実パス・ファイルシステム、同じ出力先の方式を使う。Wasm両artifactも同じ種類のファイルシステムに置く。
- 一つの条件につきwarmup 1回＋測定5組を基本とし、中央値と全sampleを保存する。結果が不確定なら、事前に定めた追加5組を一度だけ実施して全組を評価する。最速値や好都合なcohortだけを選ばない。
- coldは起動・コンパイル・準備・検索・出力を含む。検索本体、compiled module再利用、同一reactor反復は別の表にする。前計算やcopyを計測区間の外へ移して短縮と報告しない。
- wall timeと取得可能なCPU time、RSS、linear memory、保持allocation/capacityを区別する。並列slotの時間合計をwall時間の内訳として足し合わせない。
- TBLASTXは通常Node設定と既存TurboFan専用設定を別条件にし、各条件内で両版のflagsを一致させる。ブラウザの通常設定の結果をNodeのflags付き結果で代用しない。

### 採用判定の提案

新しい実装ラウンドの基準案であり、過去の承認や判定を書き換えない。

1. 対象のraw output、必要な内部状態・順序・浮動小数点bitsが一致すること。NCBIは外部の検証oracleだけに使う。
2. 候補が対象とする事前指定の各主fixtureで、全体時間の中央値5%以上の短縮を目標とする。小さい局所改善は正しさと非回帰を満たした部品として保留できるが、対象全体の目標達成とはしない。保留部品の統合は新しい全体測定を要する。
3. 各対照の悪化は時間`max(5%, 50 ms)`、RSS`max(10%, 16 MiB)`以内を基準案とする。linear memoryの既存上限、API budgetも維持する。
4. `Wasm時間 / 同一版Native時間`をn1とn8で示す。「Native並み」は事前指定の主fixture・runtimeにおける検索本体の比率1.20以下を目標案とし、coldと反復の実利用時間も併記する。これは達成予測ではない。両者の絶対時間も示し、Nativeの悪化による比率改善を認めない。
5. 出力差、trap、thread契約違反、候補由来のtimeoutは不採用。未実行・不確定を合格にしない。各候補を単独で評価してから統合する。

### parityと既存課題の境界

- 現行[AGENTS.md](../AGENTS.md)のPR5 Gate Aと登録済みNCBI platform fingerprintのGate Bを維持する。platform-local NCBI値で凍結expectedを更新しない。
- megablastの修正済み挙動と凍結Sakai契約との差は、現行sourceの比較と正式release認証を分けて記録する。baseline一致だけでNCBI互換を宣言しない。
- TBLASTXの非既定local-subject `db_gencode`だけは既存例外に従う。code 4回帰では対応するNCBI DB oracleを使い、それ以外の差を許容しない。
- 既存のTBLASTX reuse時間の失敗と[reactor memory課題](wasm_reactor_memory_followup_20260914.md)は維持する。新候補について固定回数のbaseline/candidate比較で追加悪化を確認する。観測期間を伸ばし続けて合格にしない。
- 全コード変更に実在NCBI sourceのpath・行番号・snippetを付ける。pruning・ordering・浮動小数点変更とrelease向け主張は、採用前に`ncbi_parity_auditor`の独立read-only確認を受ける。

## 4. megablast実装計画

### M0: コピーが実際に発生する範囲を確認

対象: [greedy.rs](../LOSAT/src/algorithm/blastn/alignment/greedy.rs)の`blast_greedy_align`、`NonAffineGreedyRow`、`GreedyNonAffineMem`。

非アフィン経路ではbase行・pool行の`to_vec()`と、終了時の`persist()`による書き戻しがある。主fixtureで各経路の呼出数、行数、確保数、コピーbytes、最大pool長、再試行数を測る。非ゼロgapのaffine経路へ効果を一般化しない。

NCBI owner: `core/greedy_align.c:71–77,100–125,666–678`のpool rewind、領域取得、行ポインタ設定。呼出側のscratch寿命も確認する。

### M1: 行の値を既存scratchに置いたまま計算する

1. 行descriptorにorigin、長さ、base/pool内の位置を保持し、行ごとの所有`Vec<i32>`をなくす。既存の二つのbase行とtraceback poolを唯一の値の所有者とする。
2. `get/set`をその領域へ向け、pool→Vec→poolの往復コピーを削除する。pool拡張をまたぐ参照を保持せず、整数位置から必要時に参照する。
3. rewindは使用位置だけを戻し、NCBIが残す既存cell値を保持する。score-onlyの二行交替、tracebackでの行保存、追加行の確保時点を維持する。
4. 正常終了、非収束、fence、再試行の全returnを調べ、次の呼出しから見えるscratch内容がbaseline/NCBIの契約と一致することを確認する。単なる`persist`削除で終えない。

検証: fresh/reused scratchの連続呼出し、forward/reverse、score-only/traceback、空・極短配列、非収束→拡張retry、fence、長さとdiagonal幅の境界。score・座標・edit scriptに加え、次回から参照可能なscratch内容を比較する。非ゼロgap経路も非回帰確認する。

完了条件: 行別allocationと対象の往復copyが消え、raw一致と全体時間・memory条件を満たす。領域を大きく先取りすることでallocation数だけを減らした結果は別途memory評価する。

### M2: 最初の不一致を探す処理を高速化する（条件付き）

M1後も`find_first_mismatch_greedy`が主要費用なら、最初は非圧縮subjectの分岐だけで、複数塩基の整数比較またはWasm SIMDを試す。最初の不一致・曖昧塩基・fenceで正確に止まり、範囲内の完全なblockのみ読む。端数は既存scalar処理で扱う。圧縮subjectの位相処理は別候補とし、混ぜない。

NCBI owner: `core/greedy_align.c`の`s_FindFirstMismatch`（約313–375行）。0〜block幅前後の全長、各laneの不一致・fence、左右方向、曖昧文字を比較する。短い一致runでsetup費が勝つ場合は採用しない。

## 5. `task blastn`実装計画

### N0: DP経路と費用を固定する

対象: [blastn run.rs](../LOSAT/src/algorithm/blastn/blast_engine/run.rs)、[extension.rs](../LOSAT/src/algorithm/blastn/extension.rs)、[alignment/gapped.rs](../LOSAT/src/algorithm/blastn/alignment/gapped.rs)。

NCBI `api/blast_nucl_options.cpp:176–183`がscore-only DPとtraceback DPを選ぶ。preliminaryとfinal tracebackは必要な別段階として維持する。prepare呼出数、batch数、allocation数、先行DP数と実際に採用したDP数、cell数、順序付き包含判定時間を測る。

### N1: 成功した先行DP結果と開始位置を一緒に保持

現行`run.rs:9088`で開始位置を準備し、`9242`で再計算している。

1. 成功した先行結果に`prepare_traceback`のtupleも添える。新しいglobal cacheやHSP全体のコピーは作らない。
2. 元の順番で包含判定を行い、生き残ったHSPの先行結果があれば、そのtupleとDP結果を一緒に取り込む。
3. 先行結果がないHSPは、従来どおりその時点でprepare＋DPを行う。batch開始時に包含され、後のendpoint置換で必要になるHSPの経路を残す。prepare失敗と未計算を無理に同一のcache状態へまとめない。
4. identity検査、tree挿入、座標補正、sort/purgeの時点は維持する。

NCBI owner: `core/blast_traceback.c:436–475`、`core/blast_gapalign.c:3323`の開始位置探索、`4163`の`AdjustSubjectRange`。

検証: seedが両方0、開始位置調整、subject範囲shift、負strand、prepare失敗、包含後のendpoint置換、batch境界、先行結果あり/なし。完了条件は、成功した先行結果の取り込み時のprepare重複がなく、rawと順序が一致すること。

### N2: バッチ用配列を再利用

`jobs`、`speculative_results`、各scratch slotの中間結果をsubject処理内で再利用し、`clear`/`drain`で必要な範囲だけを更新する。結果を元のprelim indexへ戻す所有者は一つに保つ。batch16、仕事の配分、包含判定時点は変えない。

NCBI owner: `core/blast_traceback.c:403–405,509–513,583–612`。空batch、16境界、thread数より少ないjobs、繰り返す大小batch、元indexへの復元を確認する。allocation数の減少と全体時間を測り、保持capacityの増大も確認する。

### N3: 4塩基変換をquery context内で再利用

`extend_hit_ungapped_approx_ncbi`の既存式から、各開始位置の1-byte変換値を作る。最小案は各contextに約query長byteの表を持つ方式。mask適用後の同じBLASTNA値から構築し、4通りの開始位相をすべて扱う。既存の同等表現があれば優先して使う。

NCBI owner: `core/na_ungapped.c:262–349`、特に`292–304,322–334`の4塩基式。曖昧文字を通常のACGT用2-bit変換に置き換えない。既存のscore table、X-drop、exact再計算への切替を維持する。word size <11の経路は独立に保つ。

検証: BLASTNA各値からなる4文字の組合せ、位置位相0/1/2/3、短長、masked/unmasked、context境界、左右延伸、exact切替直前/一致/直後。前計算を含む時間と追加メモリで採否を決め、少数extension入力の回帰を確認する。

N1〜N3後もDPが支配する場合は、生成Wasmで残るチェック・loadの位置を特定して次の局所案を作る。過去の固定band案をそのまま再採用したり、score-only段階を省略したりしない。

## 6. TBLASTX実装計画

### X0: リンク本体の費用を分ける

対象: [sum_stats_linking/linking.rs](../LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs)。preliminaryと再評価後の両linkingを別計測し、group別にhelper訪問数、再DP回数、木のupdate/rebuild回数・伝播深さを数える。更新前後で値が同じ親ノード数も数える。

既存診断ではソートよりリンク全体の費用が大きかったが、その値を現在の速度として使わない。両linkingのNCBI ownerは`core/blast_engine.c:871,1517`。片方を削除する案は対象外。

### X1: 最大値管理木の同じ更新の伝播を打ち切る

現行`DualMaximumTree::update`は、更新後の親が変わらなくてもrootまで再計算する。

1. leafを更新し、各親で従来と同じ`DualMaximum::reduce`を実行する。
2. 新旧の親の**両channelのsumとindex全体**が等しければ、それより上への伝播を止める。同じwinner indexだけで判定しない。
3. 子の値と、必要な親の値は更新済みの状態にする。後の兄弟更新もその状態を参照する。path検証、rebuildの時点、cutoffの加算、HSP除去の順番は維持する。
4. 葉・親・rootの意味と、後勝ちの同点順位を保つ。最悪計算量は現在と同じで、同じreduceを繰り返す部分だけを減らす。

NCBI owner: `core/link_hsps.c:605–650`の最大値選択とpath検査、`906–908`のcutoff加算、`955–980`のchain取り出し。木は既存Rust実装の管理方法であり、NCBIに木があるとは主張しない。

検証: NCBI順の線形選択と木を独立に比較し、updateごとにsum/indexを確認する。空/1/2/非2冪要素数、同点、片channelだけ変化、非winner更新、winner除去、cutoff equality、両channelの勝者が異なる場合、連続updateとrebuildを含める。

完了条件: 実行したreduce数が減り、raw・chainメンバー・E-value bits・replay順が一致すること。X0で木の費用が小さければ見送り、木の局所改善をTBLASTX全体の大幅短縮と扱わない。

### X2: リンク判定のメモリアクセスを改善（条件付き）

X0/X1後のprofileが`LhHelper`と`HspLink`へのアクセスを指す場合に限定する。最初の実験案は、large-gap判定で毎回使うsum/next_largerと、条件成立時に必要な座標・HSP参照の読み方を分ける局所変更。必要なら小さな連続配列への配置を別候補として評価し、helper構築・更新費用も含める。

NCBI owner: `core/link_hsps.c:659–768,827–895`。訪問helper indexの列、next_largerのjump、strict/non-strict比較、採用predecessorを同じにする。xsumは`(h_xsum + score * lambda) - logK`の評価順を維持し、`score * lambda - logK`の先取りで丸めを変えない。

既存のdead-field除去やループ抽出と異なる変更点、生成codeで減ったload/checkを説明できなければ実装を拡大しない。Native回帰も計測する。

## 7. BLASTP実装計画

### P0: 二つのKappa経路を分ける

対象: [blastp/blast_engine.rs](../LOSAT/src/algorithm/blastp/blast_engine.rs)、[kappa.rs](../LOSAT/src/algorithm/blastp/kappa.rs)、[gapalign.rs](../LOSAT/src/algorithm/blastp/gapalign.rs)、`core/composition_adjustment/`。

- 多query: queryごとのredoでscratchを再利用している。行列の不変部分、preliminary初期化、DPが候補。
- 単一query＋複数match: matchごとにquery準備とscratchを生成し、全matchを先行redoしてから順序付き早期終了判定をする。別の計測fixtureを用意する。

行列生成数、queryコピー/組成集計/hash生成数、diagonal消去bytes、scratch確保、先行redo数と早期終了で不使用になったredo数・DP cellsを計測する。全query入力に単一queryの費用を帰属させない。

### P1: 初期行列の不変部分を一度構築

`build_redo_align_params(q_idx)`内の`build_matrix_info`は、同一検索ではmatrixとscaled ideal lambdaが共通である。

1. 現在と同じ引数・演算で初期行列を最初に必要となる境界で一度構築する。no-hitやunsupported入力のerror発生条件を変えないよう、無条件に関数先頭へ移動しない。
2. 最小変更では完成済みの値をquery用paramsへ複製し、対数・丸めを伴う再生成をなくす。参照共有への型変更は、残るcopyが実測上問題となる場合だけ次の候補にする。
3. query長、mutable callback context、pair別の組成補正行列は従来の所有者に残す。可変workspaceをスレッド間で共有しない。

NCBI owner: `core/blast_kappa.c:2228–2236`、`composition_adjustment/composition_adjustment.c`の`Blast_Int4MatrixFromFreq`。初期行列全要素、λ bits、組成補正設定、no-hit、unsupported設定のerrorを検証する。多queryでも実測し、生成call数をquery数依存から検索あたり必要回数へ減らす。

### P2: 単一query経路の準備とscratchを再利用

不変の`build_query_workspace`結果を一度作り、queryコピー、組成集計、word hashを各matchで作り直さない。DP/composition/range用scratchは検索内の独立した作業単位に持たせ、NCBIの各match開始時に必要な状態だけをresetする。

NCBI owner: `core/blast_kappa.c:2309`以降の`s_GetQueryInfo`、`3329–3334,3493–3503`の作業領域。Rayonの`map_init/for_each_init`回数はthread数とは限らないため、実際の生成回数とretained memoryを測る。永久worker poolやcross-search cacheは導入しない。

検証: 同じquery＋異なるsubject、大小matchの連続、fresh/reused、retry、組成補正、no-hit、エラー後の次の検索。多queryの主fixtureは適用外対照として保護する。

### P3: subjectごとのdiagonal全消去を削減（状態証明が前提）

`prepare_independent_subject`の全`fill`費用を測った後、NCBIのoffset方式でscratchを再利用できるか調べる。全消去を単に削除する実装はしない。

NCBI owner: `core/blast_extend.c:162–184`の`Blast_ExtendWordExit`。scratchが受け取るsubject順、前回のoffset区間、rollover、last_hit/flag全状態を調べる。仕事の順番を仮定せず、境界の再初期化が必要な場合を明示する。検索候補列がfresh-stateと同じになる証明ができなければP3は保留する。

検証: subject順の変化、scratch間の分配、短長交互、極短subject、window境界、offset rollover、同一diagonalへの連続hit。clear量だけでなく、新しい管理branchの費用まで含めて採否を決める。

### P4: 単一queryの先行redoを有界化（条件付き）

P0で不使用redoのDP費用が大きい場合に限る。元のmatch順に小さな範囲を処理し、結果を同じ順にheapへ反映してから次範囲を判定する。新しいheap状態でNCBIの既存early termination条件が成立する候補のredoを省けるか検証する。

NCBI owner: `core/blast_kappa.c:3525–3539`。新しいpruning条件、近似score、感度変更は加えない。結果が必要になった時点のerror伝播も維持する。少数のbatch候補を事前に固定し、多数を測って最速だけを選ばない。

検証: heapが空/満杯、cutoff直前/同値/直後、同score/E-value、後続matchでheap更新、早期終了で不要になるerror、n1/n2/n4/n8。省いたredo/cell数とwall改善の双方を示す。

### P5: DP kernelは新しいprofileを根拠に選ぶ

P1〜P4を評価してもDPが支配するなら、既存の行単位書き込み最適化後のWasmをprofileし、cellあたりに残るscore lookup、subject境界処理、band更新を特定する。一つの不要演算を指定した局所差分として実装する。NCBI `core/blast_gapalign.c`のALIGN_EXを参照し、cell集合、X-drop、tie-break、fence、edit scriptを維持する。一般的な「DPをSIMD化」「調整行列をcache」といった未証明の案を先に採用しない。

## 8. 実装単位と最終検証

各M/N/X/P候補は別差分・別判定とし、担当ファイルとNCBI ownerを記録する。共有`blastn`ファイルではmegablast/task blastnの両方を確認し、共有protein kernelではblastp/tblastxの該当経路を確認する。

候補ごとの順序:

1. sourceと実測で仮説を固定する。
2. 局所実装と必要なNCBI境界unitを完成させる。
3. focused raw比較を行う。差があればその相違の修正を先に完了する。
4. 対象program/taskの比較と、影響するthread/targetの検証を行う。
5. 診断なしA/Bで時間・memoryを評価し、ACCEPTED / REJECTED / INCONCLUSIVE / SKIPPEDを記録する。複数候補の効果を足し算しない。

統合版では採用候補を合わせて再測定し、Rustの`cargo test`、`cargo clippy`、`cargo fmt --check`と対象比較を実行する。既存の無関係なformat差分を一括修正しない。

既存harnessを優先する:

- `LOSAT/tests/compare_blastn_parity.py`と各programのcomparison scripts。
- `LOSAT/tests/benchmark_wasm_threading.py`（baseline/candidateとrunnerパスを明示）。
- `LOSAT/tests/benchmark_wasi_reuse.js`、`check_wasm_threading.py`、`check_wasi_reactor.js`。
- `wasm_performance.py`はmanifestの現在の契約を確認して使う。同scriptのwarm経路はserial限定なので、threaded reuseの証拠に代用しない。

初回の大きな全matrixはS0、最終の全matrixは統合後に実施する。途中は変更に対応するfocused gateを使い、不採用候補に無関係な全suiteを繰り返さない。

最終確認範囲はNative/serial Wasm/threaded Wasm、n1/n2/n4/n8のraw、outfmt 0/6/7の影響範囲、同一reactorでの反復とthread数変更。ブラウザでは通常設定の実consumerでcold/反復、serial/threaded、cancel→次回成功、worker終了、memory budgetを確認する。共有runtimeを変更する必要が出た場合は独立した候補に分け、既存scheduler/API ownerを維持する。

成果物は新しいrun directoryにbaseline/candidateのhash、argv、raw出力、差分、全sample、診断counter、適用外・失敗・未実行の一覧、採用理由を保存する。[既存evidence template](../.agents/skills/verify-ncbi-parity-and-speed/references/evidence.md)を利用する。最後に4対象ごとに「削れた処理」「全体時間」「Nativeとの差」「残る制約」を報告する。

## 9. この計画作成時の検証

- 現行LOSATソース、NCBIローカルソース、既存の採用/不採用記録、比較harnessの入口を読解した。
- 変更は本計画書のみ。検索コード、fixture、既存証拠は変更していない。
- 新規build、parity試験、速度測定、実装候補の独立監査は未実施。上記はそれらを実行するための提案であり、効果や合格を宣言するものではない。
