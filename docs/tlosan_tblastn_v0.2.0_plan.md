# TLOSAN（TBLASTN）v0.2.0 総合計画

状態: 実装前の提案。TLOSAN は「LOSAT による NCBI TBLASTN の Rust 実装」を指す。
この文書は既存の v0.1.0 認証を拡張したという主張ではない。

実装作業ブランチ: **`feature/tlosan-tblastn-v0.2.0`**。2026-09-23 に
`origin/main` (`15ea61c101ca401a6742896f2cb356ed60406350`) から作成した。
依頼時点で `origin/dev` は存在せず、依頼者が `origin/main` を基点に指定した。
各実装セッションはこのブランチを使用し、開始時に現在のブランチと作業ツリーを
確認する。別の基点やブランチに切り替えて作業を進めない。

セッション別の実行指示は
[TBLASTN v0.2.0 セッションプロンプト](tlosan_tblastn_v0.2.0_sessions/README.md)
に分けた。各セッションは前段階の証拠と未解決事項を確認してから着手する。

## 1. 達成条件と範囲

公開コマンドは既存 CLI と同じ `LOSAT tblastn -query proteins.faa -subject genome.fna`。
タンパク質クエリを、探索時に六つの読み枠へ翻訳した塩基配列に照合する。
v0.2.0 の完了条件は次のすべてとする。

1. `-task tblastn` の NCBI ソースで定義されたローカル `-subject` 動作を Rust で実装する。
   `tblastn-fast`、PSI-TBLASTN、データベース検索など、別の経路はそれぞれソースを
   追跡して実装・認証するまで明示的に拒否する。既存の未対応オプションも黙って無視しない。
2. NCBI `gc.prt` に存在する **全 27 ID** を扱う:
   `1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,32,33`。
   NCBI 2.17.0+ の TBLASTN CLI はこのうち 32 以外の 26 ID だけを受理する。
   ユーザーの「すべて」という要件に従って LOSAT TBLASTN の CLI は 32 も受理する。
   この NCBI CLI との差を製品判断に明記し、32 は NCBI ソースと比較専用の
   NCBI API オラクルで検証する。無効な ID を標準コード 1 に置き換えない。
3. `-outfmt 0`、`6`、`7` をすべて実装する。デフォルト設定、対象として宣言した
   オプション・データで NCBI の出力バイト、行順、ヘッダー・フッターまで一致させる。
   非標準の subject genetic code による差だけは第 2 節の限定例外として分類する。
   `0` のアラインメント本文・統計部、`7` のヒットなしクエリも対象とする。
4. ネイティブの `-num_threads 1,2,4,8` で LOSAT の出力が同一であり、
   コード 1 は NCBI ローカルとバイト一致、非標準コードは第 2 節の限定契約に
   一致する。実際に複数の独立した仕事があるケースで
   並列経路が動いた証拠を残す。WASI serial と WASI threads を v0.2.0 の配布対象に
   含める場合も同じ入力でネイティブとのバイト一致を必須とする。
5. NCBI 実行ファイルはテスト・比較専用とし、LOSAT の実行時・ビルド・フォールバック
   経路から呼ばない。Rust の各移植箇所には、`AGENTS.md` に従い、直上に NCBI の
   ファイル・行番号と該当 C/C++ 断片を記す。

NCBI `tblastn_options.cpp:53-85` は word threshold 13、sum statistics、組成補正
`eCompositionMatrixAdjust`、subject genetic code の既定値を定義する。NCBI
`blast_engine.c:747-888` は被検索配列の翻訳、六つのフレームの探索・HSP 統合、
linking/E-value 計算のタイミングを定義する。既存の BLASTP/TBLASTX 関数を使う場合も、
この順序と入力状態が等しいことを証明する。

最初の基準プロファイルは `-task tblastn`、`-outfmt 0`、BLOSUM62、
word size 3、threshold 13、window size 40、gap 11/1、E-value 10、
SEG `12 2.2 2.5`、組成補正 2、`-db_gencode 1`、`-num_threads 1`。
NCBI `blast_prot_options.cpp:85-141`、`tblastn_options.cpp:53-85`、
`format_flags.cpp:37` と公式オプション表を根拠とし、着手時に実行ファイルの
`-help` と基準出力で確定する。`-max_intron_length 0` は
`blast_parameters.c:793-815` で既定の linking 長を選ぶ値であり、
「linking 無効」と解釈しない。

## 2. 非標準の被検索側遺伝暗号: 明示された製品要件

ユーザーは、TLOSATX と同様、TLOSAN のローカル `-subject` でも指定した
`-db_gencode` を尊重するよう明示した。全 27 ID について、翻訳、候補探索、
再評価、スコア、統計、座標、`outfmt 0/6/7` の表示まで選択した subject code
を一貫して使う。NCBI ローカル `-subject` が非標準コードを標準コード 1 に
戻す場合も、その動作を LOSAT に移植しない。

これは現行 `AGENTS.md` の TBLASTX にだけ認められた例外を TBLASTN へ
自動的に拡張することではない。実装の最初に、今回の明示的な要件に基づく
TBLASTN 限定の製品判断を記録し、`AGENTS.md` に同じ範囲と
コード 32 の CLI 受理を明記する。
例外は subject genetic code が非標準のローカル検索に限定し、NCBI の
実行順、候補集合の決め方、スコアリング、linking、フィルター、統計、
出力書式に別の相違を許さない。

NCBI CLI が受理する 26 ID は `-subject` と `-db` を、同じ入力・コード・
オプションで新規生成する。コード 32 は NCBI の `gc.prt:340-347` と
`FindGeneticCode(32)` を使う比較専用の NCBI C++ API ハーネスを用意し、
翻訳だけでなく TBLASTN の HSP・統計・表示まで照合する。NCBI 実行コードを
LOSAT の実行時・ビルド・配布物へ含めない。API ハーネスで出力契約を
検証できなければ、コード 32 を認証済みとしてリリースしない。

NCBI `-subject` はローカル検索固有の挙動とコード 1 の完全一致を確認する
オラクルとし、`-db` は指定コードが実際に適用された
翻訳・探索・報告の参照とする。`-db` と `-subject` はデータベース統計や
見出しが異なり得るため、生バイトを互いに直接比較しない。非標準コードの
LOSAT と NCBI ローカルとの差は、同じコード 1 対照と `-db` 結果を使い、
指定コードの翻訳に起因する項目だけと証明する。差を数えただけで合格にしない。

根拠として、NCBI の引数設定は `blast_args.cpp:1029-1056`、ローカル subject
の構築は `blast_args.cpp:2538-2557`、subject の genetic code 読み出しは
`blast_setup_cxx.cpp:800-810` にある。実行ファイルの現行挙動は最初の
比較で確定し、ソースからの推論と分けて証拠を残す。

## 3. 実装の依存順序

| 段階 | 作業 | 完了を示す成果物 |
| --- | --- | --- |
| A. 権威と基準値 | NCBI ソース・実行ファイルの版と SHA、LOSAT commit、CLI・出力オプションを固定。既存 v0.1.0 回帰とローカル遺伝暗号挙動を記録。TBLASTN のローカル遺伝暗号とコード 32 CLI 受理に限る製品判断・`AGENTS.md` を明示要件に合わせて更新。 | ソース対応表、再実行可能なコマンド、入力 SHA-256、差分分類、承認済みの境界記録 |
| B. CLI・遺伝暗号 | `tblastn` サブコマンド、NCBI 既定値、27 ID の検証を追加。TBLASTN の引数検証は 32 も受理し、既存 TBLASTX の CLI 契約は保持。`utils/genetic_code.rs` と `core/gencode_singleton.rs` の重複表と暗黙のコード 1 フォールバックを整理し、64 コドン × 全 ID を NCBI 表と照合。 | 全 ID の翻訳・不正 ID 拒否、既存 TBLASTX 回帰の維持 |
| C. 探索本体 | BLASTP のタンパク質クエリ符号化・lookup・伸長と TBLASTX の subject 翻訳を、NCBI TBLASTN の文脈・センチネル・フレーム順に接続。フレームごとの初期 HSP、統合、再評価、削除を NCBI のタイミングで実行。 | 各段階の先頭相違を説明できるトレースと小型入力の一致 |
| D. 統計・連結 | 有効長、翻訳 subject の長さ換算、sum statistics、既定の intron/linking、組成補正モード 2、Kappa redo、E-value・bit score・フィルター・順位付けを移植。BLASTP の `do_link_hsps=false` を TBLASTN に流用しない。 | コード 1 の完全一致、非標準コードの差が選択した翻訳コードだけに由来する証拠 |
| E. 表示 | タンパク質 query 座標と subject 塩基座標（正負六フレーム）、ID、フレーム、アラインメント文字列、0 の本文と統計、6 の各列、7 の各クエリ見出しとフッターを実装。 | 0/6/7 のコード 1 生バイト一致、非標準コードの限定差分、ヒットなしと複数クエリの一致 |
| F. 並列・配布 | クエリまたは subject の独立単位で並列化し、NCBI の入力順・比較器・削除順で決定的に集約。必要なら内部の独立計算のみ追加並列化。Native/WASI の該当モードを比較。 | 1/2/4/8 スレッドの同一 SHA-256、並列段階の稼働証拠、WASI 比較 |
| G. 認証 | 既存 3 プログラムの回帰、TBLASTN の全行、ビルド・fmt・clippy・適用可能な Wasm ゲート、性能測定、独立した読み取り専用監査。 | v0.2.0 の適用範囲・例外・証拠・未対応を記した認証記録 |

各段階の具体的な作業指示と引き継ぎ条件は、上記のセッションプロンプトに記す。
段階を完了できない場合は、満たしていない合格条件を明示して次段階へ進まない。

主な LOSAT 変更候補は `src/cli.rs`、`src/algorithm/mod.rs`、新設の
`src/algorithm/tblastn/`、`src/blastinput/`、`src/api/`、`src/utils/genetic_code.rs`、
`src/core/gencode_singleton.rs`、`src/report/`、`src/utils/threading.rs`。
BLASTP や TBLASTX の既存エンジンを大規模に分岐させる前に、共用する演算の
NCBI 引数と状態が等しい箇所だけを抽出する。`outfmt 0` は既存の簡易 pairwise
writer ではなく、NCBI の TBLASTN 表示経路を権威として照合する。

## 4. テストデータと比較表

テストマニフェストに case ID、入力、genetic code、全オプション、target、
threads、outfmt、NCBI 版、期待する契約区分を保持する。出力は毎回新しい一時
ディレクトリに生成し、古い `tests/*_out` を無条件に上書きしない。

| 層 | 入力・狙い | 必須確認 |
| --- | --- | --- |
| 全遺伝暗号の単位検証 | NCBI `gc.prt` の全 27 ID ごとに 64 コドン、曖昧塩基、終止を照合する小型データ | 27 表を全件比較。コード 32 は NCBI API ハーネスも使用。不正 ID は明示的エラーで、コード 1 に置換しない |
| 六フレーム境界 | 同じ既知タンパク質に対応する塩基領域を +1/+2/+3/−1/−2/−3 に配置。末端、部分コドン、`N`、終止、低複雑度、no-hit、同点 HSP も個別に作る | クエリ AA 座標、subject nt 座標・向き、フレーム、SEG、候補と表示 |
| 標準コードの実データ | `LOSAT/tests/fasta/AvCLPV.faa`（120配列）対 `AvCLPV.fasta`（約422 KB）。まず単一タンパク質、次に全件。さらに `PsCLPV.fasta` と組み合わせる | 自己・異種、複数クエリ、複数 HSP、並列集約 |
| 非標準コードの実データ | `AP027131.faa`（595配列）対 `AP027131.fasta`（約672 KB）、`AP027133.fasta`（約615 KB）。既存のコード 4 比較資産を利用 | コード 4 の subject 翻訳とローカル NCBI 境界。コード 1/4 の識別コドンを含む最小化例も保持 |
| コード 32 | `gc.prt` のコード 32 で標準コードと翻訳結果が異なるコドンを含む合成 subject。比較専用 NCBI C++ API ハーネスと対にする | TBLASTN CLI の受理、六フレーム、HSP・統計・0/6/7 表示をソース実装と比較 |
| 規模・性能 | 既存 `EDL933.fna` / `Sakai.fna` を subject 候補とし、対応するタンパク質クエリを既存の注釈または明示した取得元からチェックサム付きで用意 | 長配列のチャンク境界、ピークメモリ、現実的な候補数、並列稼働 |

合成ケースの正解出力は手書きしない。固定した NCBI 実行ファイルで生成する。
高速な基本マトリクスは `27 コード × 3 outfmt × threads {1,4}` の 162 条件とし、
NCBI CLI が受理する 26 コードは対応ケースとまずバイト比較する。
コード 32 は NCBI API ハーネスの対応ケースと比較する。非標準コードによる差は
第 2 節の限定例外に照らし、コード以外の相違を拒否する。コードごとに
変化するコドンを含むケースを置き、単なる引数受理試験にしない。
代表的なコード 1/4/11 と実データでは threads 2/8、複数 subject、
同点、no-hit を追加する。
並列化で順序が変わり得るケースは複数回実行し、全 SHA-256 の一致を確認する。

比較はまず生バイトで行う。差があれば、最初の不一致を query/subject ID、
フレーム、開始終了座標、raw score、bit score、E-value、identity、gaps、
HSP 順序、出力文字列へ分解する。差分のある全ケースを分類してから同じ
NCBI 呼び出し段階を一括修正し、候補集合を変える近似や隠れた fallback を入れない。

## 5. 並列・Wasm・性能の合格条件

- NCBI `-subject` はスレッド数を実質 1 に制限し得るため、NCBI の
  `-subject -num_threads 4` を「NCBI の 4 スレッド性能」の証拠にしない。
  コード 1 では同じローカル subject ケースを出力オラクルとし、
  LOSAT の 1/2/4/8 の生バイトが一致することを確認する。非標準コードでは
  各スレッドの LOSAT 出力同士のバイト一致と第 2 節の限定契約を確認する。
- 並列時もフレーム内・フレーム間 HSP 統合、linking、Kappa redo、hitlist の
  NCBI 順を維持する。1 クエリ・1 subject しかない場合に無理に処理を
  分割して出力を変えない。複数クエリ/subject で実際の並列仕事が走ることを
  別の計測ケースで示す。
- plain `wasm32-wasip1` は serial。実スレッドの検証対象は
  `wasm32-wasip1-threads` と `wasm-threads`。両者を混同しない。
- 正式な速度測定は同じ release build・入力・出力先・スレッド数を用い、
  untimed warmup 1 回と timed 3 回の中央値・全範囲を報告する。
  NCBI の速度測定では事前に `makeblastdb` した `-db` を使い、DB 作成時間・
  コマンド・版・入力 SHA-256 を別に記録する。速度向上でバイト一致を
  崩した経路は採用しない。

## 6. 完了判定と現時点のリスク

v0.2.0 TLOSAN の認証は、全 27 ID と 0/6/7 の宣言済みマトリクス、
1/2/4/8 スレッドの一致、既存 BLASTN/BLASTP/TBLASTX の回帰、適用可能な
Native/Wasm 出力、ソース参照コメント、性能の出力 SHA 記録が揃って初めて行う。
発売向けのパリティまたは性能主張には `ncbi_parity_auditor` の独立監査を通す。

最初の必須作業は、明示された TBLASTN のローカル subject genetic code
要件を製品判断と `AGENTS.md` に反映し、NCBI `-subject` / `-db` の現行差を
各コードで特徴付け、コード 32 の比較専用 API オラクルを成立させることである。
主な技術リスクは TBLASTN 既定の
uneven-gap linking と組成補正の組合せ、
`outfmt 0` の全ヘッダー・統計表示、並列時の同点 HSP 順序である。
これらの差を未解決のまま「TBLASTN 対応」と記載しない。

### NCBI の主要ソース

- `c++/src/algo/blast/blastinput/tblastn_args.cpp:45-132`
- `c++/src/algo/blast/blastinput/blast_args.cpp:997-1056,2538-2557`
- `c++/src/algo/blast/api/tblastn_options.cpp:46-85`
- `c++/src/algo/blast/api/blast_objmgr_tools.cpp:110-145`
- `c++/src/algo/blast/api/blast_setup_cxx.cpp:800-810`
- `c++/src/algo/blast/core/blast_engine.c:747-888`
- `c++/src/algo/blast/core/blast_parameters.c:781-817,1028-1032`
- `c++/src/algo/blast/core/link_hsps.c:266-318,1765-1810`
- `c++/src/algo/blast/core/blast_kappa.c:2423-2453`
- `c++/src/algo/blast/core/blast_hits.c:1084-1143`
- `c++/src/objtools/align_format/tabular.cpp:907-1095`
- `c++/src/objtools/align_format/showalign.cpp:310-335`
- `c++/src/objtools/align_format/format_flags.cpp:37`
- `c++/src/objects/seqfeat/gc.prt:105-357`

公式 CLI オプション表: https://www.ncbi.nlm.nih.gov/sites/books/NBK279684/table/appendices.T.tblastn_application_options/

ソース位置は着手時に凍結した NCBI checkout で再確認する。
