# LOSAT Splign 実装計画（ペアワイズ・mRNA・エキソン表）

作成日: 2026-10-04。状態: **計画。実装・認証は未実施**。

NCBI Splign は、cDNA（mRNA）をゲノム配列にスプライス整列し、エキソン・イントロン
境界を求めるツールである。本書は Splign を LOSAT に純 Rust で移植する計画を定める。
最初の範囲は、ペアワイズモード（`-query`/`-subj`）、mRNA モード、エキソン表出力に
限る。この範囲では NCBI の出力をバイト単位で再現する。

Splign は NCBI BLAST+ ではなく NCBI C++ Toolkit の一部である。そのため、BLAST+ だけを
挙動の根拠とする現行の [AGENTS.md](../AGENTS.md) と
[PD-PURE-RUST-RUNTIME-AUTHORITY](product_decisions/PD-PURE-RUST-RUNTIME-AUTHORITY.md)
に、Splign に限った権威の追加が必要になる（第 3 節）。

## 1. 目的と完了の定義

公開コマンドは `LOSAT splign -query mrna.fa -subj genome.fa [options]` とする。
初期範囲の完了条件は次のすべてである。

1. **出力**: 第 4 節の範囲で、標準出力のエキソン表と `-log` ファイル（既定
   `splign.log`）が固定オラクルと生バイトで一致し、終了コードも一致する。
2. **経路**: NCBI のペアワイズ経路を、同じ順序・同じ入力状態で Rust に移植する。
   経路は megablast によるヒット取得、区画（compartment）探索、スプライス DP、
   後処理、向きの判定、整形からなる。外部で計算したヒットを受け取る入口は
   テスト専用に限る。
3. **範囲外の拒否**: 第 4.3 節の項目は黙って無視せず、明示的なエラーで終了する。
4. **独立性**: NCBI の実行ファイル・ライブラリ・FFI・subprocess を、LOSAT の
   実行時・ビルド・フォールバックに使わない。オラクルとハーネスは比較専用とする。
5. **参照コメント**: 移植した各 Rust 関数の直上に、固定コミット上の NCBI ファイル・
   行番号と該当断片を記す（AGENTS.md 第 4 項）。
6. **既存契約の維持**: BLASTN、BLASTP、TBLASTX、TBLASTN、進行中の BLASTX の出力・
   性能の契約を変えない。blastn エンジンに手を入れる段階 E では、既存の BLASTN
   gate を通す。
7. **決定性**: 同じ入力を繰り返し実行しても出力が同じである。

## 2. 決定事項

| ID | 決定 | 決定者 | 参照 |
|---|---|---|---|
| U-1 | Splign に限り、NCBI C++ Toolkit の固定コミットを挙動の根拠として認める | ユーザー（2026-10-04） | 3.1 |
| U-2 | 最初の範囲は「ペアワイズモード・mRNA モード・エキソン表出力」 | ユーザー（2026-10-04） | 4 |
| U-3 | バッチモードの配列取得（LDS2、BLAST DB、GenBank）は現時点では予定しない | ユーザー（2026-10-04） | 4.3 |
| U-4 | 計画書を `docs/` に保存する | ユーザー（2026-10-04） | 本書 |
| D-1 | 固定コミットは `release/29.3.0`（`b174d53e550125070f4e0afeb3726aec87e57e0c`） | 推奨案を採用 | 3.2 |
| D-2 | 合否判定のオラクルは、固定コミットから自前ビルドした splign とする。NCBI FTP の配布バイナリは診断専用 | 推奨案を採用 | 3.3 |
| D-3 | 段階ごとの比較専用ハーネス H1〜H4 を、固定コミットのライブラリにリンクして作る | 推奨案を採用 | 3.4 |
| D-4 | 初期範囲で受け付けるオプションは 4.2 の表の通りとし、それ以外は明示的に拒否する | 推奨案を採用 | 4.2, 4.3 |
| D-5 | megablast のヒットは、LOSAT blastn エンジンにメモリ内 API を追加して得る。`CBl2Seq` と同じ入力状態であることを証明する | 推奨案を採用 | 段階 E |
| D-6 | 対象環境は native Linux x86_64 の単一スレッドとする。Wasm と他の OS は初期範囲外。性能目標は設けず、計測だけを行う | 推奨案を採用 | 9 |
| D-7 | 実装ブランチは、着手時に `origin/main` から `feature/splign-pairwise-mrna` を作る。v0.2.0 のリリース作業とは独立に進める | 推奨案を採用 | 10 |
| D-8 | 引数エラーと stderr の診断文は、成功・失敗の別と終了コードの一致を求め、文言の一致は求めない | 推奨案をユーザーが承認（2026-10-04） | 4.4 |
| D-9 | Rust の配置は第 6 節の通りとする | 推奨案を採用 | 6 |
| D-10 | NCBI 側の欠陥（クラッシュ、未定義動作、決定的なバグ）は、既存の NCBI 欠陥方針に従って集め、まとめてユーザーに判断を仰ぐ | 既存方針 | 12 |

## 3. 権威と固定コミット

### 3.1 製品判断（段階 A で記録する）

段階 A で `docs/product_decisions/PD-SPLIGN-TOOLKIT-AUTHORITY.md` を作り、AGENTS.md に
同じ範囲の限定条項を追記する。記録する内容は次の通り。

- **対象**: `LOSAT splign` サブコマンドと、そのために移植するコード（第 6 節の配置）。
- **権威**: `ncbi/ncbi-cxx-toolkit-public` のコミット
  `b174d53e550125070f4e0afeb3726aec87e57e0c`（タグ `release/29.3.0`）の C++ ソース。
  Splign 内部の megablast については、従来どおり BLAST+ 2.17.0 のソースが権威である。
  固定コミットの BLAST core/api が BLAST+ 2.17.0 と同一であることを証拠として残す（3.2）。
- **他プログラムに波及させない**: BLASTN、BLASTP、TBLASTX、TBLASTN、BLASTX の権威は
  BLAST+ のまま。Toolkit の他のコード（別のアライナーなど）を根拠に機能を追加する
  ことは認めない。
- **オラクル**: 固定コミットからビルドした splign と、比較専用ハーネス。実行時や
  ビルドの依存にしない。
- **引用形式**: `ncbi-cxx-toolkit-public@b174d53e:src/algo/align/splign/splign.cpp:1162-1300`
  のように、コミット・パス・行番号を記す。
- **固定コミットの変更**: PD の新しい版として扱い、全 fixture の再凍結を伴う。

### 3.2 固定コミットの選定根拠（2026-10-04 に確認）

- BLAST+ 2.17.0 のソース配布には Splign 本体がない。
  `c++/src/algo/align/Makefile.in` と `c++/src/app/splign/Makefile.in` は空の
  スケルトンで、`include/algo/align` は存在しない。
- 固定コミットの BLAST core/api（`src/algo/blast/{core,api}` と
  `include/algo/blast/{core,api}`、計 298 ファイル）は、BLAST+ 2.17.0 の配布ソース
  （`/home/kawato/tools/ncbi-blast-source/ncbi-blast-2.17.0+-src/c++`）と、`$Id$`
  行の展開の差を除いて一致する。したがって、オラクル内部の megablast は LOSAT の
  既存権威と同じ実装である。
- 他のタグとの比較:
  - `release/29.0.0`（2025-01-30）は BLAST 2.16.1。
  - `release/29.1.0`（2025-05-02）と `release/29.2.0` は BLAST 2.17.0 だが、
    `blast_gapalign.c`（`calloc` の引数順、ループの範囲チェックの順序）と
    `blast_node.cpp`（`n->Detach()`）が 2.17.0 の配布と異なる。
  - `release/29.3.0` 以降は一致する。
  - BLAST+ 2.17.0 の配布は `NCBI_PRODUCTION_VER 20250122`、
    `NCBI_DEVELOPMENT_VER 20241129` である。
- 2025-02-25 の Splign セグフォルト修正（GP-39840。`splign.cpp:2283`、`size_t` の
  アンダーフロー）を含む。
- 第 5 節の Splign 経路のファイルは、調査時点の main
  （`a27f87f96ef5df4bcb4c937d44e2e803463cd07d`、2026-10-03）と同一である。差が
  あるのは経路外の `align_compare.{cpp,hpp}`、`align_sort.cpp` と、テスト用の
  `splign/unit_test/Makefile.unit_test_splign.app` だけである。

### 3.3 オラクル

- **合否判定用**: 固定コミットから Release 構成でビルドした `splign`。ビルド手順、
  コンパイラ、依存ライブラリの版、バイナリの sha256、`-version-full` の出力を
  段階 A の証拠に残す。
- **診断専用**: NCBI FTP の `genomes/TOOLS/splign/linux-i64/splign.tar.gz`
  （Splign 2.1.0、Build-Date 2024-02-05、GCC 7.3、tarball の sha256
  `0ddaae214b8420560faff933078ce0f97d9d2519b1c20957d4fe1ea829948e90`）。GP-39840 の
  修正前で、内部の BLAST も古いので、合否判定には使わない。
- 版番号はどちらも 2.1.0 で区別できない。証拠には必ずビルド元のコミットか sha256 を記す。

### 3.4 比較専用ハーネス

固定コミットのライブラリにリンクする小さな C++ プログラムを `tools/splign_oracle/` に
置く。LOSAT のビルドには含めない。

| ID | 出力 | 用途 |
|---|---|---|
| H1 | `x_SetupBlastOptions` と同じ設定の `CBl2Seq` の結果。Seq-align ごとの順序、鎖、座標、Dense-seg、`num_ident`、`score`、`bit_score`、`e_value` と、`CBlastTabular` への変換・`FlipStrands` 後の値 | 段階 C・D の入力、段階 E の合否 |
| H2 | `CCompartmentAccessor` の結果（区画の範囲、鎖、状態、所属するヒット） | 段階 C の合否 |
| H3 | `CSplicedAligner16` の transcript とスコア（配列、パターン、端の自由度を指定） | 段階 B の合否 |
| H4 | `CSplign::Run` 後、整形前の区画と segment（範囲、identity、長さ、注釈、詳細、エキソンかギャップか、スコア） | 段階 D の診断 |

素の NW とバンド DP は、Toolkit の `nw_aligner` アプリ（固定コミットの
`src/app/nw_aligner`）でも照合する。このアプリはスプライス DP を扱えない。

## 4. 対応範囲

### 4.1 公開コマンド

`LOSAT splign -query <mRNA FASTA> -subj <genomic FASTA> [options]`。引数名は NCBI
splign と同じにする（単一ダッシュ、`-subj`、`-W`）。NCBI は
`CFastaReader(fAssumeNuc | fOneSeq)` で各 FASTA の最初の配列だけを読み、Bioseq の
最後の Seq-id を使う（`splign_app.cpp:573-595`）。LOSAT も同じにする。

### 4.2 受け付けるオプション

行番号は固定コミット上のもの。`cmdargs` は `src/algo/align/splign/splign_cmdargs.cpp`、
`app` は `src/app/splign/splign_app.cpp`、`splign` は `src/algo/align/splign/splign.cpp`。

| オプション | 既定値 | NCBI の定義 |
|---|---|---|
| `-query`, `-subj` | 必須 | app:138, 149 |
| `-W` | 28 | app:166 |
| `-mask_ranges` | なし | app:171 |
| `-direction` | `default`（mRNA では `auto`） | app:178, 1042 |
| `-log` | `splign.log` | app:190 |
| `-type` | `mrna`（`est` は拒否） | cmdargs:51 |
| スコアの上書き 8 種（`-match_score`、`-mismatch_score`、`-gap_opening_score`、`-gap_extension_score`、`-gt_ag_splice_score`、`-gc_ag_splice_score`、`-at_ac_splice_score`、`-non_consensus_splice_score`） | mRNA の既定値: 1000、-1044、-3070、-173、-4270、-5314、-6358、-7395 | cmdargs:58-112、splign:127-157 |
| `-compartment_penalty` | 0.55 | cmdargs:115、splign:596 |
| `-min_compartment_idty` | 0.70 | cmdargs:125、splign:534 |
| `-min_singleton_idty` | 未指定なら `min_compartment_idty` | cmdargs:132, 260-265 |
| `-min_singleton_idty_bps` | 9999999 | cmdargs:139 |
| `-min_exon_idty` | 0.75 | cmdargs:149、splign:481 |
| `-min_polya_ext_idty` | 1.0 | cmdargs:157、splign:490 |
| `-min_polya_len` | 1 | cmdargs:166、splign:499 |
| `-max_intron` | 1200000 | cmdargs:173、`compartment_finder.hpp:215` |
| `-min_hole_len` | 0 | cmdargs:181、splign:508 |
| `-trim_holes_to_codons` | false | cmdargs:190、splign:517 |
| `-max_space` | 4096（MB） | cmdargs:197、`nw_aligner.hpp:184` |
| `-max_part_exon_ident_drop` | 0.25 | cmdargs:206、splign:547 |

値の範囲の制約は `cmdargs:221-247` と `app:208` に従う。`ArgsToSplign`（cmdargs:252）
が値を `CSplign` に設定する順序もそのまま移植する。

### 4.3 明示的に拒否するもの

| 項目 | 理由 |
|---|---|
| `-hits`, `-comps` | バッチモード。U-2 の範囲外 |
| `-mklds`, `-ldsdir`, `-blastdb`, `-genbank`, `-query_id`, `-subj_id` | 配列取得の基盤が必要。U-3 により予定しない |
| `-disc` | dc-megablast が LOSAT で未対応 |
| `-type est` | 初期範囲外（第 13 節の後続候補） |
| `-asn`, `-aln` | 初期範囲外（後続候補） |
| `-test` | NCBI の開発用モード（`20_28_90_cut20`、`20_28_plus`）。初期範囲外 |

これらが指定されたら、処理を始める前に非 0 で終了し、未対応であることを明示する。

### 4.4 出力の契約

- **標準出力**: `CSplignFormatter::AsExonTable`。flags は
  `eTF_NoExonScores | eTF_UseFastaStyleIds`（app:997, 1143）。
- **`-log` ファイル**: `x_LogCompartmentStatus`（app:471）が書く区画状態の行。既定で
  作業ディレクトリに `splign.log` を作る点も NCBI と同じにする。
- **終了コード**: NCBI と一致させる。
- **stderr**: 成功・失敗の別と終了コードの一致を求める。NCBI の診断フレームワークが
  出す文言（引数解析の usage 表示、`CFastaReader` の警告など）の一致は求めない（D-8）。
  差は証拠に記録する。

## 5. NCBI の処理経路と移植対象

| 処理 | NCBI（固定コミット） | 主な関数と行 | Rust の配置 | 段階 |
|---|---|---|---|---|
| 引数・入力 | `src/app/splign/splign_app.cpp`、`splign_cmdargs.cpp` | `Init` 95、`Run` 597、`x_ReadFastaSetId` 573、`SetupArgDescriptions` 47、`ArgsToSplign` 252 | `src/algorithm/splign/args.rs`、`app.rs` | F |
| ヒット取得 | `splign_app.cpp` | `x_SetupBlastOptions` 518、`x_GetBl2SeqHits` 848 | blastn のメモリ内 API と `src/algorithm/splign/hits.rs` | E |
| ヒットの表現 | `src/algo/align/util/{align_shadow,blast_tabular}.cpp`、`include/algo/align/util/{align_shadow,blast_tabular,hit_comparator,hit_filter}.hpp` | `CBlastTabular(const CSeq_align&)` 51、`FlipStrands`、`CHitFilter::s_GetSpan` | `src/align/util/` | C |
| 区画探索 | `include/algo/align/util/compartment_finder.hpp` | `CCompartmentFinder::Run` 690、`CCompartmentAccessor::Run` 1404 | `src/align/util/compartment_finder.rs` | C |
| 区画ごとの整列 | `src/algo/align/splign/splign.cpp` | `Run` 1162、`x_SplitQualifyingHits` 757、`x_SetPattern` 812、`x_RunOnCompartment` 1430、`x_Run` 1935、`x_ProcessTermSegm` 3012、`x_GetGenomicExtent` 3081、`IsPolyA` 1380、`s_TestPolyA` 1392、`x_FinalizeAlignedCompartment` 1147、`x_IsInGap` 606、`x_LoadSequence` 618 | `src/algorithm/splign/splign.rs` | D |
| DP | `src/algo/align/nw/{nw_aligner,nw_band_aligner,nw_spliced_aligner,nw_spliced_aligner16}.cpp` | `CNWAligner::Run` 503、`x_Run` 533、`SetPattern` 845、`x_CheckMemoryLimit` 1098、`CBandAligner::x_Align` 110、`CSplicedAligner16::x_Align` 137 | `src/align/nw/` | B |
| segment 化と改善 | `src/algo/align/nw/nw_formatter.cpp`、`src/algo/align/splign/splign_exon_trim.cpp` | `CNWFormatter::MakeSegments` 864、`SSegment::ImproveFromLeft1` 222、`ImproveFromLeft` 357、`ImproveFromRight1` 496、`ImproveFromRight` 628、`Update` 798、`s_IsConsensusSplice` 836、`CSplignTrim` 46-614 | `src/align/nw/formatter.rs`、`src/algorithm/splign/exon_trim.rs` | B, D |
| 向きと ORF | `splign_app.cpp`、`splign.cpp`、`src/algo/sequence/orf.cpp` | `x_ProcessPair` 993（向きの分岐 1042-1141）、`GetCds` 1080、`COrf::FindOrfs` 336-380 | `src/algorithm/splign/app.rs`、`orf.rs` | D |
| 整形 | `src/algo/align/splign/splign_formatter.cpp`、`splign_app.cpp` | `AsExonTable` 100、`x_LogCompartmentStatus` 471 | `src/algorithm/splign/formatter.rs` | D |

経路の C++ 行数（ヘッダ込み、調査時点の main で計数。経路のファイルは固定コミットと同一）:

| 部分 | 行数 |
|---|---|
| splign 本体・後処理・整形・引数 | 6,300 |
| NW／スプライス DP（`nw_aligner`、`nw_band_aligner`、`nw_spliced_aligner{,16}`、`nw_formatter`） | 5,775 |
| ヒット処理・区画探索（`compartment_finder.hpp` 1,636、`hit_filter.hpp` 1,063 ほか） | 4,833 |
| アプリ（CLI） | 1,399 |
| ORF 探索 | 685 |
| 計 | 約 19,000 |

`-comps` と compart ツール用の `compart_matching`（2,557 行）は含まない。

段階 A の棚卸しでは、経路上の C++ ライブラリの挙動も移植対象として列挙する。

- `CFastaReader`: ID の解析、大文字・小文字、IUPAC の曖昧文字、`-` や N の連続の扱い、空配列。
- `CSeqVector` の IUPAC 表現。
- `CSeqMap` によるギャップ判定（`x_IsInGap` 606）。
- `CSeq_id::GetSeqIdString` と `AsFastaString` の書式。
- ostream による float・double の書式。例: log のスコアは既定の精度 6 で `761.788`。
- `std::stable_sort`。
- 整数型の幅と、例外から区画の状態への対応。

## 6. Rust の構成と再利用の方針

- `LOSAT/src/align/nw/`: `CNWAligner`、`CBandAligner`、`CSplicedAligner`、
  `CSplicedAligner16`、`CNWFormatter`（`SSegment` を含む）。
- `LOSAT/src/align/util/`: `CAlignShadow`、`CBlastTabular`、`CHitComparator`、
  `CHitFilter`（経路で使う関数）、`CCompartmentFinder`、`CCompartmentAccessor`。
- `LOSAT/src/algorithm/splign/`: `CSplign`、`CSplignTrim`、`splign_util`、
  `CSplignFormatter`（初期範囲は `AsExonTable` のみ）、引数、アプリ層（向きの判定、
  log 出力）、`COrf::FindOrfs`。
- `LOSAT/src/algorithm/blastn/`: メモリ内で HSP を返す API（段階 E）。既存の出力経路は
  変えない。
- `LOSAT/src/main.rs`、`cli.rs`: `splign` サブコマンド。
- FASTA は既存の rust-bio の読み込みを使ってもよいが、第 5 節の棚卸しで確定した
  `CFastaReader` の意味に合わせる層を挟む。
- 名前は NCBI の関数名に対応させる（AGENTS.md の慣例）。C++ のテンプレートは、
  実際に使われる型（`CBlastTabular`）に具体化してよい。
- NCBI の `stable_sort` には Rust の安定ソート（`sort_by`）を使い、比較関数の
  フィールドと向きを一致させる（PD-PURE-RUST-RUNTIME-AUTHORITY のソート方針）。

## 7. 数値と順序で守ること

- **DP のスコアは整数**（`TScore`）。イントロン長のペナルティ `ilen >> ibs`
  （`nw_spliced_aligner16.cpp:425, 466`）を含め、型の幅とオーバーフローの挙動を合わせる。
- **浮動小数点の型**: `CBlastTabular` の identity は `float(matches / aln_len)`
  （f32）。区画 DP の点数は double で、`compartment_finder.hpp:593-596` の定数
  （`kPenaltyPerIntronBase`、`kPenaltyPerIntronPos`）を含む。f32 と f64 を NCBI と
  同じ箇所で使う。
- **メモリ上限**: `x_CheckMemoryLimit`（`nw_aligner.cpp:1098`、
  `nw_band_aligner.cpp:446`）の double 計算と、例外から区画のエラー状態への対応を
  再現する。これは出力行に影響する。
- **ソート**: 経路上のソートはすべて `stable_sort` である
  （`compartment_finder.hpp:681, 705, 1093, 1133, 1343, 1419`、
  `hit_filter.hpp:161, 562, 932`、`splign.cpp:820`、`nw_aligner.cpp:594`）。入力の
  順序が同点の並びを決めるので、ヒットの順序は H1 と一致させる。
- **例外の分類**: `eNoHits`、`eNoAlignment`、`eMemoryLimit` などの種類が、区画の状態と
  log の行を決める（`splign.cpp:1283-1297`）。Rust では同じ分類を持つ `Result` で表す。
- 段階 A の棚卸しで、経路上の libm 呼び出し（あれば）と浮動小数点の出力書式を列挙する。

## 8. テストデータ

fixture、凍結した出力、生成スクリプトは `LOSAT/tests/splign/` に置き、manifest は
`LOSAT/tests/splign_parity_manifest.tsv` とする。

### 8.1 合成 fixture

生成スクリプトと乱数の種を固定して保存する。

| ID | 内容 | 主に通す経路 |
|---|---|---|
| S1 | プラス鎖、GT-AG のイントロン 4 個、poly-A 付き | 基本経路、poly-A の検出 |
| S2 | ゲノム上でマイナス鎖に置いた遺伝子 | 鎖の反転、座標 |
| S3 | GC-AG、AT-AC、非コンセンサスのイントロン | スプライス 4 種、`auto` の逆向き整列 |
| S4 | エキソン境界付近のミスマッチと挿入・欠失 | `ImproveFrom*`、`CSplignTrim` |
| S5 | 13 bp・26 bp 前後の短い末端エキソンと短い内部エキソン | `m_MinPatternHitLength`、末端 segment の処理 |
| S6 | ゲノム上の 2 コピー（パラログと、イントロンのない偽遺伝子） | 複数区画、`compartment_penalty` の境界 |
| S7 | エキソン近くに N の連続（`CFastaReader` が扱うなら `-` も）を含むゲノム | `x_IsInGap` |
| S8 | ヒットなし、または `min_compartment_idty` 未満 | 空・エラーの区画、log、終了コード |
| S9 | `max_intron` 付近のイントロン、`max_space` 超過 | 区画の分割、メモリ上限のエラー状態 |
| S10 | 逆鎖に長い ORF を持つ mRNA、先頭に poly-T を持つ mRNA | `auto` の向きの選択 |
| S11 | `-mask_ranges` の指定 | クエリのハードマスク |
| S12 | 4.2 の各オプションを既定値以外にする掃引 | 引数から `CSplign` への設定 |

### 8.2 実データ fixture

段階 A で accession.version とゲノム領域を確定し、取得元 URL と sha256 を記録する。
ゲノム側は遺伝子の前後を含む部分配列とし、オラクルと LOSAT の実行時間を抑える。

| ID | 選び方 |
|---|---|
| R1 | プラス鎖でエキソンの少ない短いヒト遺伝子の RefSeq mRNA と、対応するゲノム領域 |
| R2 | エキソンが 10 個以上の遺伝子 |
| R3 | マイナス鎖の遺伝子 |
| R4 | U12 型（AT-AC）イントロンを持つことが知られた遺伝子 |
| R5 | 近傍に偽遺伝子やパラログを持つ遺伝子（複数区画） |
| R6 | ヒト以外（イントロンの短い生物）の遺伝子 |

### 8.3 NCBI の単体テスト

`src/algo/align/splign/unit_test/unit_test_splign.cpp`（376 行）は、GenBank から gi で
配列を取り、`mrna_in.asn` と `est_in.asn` のヒットから `*_expected.asn` と比較する。
ネットワークでの取得が前提で、期待値は ASN（初期範囲外）なので、参考資料として扱う。

### 8.4 凍結するもの

各ケースについて、コマンド、オラクルの sha256、標準出力、log、終了コード、stderr、
H1・H2・H4 のダンプを凍結する。

## 9. 比較と合格条件

- **最終の合否**: 標準出力と log の生バイト比較と、終了コードの一致で判定する。
  正規化や並べ替えはしない。
- **部品段階の合否**: 段階 B〜D は、ハーネスの出力（H1〜H4、`nw_aligner`）との完全一致を
  条件とする。構造化した差分は診断にだけ使う。
- **決定性**: 同じ入力を 3 回実行し、同じ出力になる。
- **既存 gate**: 段階 E と F で、既存の BLASTN gate と純 Rust 境界チェック
  （`LOSAT/tests/check_pure_rust_runtime_boundary.py`）を通す。
- **性能**: オラクルと LOSAT の実行時間を参考として記録する。標準の手順（ウォーム
  アップ 1 回のあと 3 回計測し、中央値と範囲を報告）に従う。合否には使わない。

## 10. 段階別の実装計画

各段階の証拠は `docs/evidence/splign_stage_<段階>/` に置く。

### A. 権威・オラクル・棚卸し・fixture

1. PD の記録（3.1）と、AGENTS.md への限定条項の追記。
2. 固定コミットのソースを取得し、オラクルの splign をビルドする。ビルドの記録を残す（3.3）。
3. ハーネス H1〜H4 を作る（3.4）。
4. 経路を棚卸しする。初期範囲のオプションで到達する全関数と C++ ライブラリの挙動
   （第 5 節）を、「移植する／経路外／拒否する」に分類し、ファイル・行番号付きで
   `call_path_inventory.tsv` に記録する。後から監査で見つかった関数を一つずつ
   足すのではなく、ここで全量を洗い出す。
5. fixture（第 8 節）を作り、オラクルの出力を凍結する。
6. 診断として、LOSAT の `blastn -task megablast -dust no -word_size 28 -subject` を H1
   と比べ、段階 E の作業量を見積もる。

完了条件: PD と AGENTS.md の更新、オラクルのビルド記録、H1〜H4 の動作、棚卸し表、
凍結した fixture がそろっていること。

### B. NW の核

`CNWAligner`（`SetSequences`、`SetEndSpaceFree`、`SetPattern`、ガイド付きの `x_Run` の
単一スレッド経路、`x_Align`、`x_DoBackTrace`、`GetTranscript`、`ScoreFromTranscript`、
メモリ上限の判定）、`CBandAligner`、`CSplicedAligner`、`CSplicedAligner16`、
`CNWFormatter` と `SSegment` を移植する。

完了条件: 単体テストが `nw_aligner` の素の NW・バンド DP と H3 のスプライス DP に、
全ケースで一致すること。

### C. ヒットと区画

`CAlignShadow`、`CBlastTabular`（Seq-align 相当からの変換と `FlipStrands`）、
`CHitComparator`、`CHitFilter` のうち経路で使う関数、`CCompartmentFinder`、
`CCompartmentAccessor` を移植する。

完了条件: H1 のダンプを入力にして、全 fixture で H2 と一致すること。

### D. CSplign 本体・ORF・整形

第 5 節の「区画ごとの整列」「向きと ORF」「整形」の行と、`CSplignTrim`・
`splign_util` を移植する。テスト専用の入口で H1 のヒットと FASTA を受け取り、
エキソン表と log を出力する。

完了条件: 全 fixture で標準出力と log がオラクルと生バイト一致し、H4 とも一致すること。
オラクルは内部で `CBl2Seq` の同じヒットを使うので、この段階は megablast と切り離して
判定できる。

### E. megablast のヒット源

LOSAT blastn エンジンに、HSP を edit script ごと NCBI の結果順で返すメモリ内 API を
追加する。`CBlastOptionsFactory::Create(eMegablast)` の既定値に `SetWordSize(W)`、
`SetMaskAtHash(true)`、`SetDustFiltering(false)` を加えた状態（app:518-546）と、
LOSAT の入力状態が等しいことを棚卸しで示す。blastn アプリの `-task megablast` の
既定値との差も、この棚卸しで洗い出す。

完了条件: 全 fixture で H1 と完全に一致すること（順序、座標、Dense-seg、`num_ident`、
スコア、E 値）。既存の BLASTN gate に合格し、既存の出力が変わらないこと。
`blastn -dust no` の CLI ケースを BLASTN manifest に追加することも推奨する（必須ではない）。

### F. CLI の統合

`LOSAT splign` の引数解析（名前、型、制約、既定値、排他）、範囲外の拒否、FASTA の
意味、log ファイルの書き出し、終了コードを実装する。

完了条件: 全 fixture で生バイトが一致すること。拒否の契約、反復実行での決定性、
純 Rust 境界チェック、README・CHANGELOG・`--help` の更新もそろっていること。

### G. 独立監査と認証

`ncbi_parity_auditor` による独立した確認を受け、`docs/release/splign_certification.md`
に認証の範囲（第 4 節）と証拠を記録する。認証の範囲外を互換と主張しない。

## 11. 並行実行

- B と C は新しいモジュールで重ならないので、別の worktree で並行して進めてよい。
  どちらも段階 A の H1〜H3 が前提である。
- E は blastn エンジン（`run.rs` など）を変更する。他の blastn エンジン作業
  （エンジンの共有部分に触れる BLASTX の作業を含む）と同時に進めない。A の完了後なら、
  B〜D と並行してよい。
- D は B と C の後、F は D と E の後、G は F の後に行う。
- 調査、比較、分類などの機械的な作業は Sonnet のエージェントに任せる（既存の運用方針）。

## 12. 主なリスクと対処

| リスク | 対処 |
|---|---|
| LOSAT の megablast と `CBl2Seq` のヒットの差。既存の差（PD-BLASTN-HSP-CANONICALIZATION の同点のクラス、2026-09-13 計画時点で未解決の NZ 自己比較 1 行）と、未認証の `-dust no` が Splign に波及する | 段階 D を H1 の入力で切り離して判定する。段階 E で H1 と完全一致を確認する。差が残る場合は原因を分類し、まとめてユーザーに判断を仰ぐ |
| `CBl2Seq` と blastn アプリでオプションの既定値が違う | 段階 A と E の棚卸しで、オプション構造体の全項目を比べる |
| C++ ライブラリの挙動（`CFastaReader` の ID とギャップ、`CSeqMap`、IUPAC、ostream の数値書式） | 段階 A で棚卸しし、S7 などの fixture で固定する |
| メモリと時間。O(N×M) の DP は最大 4 GiB（`max_space` の既定値） | fixture の大きさを抑える。メモリ上限の判定とエラー状態は S9 で再現する |
| オラクルのビルドが難しい（SQLite3、LMDB などの依存） | 段階 A で手順を確定して記録する。必要なら依存の少ない構成でビルドし、出力が変わらないことを確認する |
| 移植中に NCBI 側の欠陥が見つかる（GP-39840 と同種の `size_t` のアンダーフローなど） | D-10。既存の NCBI 欠陥方針に沿って集め、まとめて判断を仰ぐ |
| `-trim_holes_to_codons` は mRNA の CDS 情報（`CBioseq_Handle`）を使う。FASTA 入力では情報がない | 段階 A で、FASTA 入力時の挙動（何もしないかどうか）を確認し、fixture で固定する |
| 範囲の膨張 | 第 13 節の項目は、初期範囲の認証が終わるまで着手しない |

## 13. 初期範囲の後の拡張候補（予定は未定）

- `-type est`（スコアの既定値と、`-direction` の既定が `both` になる点）。
- `-aln`（ペアワイズ表示）。
- `-asn`（Seq-align の ASN.1 テキスト出力と `CScoreBuilderBase` のスコア）。
- `-disc`（先に dc-megablast の実装が必要）。
- Wasm と他の native プラットフォーム。
- 性能の改善。

## 14. 調査の記録（2026-10-04）

- **配布バイナリでの確認**: NCBI FTP の配布バイナリを、合成データ（ランダム配列に
  エキソン 5 個、GT-AG イントロン 4 個、poly-A 25 塩基）で 2 回実行し、標準出力が一致
  した。既定の `auto` では、ORF の判定の結果として両方の向き（`-1` と `+1`）の
  モデルが出力された。
- **LOSAT の現状**:
  - CLI は blastn、blastp、tblastx、tblastn。blastn の `-task` は megablast と blastn
    だけで、dc-megablast は未対応。
  - BLASTN の認証（`docs/release/blastn_v0.1.0_certification.md`）は、DUST on、両鎖、
    native serial の 14 ケースが範囲である。
  - 大域整列（NW）、ASN.1 出力、ORF 探索、poly-A 検出はない。
  - blastn にはメモリ上で HSP を返す API がない（`run()` はファイルに書き出す）。
  - `src/` は約 137,300 行、`algorithm/blastn/` は 29,644 行。
- **調査資料の退避先**（リポジトリ外）: `/home/kawato/losat-splign-research-20261004/`
  - `NOTES.md`
  - 固定コミットの Splign 経路と BLAST core/api の抜粋（sha256 `36fc8cb7…`）
  - 現行 main の align 関連の抜粋（sha256 `b95f9a01…`）
  - 配布バイナリの tarball
  - 合成 fixture と出力
