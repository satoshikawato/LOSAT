# LOSAT Magic-BLAST 実装計画（`-subject`・短いリード・SAM/tabular）

作成日: 2026-10-04。状態: **計画。実装・認証は未実施**。着手は BLASTX v0.2.0 と
Splign の後とする（U-7）。

NCBI Magic-BLAST は、RNA-seq や DNA-seq のリードをゲノムや転写産物に整列する
マッパーで、BLAST の core の上に作られている。リードから lookup table を作り、参照を
走査して種を見つけ、jumper と呼ぶ表駆動の伸長で整列し、HSP をつないでスプライス整列と
ペアを作る。本書は Magic-BLAST を LOSAT に純 Rust で移植する計画を定める。最初の範囲は、
参照を `-subject` の FASTA で与え、リードを FASTA・FASTQ・FASTC で読み、SAM と
tabular で出力する場合に限る。この範囲では NCBI の出力をバイト単位で再現する。

Magic-BLAST の検索部分（core、api、blastinput）は BLAST+ 2.17.0 のソースに含まれる。
一方、アプリ層（SAM・tabular の整形、スレッド、未整列リードの出力）は BLAST+ の
配布に含まれず、NCBI C++ Toolkit と Magic-BLAST のソース配布にだけある。そのため
現行の [AGENTS.md](../AGENTS.md) と
[PD-PURE-RUST-RUNTIME-AUTHORITY](product_decisions/PD-PURE-RUST-RUNTIME-AUTHORITY.md)
に、アプリ層に限った権威の追加が必要になる（第 3 節）。

## 1. 目的と完了の定義

公開コマンドは `LOSAT magicblast -query reads.fq -subject genome.fa [options]` とする。
初期範囲の完了条件は次のすべてである。

1. **出力**: 第 4 節の範囲で、本体の出力（標準出力または `-out`）と
   `-out_unaligned` のファイルが固定オラクルと生バイトで一致し、終了コードも一致する。
   コマンドライン文字列の 2 箇所だけは 4.4 の規則で比べる（D-3）。
2. **経路**: NCBI の経路を、同じ順序・同じ入力状態で Rust に移植する。経路は、リードの
   読み込みとバッチ分割、リードの品質フィルタ、ハッシュ lookup table と参照の単語数
   フィルタ、参照の走査、jumper 伸長、mapper の HSP writer（チェイン化、スプライス、
   ペアの選択）、最終処理、結果の組み立て、整形からなる。外部で計算した中間結果を
   受け取る入口はテスト専用に限る。
3. **範囲外の拒否**: 第 4.3 節の項目は黙って無視せず、明示的なエラーで終了する。
4. **独立性**: NCBI の実行ファイル・ライブラリ・FFI・subprocess を、LOSAT の
   実行時・ビルド・フォールバックに使わない。オラクルとハーネスは比較専用とする。
5. **参照コメント**: 移植した各 Rust 関数の直上に、NCBI のファイル・行番号と該当断片を
   記す（AGENTS.md 第 4 項）。
6. **既存契約の維持**: 新しいモジュールとして作り、BLASTN、BLASTP、TBLASTX、TBLASTN、
   BLASTX のエンジンは変更しない。共有部品（塩基の符号化、FASTA の読み込み、Seq-id の
   書式、スレッドプール）に手を入れる場合は、該当する既存 gate を通す。
7. **スレッド数と決定性**: 同じ入力なら、`-num_threads` の値や実行の回数によらず同じ
   出力になる（U-2、4.5）。
8. **大きなゲノム**: native Linux で、ヒトゲノム（GRCh38）全体を `-subject` に与えた
   実行が完走する（U-5、8.4）。

## 2. 決定事項

| ID | 決定 | 決定者 | 参照 |
|---|---|---|---|
| U-1 | Magic-BLAST のアプリ層に限り、NCBI C++ Toolkit の固定コミットを挙動の根拠として認める | ユーザー（2026-10-04） | 3.1 |
| U-2 | LOSAT は `-num_threads` によらず入力順で出力する。合否は NCBI の `-num_threads 1` と比べる | ユーザー（2026-10-04） | 4.5 |
| U-3 | バッチの切り方を NCBI と同じにする | ユーザー（2026-10-04） | 4.5 |
| U-4 | SRA 入力は対象外とし、指定されたら明示的なエラーにする | ユーザー（2026-10-04） | 4.3 |
| U-5 | Wasm で大きなゲノムを扱えなくてよい。native で扱えればよい | ユーザー（2026-10-04） | 1, 8.4 |
| U-6 | 既定値と挙動はソースに従う。ヘルプ文・文書・論文の記述とは合わせない | ユーザー（2026-10-04） | 4.2, 12 |
| U-7 | 着手の順は BLASTX v0.2.0、Splign、Magic-BLAST とする | ユーザー（2026-10-04） | 10 |
| U-8 | 最初の範囲は、参照を `-subject` で与える場合に限る。`-db` は後続とする | ユーザー（2026-10-04） | 4.3, 14 |
| U-9 | 計画書を `docs/` に Splign 計画書と同じ形式で保存する | ユーザー（2026-10-04） | 本書 |
| U-10 | N-1: worker thread の中の入力エラーは、出力のバイト列を NCBI と同じにしたうえで終了コード 1 にする（承認済みの例外） | ユーザー（2026-10-04、推奨案） | 12 |
| U-11 | N-2: 出力の書き込みの失敗は終了コード 6 にする（承認済みの例外） | ユーザー（2026-10-04、推奨案） | 12 |
| U-12 | N-3: 握りつぶされる例外は、段階 A で再現する入力を探し、なければ LOSAT では明示的なエラーにする | ユーザー（2026-10-04、推奨案） | 12 |
| U-13 | N-4〜N-7: 決定的な NCBI の挙動として、すべて再現する | ユーザー（2026-10-04、推奨案） | 12 |
| D-1 | PD は Splign と分けて `PD-MAGICBLAST-TOOLKIT-AUTHORITY` とする。固定コミットは Splign と同じ `release/29.3.0`（`b174d53e550125070f4e0afeb3726aec87e57e0c`）とする | 推奨案を採用（ユーザーが推奨に従うと回答） | 3.1, 3.2 |
| D-2 | NCBI の複数スレッドの実行とは、バッチ単位の並び順の差だけを許す。バッチの中身と、バッチ内の順序の差は許さない | 推奨案を採用 | 4.5 |
| D-3 | SAM の `@PG` 行の `CL:` と tabular の 2 行目は、NCBI と同じ規則で LOSAT 自身の argv から作る。比較ではこの値だけを両側で置き換える | 推奨案を採用 | 4.4 |
| D-4 | 合否判定のオラクルは、固定コミットから VDB 付きで自前ビルドした magicblast とする。NCBI FTP の 1.7.2 配布バイナリは診断専用 | 推奨案を採用 | 3.3 |
| D-5 | 段階ごとの比較専用ハーネス H1〜H5 を、固定コミットのライブラリにリンクして作る | 推奨案を採用 | 3.4 |
| D-6 | 初期範囲で受け付けるオプションは 4.2 の表の通りとし、それ以外は明示的に拒否する | 推奨案を採用 | 4.2, 4.3 |
| D-7 | 対象環境は native Linux x86_64 とし、単一スレッドと複数スレッドの両方を扱う。Wasm と他の OS は初期範囲外。性能目標は設けず、計測だけを行う | 推奨案を採用 | 9 |
| D-8 | NCBI が読む環境変数のうち、バッチを決める `BATCH_SIZE` と `BATCH_NUM_SEQS` だけを同じ意味で移植する。開発用のスイッチ（`MAPPER_*`、`BTOP_NO_SPLICE_SIGNALS`）は移植せず、オラクルはこれらを未設定で実行する | 推奨案を採用 | 4.5, 7 |
| D-9 | ソートは PD-PURE-RUST-RUNTIME-AUTHORITY の方針に従う（NCBI の比較関数と Rust の安定ソート、同点は入力順。glibc の `qsort` は模倣しない） | 既存方針 | 7 |
| D-10 | stderr の診断文は、成功・失敗の別と終了コードの一致を求め、文言の一致は求めない | Splign の D-8 と同じ方針を採用 | 4.4 |
| D-11 | NCBI 側の欠陥は、既存の NCBI 欠陥方針に従って第 12 節に集め、まとめてユーザーに判断を仰ぐ。計画時点の候補は U-10〜U-13 で決定済み。実装中に新しい候補が見つかったら、同じ方法で集めて判断を仰ぐ | 既存方針 | 12 |
| D-12 | 実装ブランチは、着手時に `origin/main` から `feature/magicblast-subject` を作る | 推奨案を採用 | 10 |
| D-13 | Rust の配置は第 6 節の通りとする | 推奨案を採用 | 6 |
| D-14 | 固定コミットのソースの取得とビルドは Splign と共有する（同じ checkout から splign と magicblast を作る） | 推奨案を採用 | 3.3 |

## 3. 権威と固定コミット

### 3.1 製品判断（段階 A で記録する）

段階 A で次の 2 つの PD を作り、AGENTS.md に同じ範囲の限定条項を追記する。

**`docs/product_decisions/PD-MAGICBLAST-TOOLKIT-AUTHORITY.md`**

- **対象**: `LOSAT magicblast` サブコマンドと、そのために移植するコード（第 6 節の配置）。
- **権威**:
  - core、api、blastinput は、従来どおり BLAST+ 2.17.0 のソースとする。
  - アプリ層 `src/app/magicblast/`（`magicblast_app.cpp`、`magicblast_util.{cpp,hpp}`、
    `magicblast_thread.{cpp,hpp}`）と、アプリ層だけが使う
    `include/algo/sequence/consensus_splice.hpp` は、`ncbi/ncbi-cxx-toolkit-public` の
    コミット `b174d53e550125070f4e0afeb3726aec87e57e0c`（タグ `release/29.3.0`）とする。
  - 固定コミットの core・api と blastinput は BLAST+ 2.17.0 と同一なので、2 つの権威は
    食い違わない（3.2）。
- **他プログラムに波及させない**: BLASTN、BLASTP、TBLASTX、TBLASTN、BLASTX の権威は
  BLAST+ のまま。アプリ層を根拠に、他のプログラムへ機能を足すことは認めない。
  `blast_sra_input`（SRA 入力）は権威に含めない（U-4）。
- **オラクル**: 固定コミットからビルドした magicblast と、比較専用ハーネス（3.3、3.4）。
  実行時やビルドの依存にしない。
- **引用形式**: アプリ層は
  `ncbi-cxx-toolkit-public@b174d53e:src/app/magicblast/magicblast_util.cpp:848-896`
  のように、コミット・パス・行番号を記す。core・api・blastinput は従来どおり
  BLAST+ 2.17.0 のパスと行番号を記す。
- **固定コミットの変更**: PD の新しい版として扱い、全 fixture の再凍結を伴う。
  Splign の PD と同時に変える必要はないが、同じコミットを使う間はビルドを共有する。

**`docs/product_decisions/PD-MAGICBLAST-BATCH-ORDER-AND-CMDLINE.md`**

- 4.5 のバッチと出力順の取り決め（U-2、U-3、D-2、D-8）。
- 4.4 のコマンドライン文字列の取り決め（D-3）。
- どちらも出力の比較方法を定めるもので、検索・スコア・整形の差は認めない。

### 3.2 固定コミットの選定根拠（2026-10-04 に確認）

- BLAST+ 2.17.0 のソース配布では、`src/app/magicblast/` と
  `src/algo/blast/blast_sra_input/` は空のスケルトン `Makefile.in` だけである。
  `scripts/projects/blast/project.lst` にも入っておらず、2.17.0 のビルドログでも
  何も作られていない。
- Magic-BLAST の最新リリースは 1.7.2（2023-04-19）である。以後のアプリ層の変更は、
  2026-07-22 に main で版番号を 1.8.0 に上げたもの（SB-4661）だけで、1.8.0 は未リリース。
- 固定コミットの `src/app/magicblast/` の 8 ファイルは、1.7.2 のソース配布と同一である。
  差は `Makefile.magicblast.app` のビルド情報（`WATCHERS` と `LIBS` の順）だけ。
  版の定数は 1.7.2（`magicblast_app.cpp:60-62`）。
- 固定コミットの blastinput（`src` と `include` の 46 ファイル）と core・api
  （298 ファイル）は、BLAST+ 2.17.0 と `$Id$` 行の展開の差を除いて一致する。
- 1.7.2 のソース配布に同梱のライブラリは古い（`NCBI_PRODUCTION_VER 20230304`）。
  2.17.0 との差のうち挙動に関わるものは次の通りで、1.7.2 は権威にしない。
  - `blast_hits.c`: `Blast_HSPListPurgeNullHSPs` 後の `best_evalue` の再計算。
    共通開始点・共通端点の purge での `subject.frame` の比較。
  - `blast_nalookup.c`: MB lookup の `longest_chain` の計算。
- 初期範囲の経路で Toolkit にしかないのは、アプリ層と `consensus_splice.hpp` である。
  後者はヘッダだけで、`NStr` だけを使う。
- 固定コミットは Splign 計画の D-1 と同じである。

### 3.3 オラクル

- **合否判定用**: 固定コミットから Release 構成でビルドした `magicblast`。
  - アプリは `blast_sra_input` を無条件にリンクし、VDB を必須とする
    （`Makefile.magicblast.app:5-6, 16, 18`、`CMakeLists.txt:4`）。VDB がないと、
    この target は黙ってビルド対象から外れる。VDB は SRA アクセッションの読み込みにだけ
    使われ、`-sra` を使わない比較では呼ばれない。
  - VDB は、ローカルにある bioconda の ncbi-vdb（3.3.0 または 3.4.1）を `--with-vdb` で
    指定するか、`--with-downloaded-vdb` で取得する。どちらにしたかと版を記録する。
  - 構成は、1.7.2 配布バイナリの `-version-full` が示す手順（`--with-static --with-mt
    --with-openmp --with-flat-makefile --with-vdb=... --with-static-vdb
    --with-projects=scripts/projects/magicblast/project.lst` など）を雛形にする。
  - ビルド手順、コンパイラ、依存ライブラリの版、バイナリの sha256、`-version-full` の
    出力を段階 A の証拠に残す。
  - 同じ checkout から Splign のオラクルも作る（D-14）。
- **診断専用**: NCBI FTP の `blast/executables/magicblast/1.7.2/`
  `ncbi-magicblast-1.7.2-x64-linux.tar.gz`（sha256
  `93301d1816fd87fe64bb48950be82e5bf02a45bea81d63fff78c11e9908f4604`、バイナリの sha256
  `a5206fb5b12f6264885a3f11d02edefbcc81b2a0d41b31ca5f05f68641426881`、2023-04-18 ビルド、
  GCC 7.3、revision 665961）。同梱ライブラリが古い（3.2）ので合否判定には使わない。
  段階 A で全 fixture について両者の出力を比べ、差を記録する。
- 版の表示はどちらも 1.7.2 で区別できない。証拠には必ずビルド元のコミットか sha256 を記す。
- アプリは各 worker thread の中で検索を単一スレッドで行う（`CMagicBlast` の既定の
  スレッド数は 1、`magicblast_thread.cpp:155, 169`）。経路上の OpenMP
  （`blast_nalookup.c:1860`）は `num_threads > 1` のときだけ動くので、使われない。
- オラクルは、`MAPPER_*`、`BTOP_NO_SPLICE_SIGNALS`、`BATCH_SIZE`、`BATCH_NUM_SEQS` を
  未設定（fixture が指定する場合を除く）、`LC_ALL=C` で実行し、argv と環境を記録する。

### 3.4 比較専用ハーネス

固定コミットのライブラリにリンクする小さな C++ プログラムを `tools/magicblast_oracle/`
に置く。LOSAT のビルドには含めない。公開 API で取り出せない内部状態は、ハーネス専用の
ビルドに限ってダンプ用のコードだけを足す。そのビルドの最終出力が、手を加えない
オラクルと全 fixture で一致することを確かめてから使う。

| ID | 出力 | 用途 |
|---|---|---|
| H1 | バッチごとのリード（ID、塩基、品質、ペアの印）と、リードの品質フィルタの結果（`FilterQueriesForMapping`） | 段階 B の合否 |
| H2 | ハッシュ lookup table の中身（16-mer ごとのリード上の位置）、参照の単語数フィルタで除いた語、poly-A の語の除去 | 段階 C の合否 |
| H3 | (a) `BlastNaExtendJumper` に渡る種の組（リード位置、参照位置）を到着順に。(b) subject の chunk ごとに writer へ渡る HSP list（座標、スコア、edit、スプライス信号、overhang）を到着順に | (a) 段階 C の合否、段階 D の入力。(b) 段階 D の合否、段階 E の入力 |
| H4 | `BlastHSPStreamMappingClose` 後の `BlastMappingResults`（リードごとのチェイン、HSP、ペア、count、adapter、poly-A） | 段階 E の合否、段階 F の入力 |
| H5 | `CMagicBlastResults`（Spliced-seg、"Mapper Info" の各値、フラグ、concordance、`SortAlignments` 後の順） | 段階 F の診断 |

NCBI の単体テスト `magicblast_unit_test.cpp` も固定コミットでビルドして実行し、
オラクルのビルドが正しいことの確認に使う（8.3）。

## 4. 対応範囲

### 4.1 公開コマンド

`LOSAT magicblast -query <reads> -subject <reference FASTA> [options]`。引数名は NCBI
magicblast と同じにする（単一ダッシュ）。

- **リード**: FASTA、FASTQ、FASTC。名前が `.gz` で終わるファイルは展開して読む。
  `-query` を省くと標準入力から読む（NCBI の既定値 `-`）。ペアは `-paired`
  （連続する 2 本）、`-query_mate`（2 つ目のファイル）、FASTC で与える。
- **参照**: `-subject` の FASTA（複数配列可、`.gz` 可）。`-subject_loc` で範囲を指定できる。
- **出力**: SAM（既定）と tabular。未整列リードは本体に混ぜるか、`-out_unaligned` に
  SAM・tabular・FASTA で書く。

### 4.2 受け付けるオプション

行番号は BLAST+ 2.17.0 のもの。`args` は `src/algo/blast/blastinput/blast_args.cpp`、
`margs` は `src/algo/blast/blastinput/magicblast_args.cpp`、`app` は固定コミットの
`src/app/magicblast/magicblast_app.cpp`。

| オプション | 既定値 | NCBI の定義 |
|---|---|---|
| `-query` | `-`（標準入力） | args:3427-3430。空入力の判定は app:339-343 |
| `-infmt` | `fasta`（`fastq`、`fastc` も受け付ける） | args:2020-2025, 2067-2096 |
| `-paired` | off | args:2026, 2070-2073 |
| `-query_mate` | なし（`-query` が必要） | args:2027-2032, 2098-2116 |
| `-parse_deflines` | true | args:2056-2058, 2546-2554 |
| `-lcase_masking` | off（subject の小文字だけをマスクする） | args:2013-2015, 2551-2554 |
| `-validate_seqs` | true（リードの品質フィルタ） | args:2016-2017, 2118-2120 |
| `-subject` | 必須（`-db` を使わないため） | args:2352-2357, 2536-2557 |
| `-subject_loc` | なし（すべての subject に同じ範囲を適用する） | args:2366-2378, 2548-2550 |
| `-word_size` | 18（12 以上） | margs:82-88 |
| `-gapopen` | 0 | margs:92-94 |
| `-gapextend` | 4 | margs:97-99 |
| `-penalty` | -4（0 以下） | margs:152-156 |
| `-max_intron_length` | 500000（0 以上） | margs:168-173 |
| `-score` | `0`（リード長に応じた閾値。`GetCutoffScore`、`jumper.c:4563`）。整数か `L,b,a` | args:1473-1483, 1518-1564 |
| `-max_edit_dist` | 制限なし | args:1484-1486, 1566-1568 |
| `-splice` | true | args:1487-1488 |
| `-reftype` | `genome`（`transcriptome` も受け付ける） | args:1489-1493, 1566-1577 |
| `-perc_identity` | 0（0〜100） | margs:103-107 |
| `-fr`, `-rf` | off（互いに排他。SAM 以外と組むと終了コード 1） | margs:110-132、app:363-371 |
| `-limit_lookup` | true | args:1494-1496 |
| `-max_db_word_count` | 30（2〜255） | args:1497-1503 |
| `-lookup_stride` | 0 | args:1504-1506 |
| `-out` | `-`（標準出力） | args:3441-3445 |
| `-outfmt` | `sam`（`tabular` も受け付ける） | args:3012-3019, 3052-3072 |
| `-out_unaligned` | なし（`-no_unaligned` と排他） | margs:59-65 |
| `-unaligned_fmt` | `-outfmt` と同じ（`sam`、`tabular`、`fasta`） | args:3021-3031, 3074-3086 |
| `-md_tag` | off | args:3033, 3128-3130 |
| `-no_query_id_trim` | off | args:3034-3036, 3096-3098 |
| `-no_unaligned` | off | args:3038, 3100-3102 |
| `-no_discordant` | off | args:3040-3041, 3104-3106 |
| `-tag` | なし（文字列 `\t` をタブに置き換える） | args:3043-3045, 3142-3144 |
| `-num_threads` | 1（CPU 数で頭打ちにする） | args:3148-3188、margs:204-224 |

値の範囲と排他の制約は `CArgDescriptions` の定義に従い、違反は終了コード 1 にする。
各 `Extract` 関数が値をオプション構造体に設定する順序もそのまま移植する。
`-h`、`-help`、`-version` は LOSAT の既存サブコマンドと同じ扱いとし、互換の対象外とする。

ソースとヘルプ文が食い違う点は、ソースに従う（U-6）。

- `-score` の既定の閾値: ヘルプ文（args:1476-1481）は「20 以下は長さ、30 以下は 20、
  50 以下は長さ − 10、それ以外は 40」と書く。実装の `GetCutoffScore`
  （`jumper.c:4563-4576`）は「20 以下は長さ、34 以下は 20、200 未満は
  `(Int4)(0.6 * 長さ)`、それ以外は 120」である。配布バイナリでも実装どおりに動いた。
- `-score L,b,a`: 係数は `(int)(c * 100)` で保存され、閾値は
  `(b100 + a100 * 長さ) / 100` の整数除算になる（`blast_options_local_priv.hpp:1465-1469`、
  `hspfilter_mapper.c:4758-4766`）。`a` が 0 なら `b` は無視される。
- `-reftype`: `-limit_lookup` が与えられていないときだけ、単語数フィルタの既定値
  （`genome` なら on）を決める（args:1566-1577）。`-limit_lookup` には既定値 true が
  あるので（args:1494-1496）、効果がない可能性がある。効果がある場合は、
  `transcriptome` がバッチの割り算（4.5）にも効く。段階 A で確かめ、挙動をそのまま
  移植する。

### 4.3 明示的に拒否するもの

| 項目 | 理由 |
|---|---|
| `-db` と DB 専用のオプション（`-gilist`、`-seqidlist`、`-negative_gilist`、`-negative_seqidlist`、`-taxids`、`-negative_taxids`、`-taxidlist`、`-negative_taxidlist`、`-db_soft_mask`、`-db_hard_mask`） | U-8。BLAST DB の読み込みが必要（第 14 節の後続候補） |
| `-sra`、`-sra_batch`、`-sra_cache` | U-4。VDB が必要で、純 Rust の方針と両立しない。恒久的に対象外 |
| `-infmt asn1`、`-infmt asn1b` | 初期範囲外（後続候補） |
| `-outfmt asn` | 初期範囲外（後続候補） |
| `-gzo` | 圧縮後のバイト列が zlib の実装と圧縮の細部に依存し、純 Rust の実装では生バイトの一致を示せない。初期範囲外（後続候補。比較方法の取り決めが先に要る） |

これらが指定されたら、処理を始める前に非 0 で終了し、未対応であることを明示する。

### 4.4 出力の契約

- **本体の出力**（`-out`、既定は標準出力）:
  - SAM: `PrintSAMHeader`（util:848-896）、`PrintSAM`（util:1008-1504）、
    `PrintSAMUnaligned`（util:1507-1603）。`util` は固定コミットの
    `src/app/magicblast/magicblast_util.cpp`。
  - tabular: `PrintTabularHeader`（util:396-437）、`PrintTabular`（util:441-683）、
    `PrintTabularUnaligned`（util:686-773）。
  - ヘッダは検索の前に一度だけ書く（app:383-411）。後で失敗しても、ヘッダは残る。
- **未整列リード**: 既定では、本体の同じ出力で、そのリードの整列の直後に書く。
  `-out_unaligned` を指定すると、`-unaligned_fmt` に応じたヘッダ付きで別ファイルに書く
  （FASTA にはヘッダがなく、defline 全体と 1 行の配列を書く。util:343-369）。
- **版の文字列**: tabular の 1 行目は `# MAGICBLAST 1.7.2` とする。LOSAT が既存の出力で
  NCBI の版文字列（`# BLASTN 2.17.0+` など）をそのまま出すのと同じ扱いである。SAM には
  版の文字列がない。
- **コマンドライン文字列（D-3）**: NCBI は argv の各要素（argv[0] を含む）の後に空白を
  1 つ付けて連結し（app:287-295）、SAM の `@PG` 行の `CL:` と tabular の 2 行目に書く。
  この値は NCBI 自身でも起動のしかた（`./magicblast`、PATH 経由、絶対パス）で変わる。
  LOSAT は同じ規則を自分の argv に適用する（例: `LOSAT magicblast -query r.fq -subject g.fa `）。
  比較では、SAM の `@PG` 行の `CL:` 以降（本体と `-out_unaligned` の両方）と、tabular の
  2 行目の `# ` 以降だけを、両側で同じ固定文字列に置き換える。他の行は生バイトで比べる。
- **形式に連動する内部設定**: `-outfmt` が tabular 以外のとき、NCBI のアプリは環境変数
  `MAPPER_NO_OVERLAPPED_HSP_MERGE=1` を設定する（args:3132-3140）。これにより core の
  チェイン化が変わる（`hspfilter_mapper.c:2327`）ので、SAM と tabular では整列そのものが
  変わり得る。LOSAT は、環境変数ではなく `-outfmt` から決まる内部設定として移植する。
  同様に、アプリが設定する `SEQ_ID_PREFER_ACCESSION_OVER_GI=1`（app:310-311）の効果
  （Seq-id の選び方）も移植する。
- **再現する NCBI の挙動**（決定的なので、そのまま再現する）:
  - flag 16 の記録でも QUAL を反転しない（util:1430-1433）。
  - リード ID の末尾 `.1`、`.2`、`/1`、`/2` を切るのは SAM でペアのときだけ
    （util:1205-1210、`magicblast_thread.cpp:67-68`）。
  - tabular の mate の参照列は `lcl|` 付きの `AsFastaString` で書く（util:654）。
  - tabular の compartment 列は `1:<番号>` で、番号はバッチごとに 0 から数え直す。
  - BTOP のイントロンは `^gt496ag^` の形で、マイナス鎖では逆順・相補にする。
  - MD タグの長い subject のギャップは `^` の後に `x` を並べる。
  - concordance の判定は参照の ID を比べない（`magicblast.cpp:686-696`）。
- **終了コード**: `BLAST_EXIT_*` と `CATCH_ALL`（`src/app/blast/blast_app_util.hpp:146-250`）の
  対応に従う。主なものは、成功 0、入力・引数の誤り 1、engine の誤り 3（空の subject
  など）、捕まらない例外 255。ただし、worker thread の中の入力エラーは 1、出力の
  書き込みの失敗は 6 とする（第 12 節の承認済みの例外、U-10、U-11）。
- **stderr**: 成功・失敗の別と終了コードの一致を求める。文言（usage の表示、
  `-fr` の誤りの文 `-oufmt` の綴りなど）の一致は求めない（D-10）。差は証拠に記録する。

### 4.5 スレッドとバッチ（U-2、U-3、D-2、D-8）

**NCBI の挙動**

- 各スレッドは、入力 mutex の下で次のバッチを順に取る（`magicblast_thread.cpp:112-120`）。
- バッチは、塩基数が 50,000,000（`GetQueryBatchSize`、`blast_input_aux.cpp:110-116`）か、
  配列数が 500,000（app:414）に達するまで詰める。いま読んでいるリードかペアは最後まで
  入れるので、ペアが分かれることはない（`blast_input.cpp:306-324`）。
- 各バッチは独立に、スレッド内では単一スレッドで検索する。結果は出力 mutex の下で、
  終わった順に書く（`magicblast_thread.cpp:228-239`）。
- したがって、スレッドが 2 つ以上でバッチが 2 つ以上あると、バッチの並び順が実行ごとに
  変わる。配布バイナリで、`BATCH_NUM_SEQS=100`、4 スレッドで 3 回実行し、並びは毎回
  違ったが、並べ替えた中身は同じだった。
- バッチの中身は、既定ではスレッド数に依存しない。ただし `-limit_lookup false` のときは
  `batch_size = MAX(batch_size / num_threads, 5000000)` になる（app:424-426）。
  ここでの `num_threads` は CPU 数で頭打ちにした後の値である（margs:204-224）。
- 1 バッチが 1,000 リードを超えると別の経路（`MapperWordHits`）を通る
  （`blast_engine.c:1054-1056`）。そのため、バッチの中身が変わると結果が変わり得る。

**LOSAT の規則**

1. バッチは NCBI と同じ式で切る。CPU 数での頭打ちと、`-limit_lookup false` のときの
   割り算も含める（U-3）。
2. バッチは並列に検索してよいが、出力は入力順に書く（U-2）。バッチの中身が同じなら、
   `-num_threads` の値によらず同じバイト列になる。
3. 環境変数 `BATCH_SIZE`（`blast_input_aux.cpp:86-91`）と `BATCH_NUM_SEQS`（app:417-420）を
   同じ意味で読む（D-8）。小さな fixture でも複数のバッチを作れる。

**比較の方法**

- **合否**: NCBI の `-num_threads 1` と比べる。バッチの中身がスレッド数に依存する場合
  （`-limit_lookup false`）は、N スレッドの LOSAT に対して、NCBI の 1 スレッドの実行に
  `BATCH_SIZE=max(D / N, 5000000)` を与え、同じバッチの中身にする。D は既定値
  50,000,000 か利用者の `BATCH_SIZE`、割り算は整数、N は頭打ち後の値とする。
- **NCBI の複数スレッドの実行**: 診断として比べる。許す差はバッチ単位の並び順だけ
  とする（D-2）。入力順からバッチの境界を決め、バッチごとのバイト列の集合が一致する
  ことを確かめる。

## 5. NCBI の処理経路と移植対象

`core` は `src/algo/blast/core`、`api` は `src/algo/blast/api`、`binput` は
`src/algo/blast/blastinput`（いずれも BLAST+ 2.17.0）。`app`、`thread`、`util` は固定
コミットの `src/app/magicblast/` のファイル。

| 処理 | NCBI | 主な関数と行 | Rust の配置 | 段階 |
|---|---|---|---|---|
| 引数・アプリ | `app/magicblast_app.cpp`、`magicblast_thread.cpp`、`binput/magicblast_args.cpp`、`blast_args.cpp` | `Run` app:298-467、`s_CreateInputSource` app:127-222、`s_InitializeSubject` app:226-260、`s_GetCmdlineArgs` app:287-295、`CMagicBlastThread::Main` thread:65-246、`CMagicBlastAppArgs` margs:227-288、`CMappingArgs` args:1466-1595、`CMapperQueryOptionsArgs` args:2008-2160、`CMapperFormattingArgs` args:2995-3144 | `src/algorithm/magicblast/{args,app}.rs` | G |
| リード入力・バッチ | `binput/blast_fasta_input.cpp`、`blast_input.cpp`、`blast_input_aux.cpp` | `CShortReadFastaInputSource` 506-1073、`CBlastInputOMF::GetNextSeqBatch` 294-324、`GetQueryBatchSize` 86-116 | `input.rs` | B |
| subject 入力 | `binput/blast_input_aux.cpp`、`objtools/readers/fasta.cpp` | `ReadSequencesToBlast` 222-247（`SetSubjectLocalIdMode`、`SetConvertGapsToNs`、小文字のマスク） | `subject.rs` | B |
| 設定 | `api/magicblast_options.cpp`、`core/blast_options.c`、`blast_parameters.c`、`blast_setup.c` | `SetRNAToGenomeDefaults` 70、既定値 133-217、`SReadQualityOptions` 177-200、検証 1714-1722、`s_JumperScoreBlkFill` 396-455、507-512、実効長 743-751 | `options.rs`、`setup.rs` | B |
| リードの品質フィルタ | `core/jumper.c`、`blast_filter.c` | `FilterQueriesForMapping` 4531、`s_FindDimerEntropy` 4482-4510、filter の分岐 575, 1023, 1162-1166, 1206, 1356 | `filter.rs` | B |
| lookup table | `core/blast_nalookup.c`、`lookup_wrap.c` | `BlastChooseNaLookupTable` 69-75、`BlastNaHashLookupTableNew` 2207-2278、`s_NaHashLookupScanSubjectForWordCounts` 1839-1960、`s_NaHashLookupRemovePolyAWords` 1965、lookup_wrap 120-152 | `lookup.rs` | C |
| 走査 | `core/blast_nascan.c`、`na_ungapped.c`、`blast_engine.c` | ハッシュの走査 2682-3007、`JumperNaWordFinder` 1930、`MapperWordHits` 1839-2130、`s_BlastSearchEngineCore` 702、補助の初期化 991-1062、chunk の重なり 461-471、528-536、580-584 | `scan.rs`、`engine.rs` | C |
| jumper 伸長 | `core/jumper.c`、`blast_gapalign.c` | `BlastNaExtendJumper` 3253-3820、`JumperGappedAlignmentCompressedWithTraceback` 2512、`JumperExtendRightCompressedWithTracebackOptimal` 1124、`JumperExtendLeftCompressedWithTracebackOptimal` 2110、`JumperGoodAlign` 2650、`s_ShiftGaps` 457、`JumperFindEdits` 2755、`JumperFindSpliceSignals` 2947、`s_SaveSubjectOverhangs` 2995、`GetCutoffScore` 4563、`JumperGapAlignNew` 313-336 | `jumper.rs` | D |
| HSP writer（チェイン化・スプライス・ペア） | `core/hspfilter_mapper.c`、`blast_hspstream.c`、`blast_hits.c` | `s_BlastHSPMapperSplicedPairedRun` 4635、`s_FindBestPath` 3707、`s_FindSpliceJunctions` 3518、`s_FindBestPairs` 4227、`HSPChainListInsert` 238、`BlastHSPMapperParamsNew` 4910、stream の書き込み 346-384、`BlastHSPMappingInfo` 192-345、`Blast_HSPListsMerge` 2861, 2974 | `mapper/` | E |
| 最終処理 | `core/hspfilter_mapper.c`、`blast_hspstream.c`、`spliced_hits.c` | `s_BlastHSPMapperFinal` 2402、`s_Finalize` 2295、`s_PruneChains` 1937、`s_FindAdapters` 1291、`s_FindPolyATails` 1500、`s_FindRearrangedPairs` 1865、`s_SortChains` 2137、`s_FilterChains` 2261、`BlastHSPStreamMappingClose` 209-225 | `mapper/finalize.rs`、`spliced_hits.rs` | E |
| 結果の組み立て | `api/magicblast.cpp`、`api/blast_seqalign.cpp` | `x_Run` 104-184、`s_ComputeBtopAndIdentity` 208-331、`s_CreateSeqAlign` 334、`x_CreateSeqAlignSet` 393-427、`x_BuildResultSet` 481-533、`CMagicBlastResults` 555-732、`MakeSplicedSeg` 453 | `results.rs` | F |
| 整形 | `app/magicblast_util.cpp`、`include/algo/sequence/consensus_splice.hpp` | `PrintSAMHeader` 848、`PrintSAM` 1008、`PrintSAMUnaligned` 1507、`PrintTabularHeader` 396、`PrintTabular` 441、`PrintTabularUnaligned` 686、`PrintFastaUnaligned` 343、`s_GetSpliceSiteOrientation` 921-994、`s_GetBareId` 156、`s_GetSequenceId` 180 | `format/` | F |

経路の C/C++ 行数（ヘッダ込み）:

| 部分 | 行数 |
|---|---|
| core と api の mapping 経路（`jumper.c` の呼ばれない関数と環境変数で切り替わる部分、約 2,080 行を除く） | 約 12,300 |
| 　うち `hspfilter_mapper.c` | 4,954 |
| 　うち共有ファイル内の分岐（lookup、走査、engine など） | 約 2,200 |
| 入力層（引数、短いリードの FASTA・FASTQ・FASTC） | 約 1,360 |
| アプリ層（`src/app/magicblast/` の 5 ファイル） | 2,692 |
| 計 | 約 16,400 |

索引付き DB 用の `ShortRead_IndexedWordFinder`、`blast_sra_input`、ASN.1 の入出力は
含まない。

段階 A の棚卸しでは、経路上の C++ ライブラリの挙動も移植対象として列挙する。

- `CShortReadFastaInputSource`: ID（タイトルの最初の語）、FASTQ の defline と品質、
  FASTC、ペアの印（Seqdesc の User-object `Mapping`・`has_pair`）、大文字・小文字、
  IUPAC、空のリード、改行コード、複数行の FASTA リード。
- `CFastaReader`（subject）: ID の解析、`-parse_deflines false` の `Subject_N`、
  ギャップを N に変える処理、小文字のマスク。Splign と共通の部分である。
- `CSeq_id` の解析と `BlastRank`、`GetSeqIdString(true)`、`AsFastaString`、
  `SEQ_ID_PREFER_ACCESSION_OVER_GI` の効果。
- `CSeqConvert` と `CSeqManip`（逆相補、IUPAC）。
- ostream による double の書式（tabular の `% identity` は有効数字 6 桁、例 `98.0392`）。
- `log10`（MAPQ）と `log`（二塩基エントロピー）。
- `std::list::sort`（`magicblast.cpp:533`、安定）と、経路上の `qsort`。
- gzip の展開（`CCompressionIStream`）。複数メンバーの gzip の扱いを含む。
- `GetCpuCount` によるスレッド数の頭打ち。
- 例外から終了コードへの対応（`CATCH_ALL`）。

## 6. Rust の構成と再利用の方針

- `LOSAT/src/algorithm/magicblast/`:
  - `args.rs`、`app.rs`（バッチ、スレッド、入力順の出力、終了コード）
  - `input.rs`（短いリードの FASTA・FASTQ・FASTC、ペア、gzip）、`subject.rs`
  - `options.rs`、`setup.rs`、`filter.rs`（リードの品質フィルタ）
  - `lookup.rs`（ハッシュ lookup table、単語数、poly-A）、`scan.rs`、`engine.rs`
    （chunk、word finder、`MapperWordHits`）
  - `jumper.rs`
  - `mapper/`（`hspfilter_mapper.c` を関数のまとまりで分ける: 実行、チェイン、スプライス、
    ペア、最終処理）、`spliced_hits.rs`
  - `results.rs`（`CMagicBlastResults` と、Spliced-seg に相当する構造体。整形で使う
    項目だけを持ち、ASN.1 の汎用オブジェクトは作らない）
  - `format/`（`sam.rs`、`tabular.rs`、`fasta.rs`、`consensus_splice.rs`）
- `LOSAT/src/main.rs`、`cli.rs`: `magicblast` サブコマンド。
- blastn エンジン（`run.rs` など）には手を入れない。共有ファイル（`blast_engine.c`、
  `blast_hits.c`、`blast_setup.c` など）の mapping の分岐は、magicblast モジュールの中に
  NCBI の関数名で移植する。
- 汎用の部品は既存のものを使う: 塩基の符号化と 2 bit 圧縮（`core/blast_encoding.rs`、
  `sequence/packed_nucleotide.rs`）、逆相補、Seq-id の書式（blastn・tblastn の既存
  コード）、スレッドプール（`utils/threading.rs`）、CPU 数での頭打ち（blastn の
  `GetCpuCount` 相当）。
- Splign の実装で `CFastaReader` の意味に合わせる層ができていれば、それを使う。
- gzip の展開には純 Rust の実装（例: `flate2` の `miniz_oxide` バックエンド）を使い、
  純 Rust 境界チェックで依存を確認する。
- 名前は NCBI の関数名に対応させる（AGENTS.md の慣例）。
- 大きなゲノムでは、参照を 2 bit 圧縮（ncbi2na）と曖昧塩基の範囲で一度だけ持ち、
  バッチとスレッドで読み取り専用に共有する。NCBI はバッチごとに subject を作り直す
  （`magicblast_thread.cpp:159-166`）が、出力には影響しない実装上の違いとして扱う。

## 7. 数値と順序で守ること

- **整数のスコア**: スコアは整数で、型の幅は NCBI に合わせる（Int4 など）。
  `GetCutoffScore` の `(Int4)(0.6 * len)`（`jumper.c:4563-4576`）と、`-score L,b,a` の
  `(int)(c * 100)` と整数除算（4.2）を再現する。
- **double の比較**: `JumperGoodAlign` の `100.0 * ident / len < pid`（`jumper.c:2664`）、
  重なりの比率と 0.75 の比較（`hspfilter_mapper.c:3757-3760`）、API の percent
  identity（`magicblast.cpp:330`）。
- **libm**: 二塩基エントロピーは `log` を使い、切り捨てた整数を 16 と比べる
  （`jumper.c:4482-4510`、`blast_options.c:199`）。LOSAT の統計計算と同じく Rust の
  `f64::ln` を使い、閾値の前後の fixture（M6）で差を検出する。MAPQ の `log10` は、
  ヒット数が 500 万未満では整数への丸めの境界から 0.0115 以上離れるので、ulp の差で
  値は変わらない。
- **ソート**（D-9）: 経路上の `qsort` のうち、比較関数が全順序でないのは
  `s_CompareChainsByScore`（`hspfilter_mapper.c:1918`）、`s_CompareChainsByOid`（2114）、
  `s_CompareHSPsByContextScore`（3917）、`s_ComparePairs`（4037）である。全順序なのは
  `s_CompareOffsetPairsByDiagQuery`（`jumper.c:3064`）と
  `s_CompareHSPsByContextSubjectOffset`（`hspfilter_mapper.c:3948`）。LOSAT は
  PD-PURE-RUST-RUNTIME-AUTHORITY に従い、NCBI の比較関数を Rust の安定ソートで使い、
  同点は入力順のままにする。glibc の `qsort` はメモリが足りればマージソート（安定）で
  動くので、一致する見込みである（段階 A でオラクルの環境について確かめる）。fixture で
  同点の並びの差が出たら、比較キーを足さずに分類し、D-11 に従って扱う。
- **到着順**: バッチの中では、HSP list は subject の chunk の順に writer に届く。
  `s_Finalize` は、並べ替えの前に到着順のリストで刈り込み、アダプター、poly-A、
  並び替わったペアを処理する（`hspfilter_mapper.c:2307-2317`）。同点のチェインは
  到着順に足される（165-225）。バッチの中を並列にする場合も、この順序に戻す。
- **バッチの中身**: 1,000 リードを超えるバッチの `MapperWordHits` の経路
  （`blast_engine.c:1054-1056`、`na_ungapped.c:2060-2130`）は、同じ対角線上で隣り合う
  種を捨てる。バッチの切り方を NCBI と同じにする（4.5）。
- **chunk**: subject の chunk の重なりは、最長のリードが 110 未満なら長さの 1.5 倍、
  それ以外は `DBSEQ_CHUNK_OVERLAP`（`blast_engine.c:461-471`）。リードの方が参照より
  長い場合の切り替え（`read_is_query`、`na_ungapped.c:1950`）も再現する。
- **形式に連動する設定**: `MAPPER_NO_OVERLAPPED_HSP_MERGE`（4.4）。
- **統計**: Magic-BLAST は E 値と bit score を計算しない（`blast_setup.c:507-511`）。HSP の
  E 値は 0.0 に固定される（`jumper.c:3156, 3234, 3442`、`na_ungapped.c:2283`）。
  `s_JumperScoreBlkFill` の仮の Karlin ブロックが挙動に効くかを段階 A で確かめ、効く
  部分だけを移植する。
- **例外の握りつぶし**: `CMagicBlast::x_Run` は core 以外の例外を黙って捨てる
  （`magicblast.cpp:171-184`）。第 12 節で扱う。

## 8. テストデータ

fixture、凍結した出力、生成スクリプトは `LOSAT/tests/magicblast/` に置き、manifest は
`LOSAT/tests/magicblast_parity_manifest.tsv` とする。大きな入力は、生成スクリプトと
乱数の種と sha256 だけを保存する。

### 8.1 合成 fixture

| ID | 内容 | 主に通す経路 |
|---|---|---|
| M1 | 単一リード、スプライスなし: 完全一致、ミスマッチ、2 bp の挿入・欠失、両鎖 | jumper 伸長、CIGAR、BTOP、NM |
| M2 | GT-AG のイントロンをまたぐリード（両鎖、短い端）、GC-AG、AT-AC、非コンセンサス | スプライス、`N`、`XS:A`、BTOP の `^gt..ag^` |
| M3 | 数 kb の長いリードと、スコア 50 前後の局所整列での GC-AG・AT-AC | 長いリードの経路 |
| M4 | `-paired`（交互）、`-query_mate`（2 ファイル、gz）、FASTC: concordant、逆向き、片方だけ整列、両方未整列、別の contig | ペアの選択、flag、TLEN、`-no_discordant`、`-fr`/`-rf` |
| M5 | FASTQ（品質付き、マイナス鎖を含む） | QUAL を反転しない点、ID の切り詰め |
| M6 | N や IUPAC が半分を超えるリード、二塩基エントロピーが 16 の前後のリード、poly-A、アダプター付き、15〜35 bp の短いリード | 品質フィルタ（`YF:Z:F`）、`GetCutoffScore` の境界、poly-A とアダプター |
| M7 | 参照に 30 コピー前後の反復配列、多重マップするリード | 単語数フィルタ（`-max_db_word_count` の境界）、NH、MAPQ、0x100 |
| M8 | 複数の subject（いろいろな defline、小文字、N の連続）、`-subject_loc` | `@SQ`、Seq-id、`-lcase_masking`、ギャップを N に変える処理 |
| M9 | 1,000 リードの前後のバッチ、`BATCH_NUM_SEQS` と `BATCH_SIZE` による複数のバッチ、既定の上限を超える入力（60 万リード） | `MapperWordHits`、バッチの切り方、スレッド数 1・4・8 |
| M10 | `-limit_lookup false` と `-reftype transcriptome` を複数スレッドで | バッチの割り算、`BATCH_SIZE` によるオラクル |
| M11 | 4.2 の各オプションを既定値以外にする掃引（`-score L,b,a` を含む） | 引数から設定への反映 |
| M12 | 参照より長いリード、chunk の境界をまたぐリード、chunk が複数ある長い subject | `read_is_query`、chunk の重なり |
| M13 | 空のリード入力、空の subject、ヒットなし、`-out_unaligned` と `-unaligned_fmt` の 3 種、`-no_unaligned` | 未整列の出力、終了コード |
| M14 | クエリ上で HSP が重なるリードを、SAM と tabular の両方で | `MAPPER_NO_OVERLAPPED_HSP_MERGE` |
| M15 | 4.3 の拒否項目と、引数の制約違反 | 拒否の契約、終了コード |

### 8.2 実データ fixture

段階 A で accession.version、取得元 URL、sha256 を記録する。リードは LOSAT の外で
sra-tools により取得し（fixture の準備であり、実行時の依存ではない）、種を固定して
抜き出す。

| ID | 選び方 |
|---|---|
| R1 | 分裂酵母（*S. pombe*、約 12.6 Mb）のゲノムと、RNA-seq のペアエンドリード数万組 |
| R2 | ヒトの染色体 1 本（chr21 など）と、RNA-seq のペアエンドリード |
| R3 | 細菌のゲノムと DNA-seq のリード（`-splice F`） |
| R4 | 長いリード（ONT の cDNA など）の一部 |
| R5 | 転写産物の参照（RefSeq mRNA の一部）と `-reftype transcriptome` |

### 8.3 NCBI の単体テスト

`src/algo/blast/unit_tests/api/magicblast_unit_test.cpp`（565 行、5 ケース）は、BLAST DB
`data/pombe` と ASN.1 のリードを使う。どちらも初期範囲外なので、DB を `blastdbcmd`
（オラクル側のツール）で FASTA にし、リードをペアの情報を保ったまま FASTA にして
fixture U1〜U5 とする。テストの期待値（エキソンの座標、スコア、鎖）は照合に使う。

### 8.4 大きなゲノムでの動作確認（U-5）

GRCh38 の primary assembly 全体を `-subject` に与え、RNA-seq のペアエンドリード 100 万組を
native で実行する。完走したこと、ピークメモリ、時間を記録する。オラクルも同じ入力で
完走すれば、生バイトを比べる。完走しなければ、リードの一部で比べる。

### 8.5 凍結するもの

各ケースについて、argv、環境変数、オラクルの sha256、本体の出力、`-out_unaligned` の
ファイル、終了コード、stderr、H1〜H5 のダンプを凍結する。

## 9. 比較と合格条件

- **最終の合否**: 4.4 の規則で本体の出力と `-out_unaligned` のファイルを生バイトで比べ、
  終了コードの一致を確かめる。コマンドライン文字列の置き換え以外の正規化や並べ替えは
  しない。
- **スレッド**: LOSAT の `-num_threads` 1、4、8 の出力が互いに一致し、NCBI の
  `-num_threads 1`（M10 は `BATCH_SIZE` 付き）と一致する。
- **NCBI の複数スレッドの実行**: 診断として、バッチ単位の並び順の差しかないことを
  確かめる（D-2）。
- **部品段階の合否**: 段階 B〜F は、ハーネスの出力（H1〜H5）との完全一致を条件とする。
  構造化した差分は診断にだけ使う。
- **決定性**: 同じ入力を 3 回実行し、同じ出力になる。
- **既存 gate**: 純 Rust 境界チェック（`LOSAT/tests/check_pure_rust_runtime_boundary.py`、
  gzip の依存の確認を含む）を通す。共有部品に手を入れた場合は、該当する既存 gate も通す。
- **性能**: オラクルと LOSAT の実行時間を、同じスレッド数で参考として記録する。標準の
  手順（ウォームアップ 1 回のあと 3 回計測し、中央値と範囲を報告）に従う。合否には
  使わない。

## 10. 段階別の実装計画

着手は BLASTX v0.2.0 と Splign の後とする（U-7）。各段階の証拠は
`docs/evidence/magicblast_stage_<段階>/` に置く。

### A. 権威・オラクル・棚卸し・fixture

1. PD の記録（3.1）と、AGENTS.md への限定条項の追記。
2. 固定コミットのソースを取得し（Splign と共有）、VDB 付きでオラクルの magicblast を
   ビルドする。ビルドの記録を残す（3.3）。単体テストを実行し、配布バイナリとの差を
   全 fixture で記録する。
3. ハーネス H1〜H5 を作る（3.4）。
4. 経路を棚卸しする。初期範囲のオプションで到達する全関数と C++ ライブラリの挙動
   （第 5 節）を、「移植する／経路外／拒否する」に分類し、ファイル・行番号付きで
   `call_path_inventory.tsv` に記録する。後から監査で見つかった関数を一つずつ足す
   のではなく、ここで全量を洗い出す。
5. 未確定の点を確かめる: `-reftype` の効果、短いリードの読み込みの ID の扱い、
   `GetCpuCount` の意味、オラクルの環境での `qsort` の動き、バッチの中身がリードごとの
   結果に効くかどうか、仮の Karlin ブロックの影響。
6. fixture（第 8 節）を作り、オラクルの出力を凍結する。
7. 第 12 節の決定（U-10〜U-13）を、承認済みの例外の記録と検証の証拠にする。N-3 の
   再現入力を探し、N-7 を確かめる。新しい候補が見つかったら、まとめて判断を仰ぐ。

完了条件: PD と AGENTS.md の更新、オラクルのビルド記録、H1〜H5 の動作、棚卸し表、
凍結した fixture、第 12 節の例外の記録がそろっていること。

### B. 入力・バッチ・設定・品質フィルタ

短いリードの読み込み（FASTA・FASTQ・FASTC、ペア、gzip、標準入力）、subject の読み込み、
バッチの分割（`BATCH_SIZE`、`BATCH_NUM_SEQS`、CPU 数での頭打ち）、オプションと設定、
リードの品質フィルタを移植する。

完了条件: 全 fixture で H1 と一致すること。

### C. lookup table と走査

ハッシュ lookup table、参照の単語数フィルタ、poly-A の語の除去、ハッシュの走査、
subject の chunk、`JumperNaWordFinder`、`MapperWordHits` を移植する。

完了条件: 全 fixture で H2 と H3 (a) に一致すること。

### D. jumper 伸長

`BlastNaExtendJumper` とその下の関数（第 5 節の「jumper 伸長」の行）を移植する。
テスト専用の入口で H3 (a) の種を受け取る。

完了条件: 全 fixture で H3 (b) と一致すること。

### E. HSP writer と最終処理

`hspfilter_mapper.c`、`spliced_hits.c`、stream と HSP list の mapping 部分を移植する。
テスト専用の入口で H3 (b) の HSP list を到着順に受け取る。

完了条件: 全 fixture で H4 と一致すること。

### F. 結果の組み立てと整形

`CMagicBlast` の結果の組み立て（BTOP、MD、percent identity、Spliced-seg 相当、ペア、
concordance、`SortAlignments`）と、SAM・tabular・FASTA の整形を移植する。テスト専用の
入口で H4 のチェインを受け取る。

完了条件: 全 fixture で、本体の出力と `-out_unaligned` のファイルがオラクルと 4.4 の
規則で一致し、H5 とも一致すること。オラクルは同じチェインから整形するので、この段階は
検索と切り離して判定できる。

### G. CLI・スレッド・統合

`LOSAT magicblast` の引数解析（名前、型、制約、既定値、排他）、範囲外の拒否、入力順の
出力を保つ並列化、終了コードを実装し、全段階をつなぐ。

完了条件: 全 fixture で第 9 節の最終の合否とスレッドの条件を満たすこと。拒否の契約、
反復実行での決定性、純 Rust 境界チェック、README・CHANGELOG・`--help` の更新もそろって
いること。

### H. 大きなゲノムと計測

8.4 の動作確認と、第 9 節の性能の計測を行う。

完了条件: GRCh38 全体での完走と記録、計測の記録がそろっていること。

### I. 独立監査と認証

`ncbi_parity_auditor` による独立した確認を受け、`docs/release/magicblast_certification.md`
に認証の範囲（第 4 節）と証拠を記録する。認証の範囲外を互換と主張しない。

## 11. 並行実行

- A を最初に終える。B〜F はどれも新しいモジュールで重ならない。A のハーネスのダンプを
  入力にすれば、別の worktree で並行して進めてよい。
  - B は H1、C は H2 と H3 (a)、D は H3 (a) と H3 (b)、E は H3 (b) と H4、F は H4 と H5 を使う。
- G は B〜F の後、H は G の後、I は H の後に行う。
- blastn エンジンには手を入れないので、BLASTX や他の blastn エンジン作業と直列にする
  必要はない。共有部品（符号化、FASTA、Seq-id、スレッドプール）を変える場合だけ調整する。
- 調査、比較、分類などの機械的な作業は Sonnet のエージェントに任せる（既存の運用方針）。

## 12. NCBI 側の欠陥と扱い（2026-10-04 に決定）

既存の NCBI 欠陥方針（D-11）に従う。NCBI がクラッシュ相当の動きをし、意図された正しい
結果を近い設定の NCBI 出力で確かめられる場合は、承認済みの例外として正しい結果を出す。
結果を出す決定的なバグは、バイト単位で再現する。正しい結果を定義も確認もできない場合は、
明示的に拒否する。計画時点の候補は、ユーザーがすべて推奨案で決定した（U-10〜U-13）。

| ID | 現象 | 根拠 | 決定 |
|---|---|---|---|
| N-1 | worker thread の中の例外（例: `-infmt fastq` に FASTA を与えたときの解析エラー）で、stderr にエラーを出すが終了コードは 0。出力はヘッダ（とそれまでのバッチ）だけ | `magicblast_thread.cpp`、`CThread::Wrapper`。配布バイナリで確認 | 承認済みの例外（U-10）。主スレッドで同じ解析エラーが起きた場合の終了コード 1 を、意図された結果とする。出力のバイト列は NCBI の 1 スレッドの実行と同じにする（ヘッダと、エラーの前のバッチ）。検証の証拠として、主スレッドでの解析エラーの終了コードを記録する |
| N-2 | 書き込みの失敗（`-out /dev/full`）でも終了コード 0 | 配布バイナリで確認 | 承認済みの例外（U-11）。`CATCH_ALL` の出力エラーの終了コード 6 を、意図された結果とする |
| N-3 | `CMagicBlast::x_Run` が core 以外の例外を黙って捨て、バッチのリードが未整列として出る可能性がある | `magicblast.cpp:171-184` | U-12。段階 A で再現する入力を探す。見つかったら改めて判断を仰ぐ。見つからなければ、LOSAT では同じ状況（メモリ確保の失敗など）を明示的なエラーにする |
| N-4 | `-subject_loc` に数値でない値を与えると、捕まらない例外で終了コード 255 | 配布バイナリで確認 | 再現する（U-13、終了コード 255） |
| N-5 | 空の subject ファイルで終了コード 3（`Empty CBlastQueryVector`） | 配布バイナリで確認 | 再現する（U-13） |
| N-6 | flag 16 の記録で QUAL を反転しない（SAM の仕様では反転する） | util:1430-1433。配布バイナリで確認 | 再現する（U-13、4.4） |
| N-7 | `-reftype` が効かない可能性 | args:1566-1577 | 段階 A で確かめ、挙動をそのまま再現する（U-13） |
| N-8 | DB が空のとき、単語数を数える OpenMP 経路で 0 除算 | `blast_nalookup.c:1861` | アプリでは `num_threads > 1` にならず経路外。記録のみ |

承認済みの例外（N-1、N-2）は、段階 A で PD か既存の例外の記録に、検証の証拠とともに
記す。実装中に新しい候補が見つかったら、同じ方針で集めてまとめて判断を仰ぐ。

## 13. 主なリスクと対処

| リスク | 対処 |
|---|---|
| 規模が大きい（約 16,400 行）。中心の `hspfilter_mapper.c` だけで約 5,000 行 | 段階を B〜F に分け、ハーネスのダンプで切り離して判定する。並行して進める（第 11 節） |
| バッチの中身が結果に効く。複数スレッドでは NCBI 自身の出力順が揺れる | バッチの切り方を NCBI と同じにする。合否は 1 スレッドの NCBI と比べる（4.5） |
| `qsort` の同点の並び | PD の方針どおり安定ソートにする。差が出たら分類して判断を仰ぐ（第 7 節） |
| `log` の ulp の差で二塩基エントロピーの判定が変わる | 閾値の前後の fixture（M6）で確かめる |
| ヘルプ文・文書・論文と実装の食い違い（`-score`、既定の語長など） | ソースに従う（U-6）。LOSAT の `--help` には実装どおりの説明を書く |
| オラクルのビルドが難しい（VDB、LMDB、SQLite3 などの依存） | 段階 A で手順を確定して記録する。Splign とビルドを共有する |
| 公式の配布バイナリと自前ビルドの差（同梱ライブラリの版の違い） | 合否は自前ビルドで判定し、差は診断として記録する |
| Seq-id の解析（subject の defline）と `CFastaReader` の細部 | 段階 A で棚卸しし、M8 で固定する。既存の LOSAT のコードと Splign の成果を使う |
| 大きなゲノムでのメモリと時間（NCBI はバッチごとに参照全体の単語数を数える） | 参照は 2 bit 圧縮で一度だけ持つ。8.4 で計測する。性能目標は初期範囲では設けない |
| NCBI の環境変数による開発用スイッチ | 移植せず、オラクルは未設定で実行する（D-8） |
| 範囲の膨張 | 第 14 節の項目は、初期範囲の認証が終わるまで着手しない |

## 14. 初期範囲の後の拡張候補（予定は未定）

- `-db`（BLAST DB の読み込み）と DB 専用のオプション。他のプログラムの `-db` にも使える
  共通の基盤なので、別の計画とする。
- `-infmt asn1`、`-infmt asn1b`、`-outfmt asn`。
- `-gzo`（比較方法の取り決めが先に要る）。
- Wasm（小さなゲノムに限る）と他の native プラットフォーム。
- 性能の改善。

SRA 入力は対象外とする（U-4）。

## 15. 調査の記録（2026-10-04）

- **Magic-BLAST の版**: NCBI FTP（`blast/executables/magicblast/`）には 1.0.0
  （2016-08-22）から 1.7.2（2023-04-19）までがある。1.7.0 で SAM の品質値、MAPQ、BTOP の
  スプライス信号、1.5.0 で BLAST DB v5、MD タグ、長いリード、1.6.0 で SRA のクラウド
  読み込みが入った。GitHub の `ncbi/magicblast` は文書だけで、最終更新は 2025-04-17。
  論文は Boratyn GM ほか, BMC Bioinformatics 20:405 (2019)。
- **配布バイナリでの確認**:
  - 合成データ（300 kb のランダムな参照、100 bp のリード 60,000 本）を 1 スレッドで
    2 回、4 スレッドで 2 回実行し、`@PG` 行を除く SAM が一致した（バッチは 1 つ）。
  - 小さな合成データで、ヘッダ、CIGAR、flag、QUAL、tabular、未整列の出力、
    `-fr`/`-rf`、`-subject_loc`、終了コードを確かめた（第 4 節と第 12 節）。
  - `BATCH_NUM_SEQS=100`、4 スレッドで 3 回実行し、バッチの並びは毎回違ったが、並べ替えた
    中身は同じだった。
- **LOSAT の現状**:
  - CLI は blastn、blastp、tblastx、tblastn。入力は FASTA だけで、FASTQ、gzip、標準入力、
    ペアの概念はない。BLAST DB の読み込みはない。
  - ハッシュ lookup table、jumper、HSP stream と writer の枠組み、SAM の出力はない。
    blastn のパイプラインは `run.rs` に一体で実装されている。
  - Seq-id の書式（blastn、tblastn）と、CPU 数でのスレッドの頭打ちは既存のコードがある。
  - `src/` は約 137,300 行、`algorithm/blastn/` は 29,644 行、`algorithm/tblastn/` は
    20,483 行。
- **調査資料の退避先**（リポジトリ外）: `/home/kawato/losat-magicblast-research-20261004/`
  - `notes/`: 経路、オプションと出力、Toolkit との同一性、オラクルのビルド、LOSAT の
    再利用、外部情報の調査メモ
  - `downloads/`: 1.7.2 のソース配布（sha256 `f1a5dbd7…`）、1.7.2 の Linux 配布
    （sha256 `93301d18…`）、固定コミットの Toolkit（sha256 `0dc02d6c…`）、`help.txt`
  - `comparisons/`: 1.7.2 と 2.17.0、Toolkit の比較結果
  - `toy/`: 合成の入力と配布バイナリの出力
  - `SHA256SUMS`
