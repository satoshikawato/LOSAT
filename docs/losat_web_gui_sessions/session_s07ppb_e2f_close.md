# Session S07++b — E2f の仕上げ：独立監査の第 2 回とゲート記録

## INSTRUCTION PROMPT

LOSAT の段階 E2f（S07++、BLASTN の query の batch と分割）を仕上げる。先に [セッション README](README.md) の共通規則を読み、それに従う。完了条件の正本は、総合計画書 §7 の S07++ の行である。前の指示書は [S07++](session_s07pp_e2f_blastn_query_batches.md)、権威の記録は `docs/evidence/losat_web_e2f/AUTHORITY.md`（§A〜§E）。

### このセッションの進め方（保守者の指示、2026-10-01）

- **使用量を抑える。** 機械的な作業は Agent ツールに `model: "sonnet"` で回す。例：ゲートや sweep の実行と結果の集計、NCBI との比較の実行、ゲート記録・`evidence.sha256`・`verification_cells.tsv` の行の作成（下書き）、決まった観点でのファイルの走査。メインのセッション（Opus）は、設計を伴う作業だけに使う：NCBI の移植の設計、監査の指摘の判断、修正の設計と実装、agent の成果の確認。
- **一度に 1 つ、serial に。** agent は前景（`run_in_background: false`）で 1 つずつ実行し、自己完結した指示を渡す。agent が書いたものは、メインが読んで確かめてからコミットする。エンジン側とアプリ側を並行させない（アプリ側の S09 は、このセッションの後に別のセッションで行う）。
- V-PERF で閾値を超えた case の切り分けの再計測は `perf_cases.py run … --repeat 5`（10 回ではない）。

### 現状（2026-10-01、ブランチ `feature/losat-web-gui`）

- **S07++ のコミット：**
  - エンジン：`93c58a2e4`（batch）、`bbb869a9e`（query の分割）、`74223fee3`（予備の段階の hit list）；
  - 文書：`a74be99f4`、`b951afeb2`、`7693a9c73`；
  - ゲートの実行記録：`4a5371293`（`docs/evidence/losat_web_e2f/run-20260930T163727Z/`、`head.txt` は `7693a9c73`）。
- **ゲートの結果（`7693a9c73`）：すべて通過。**
  - fmt、clippy（エンジン 4 構成、アダプタ 3 構成）、試験 830 件。
  - fixture 47 件（`-num_threads` 1・2・4、NCBI の再実行もハッシュと一致）。
  - 得点の sweep（outfmt 0・6・7 の各 same 300 / same-error 580）、word size の sweep（256 件・差 0）、切り出しの sweep（17 の設定・差 0）、title の sweep（957 / 66）。
  - batch の sweep（6 の設定・差 0。S07+ の実行ファイルは `-task blastn` で 1 件違う）。
  - `check_inputs.py`（294 件すべて期待どおり）、`split_check.py`（213 件・差 0。そのうち NCBI の分割で出力が変わる 15 件）。
  - capture（236 件・差 0）、v1 の検査（S07+ と同じ 21 件の文言の差）、V-ABI（full 460 件・凍結ハッシュ 672/676、quick 52 件）、ビルドの同一性。
  - V-PERF：最初の計測で 15 件すべて閾値の内。分割の時間：EDL933 × Sakai[1000000:1500000] の `-task blastn` で NCBI 1.56〜1.64 秒、LOSAT 1.51〜1.71 秒、出力は同じ。
- **独立監査の第 1 回（`b951afeb2` で実施）：unsupported。**
  - 指摘 1（高）：予備の段階の hit list の大きさ（`prelim_hitlist_size`）を traceback の前に守っていなかった。分割しない batch にも以前からあった差で、`74223fee3` で移植した（AUTHORITY §E）。
  - 指摘 2（低）：`split_check.py` が megablast の分割を実行していなかった。2 つの query は NCBI では 2 つの batch になるので、`7693a9c73` で直した。
  - 注記：分割した batch の lookup table を NCBI は作らない（`blast_aux_priv.cpp:206-207`）。`74223fee3` で注釈を直した。
  - 修正の後、第 1 回の再現の入力（下）と、第 1 回で一致した入力（約 170 件、全ゲノム 6 件を含む）は、すべて NCBI とバイト一致した。ただし、その入力のスクリプトは一時ディレクトリにあり、失われた。
- **残り：**
  - 第 1 回の再現を `split_check.py` に残す；
  - 独立監査の第 2 回；
  - ゲート記録（`README.md`、`evidence.sha256`）；
  - `docs/web/verification_cells.tsv` の S07++ の行；
  - README の表、計画の状態、次の指示書；
  - `main` への PR と merge（保守者は merge を承認済み）。
- **アプリ側：** `/mnt/c/Users/genom/GitHub/LOSAT-web-gui-app` に、中断した S09 の作業が未コミットで残っている（25 ファイル、試験していない、push していない）。このセッションでは触らない。

### 作業（順に、1 つずつ）

1. **第 1 回の再現を回帰の case にする（sonnet）。** `docs/evidence/losat_web_e2f/split_check.py` に、次の 4 つを決まった種で作る case として足す。
   - 期待：修正の前の実行ファイル `/home/kawato/.cache/losat-web-gui-target/s07pp-split1-native`（`bbb869a9e` の作業ツリーのビルド）では NCBI と違い、HEAD のビルドでは一致する。どちらも確かめて、記録を新しい `run-<UTC>/` に置く（既存の run は書き換えない）。
   - (a) query は EDL933[0:2,000,000]（2 つの塊）。subject は、境をまたぐ B = q[999200:1000800] と、q[0:990000) から取った 1500 文字の 560 本（置換 3%）。既定の設定、`-task blastn -outfmt 6`。NCBI は 500 の subject を出し、B を出さない。
   - (b) 同じ query。subject は X = EDL933[997000:1003000] と、Y0〜Y9 = EDL933[50000 + 80000·i : +4500]。`-max_target_seqs 1` と `3`（対照として `10` と `11`）。
   - (c) query は EDL933 と Sakai をつないだ 1 つのレコード。subject は X = つないだ配列[5510500:5516500] と、10 本の Y（つないだ配列の遠い位置の 4500 文字）。`-task megablast -max_target_seqs 1`。
   - (d) 分割しない batch。query は EDL933[0:300000]。subject は X = EDL933[100000:102000] + 25 文字の乱数 + EDL933[102000:104000] と、Y0〜Y9 = EDL933[150000 + 12000·i : +3000]。`-max_target_seqs 1`。
2. **独立監査の第 2 回（sonnet、読み取り専用）。** 役割は README の `ncbi_parity_auditor`。範囲は、第 1 回の後の変更（`74223fee3`、`7693a9c73`）と、それが触る経路である。
   - collector（`hspfilter_collector.c:83-161`）と `Blast_HitListUpdate` の heap（`blast_hits.c:3150-3300`、`s_EvalueCompareHSPLists`、`Blast_HSPListSortByEvalue`）。
   - 予備の e-value の値と、合わせた HSP の e-value。
   - `Blast_HitListMerge`（`blast_hits.c:2119-2217`）。
   - traceback が読む順と、最後の hit list の大きさ（`Blast_HSPResultsInsertHSPList` の `hitlist_size`）。
   - `-max_target_seqs`、`-subject_besthit`、`-max_hsps` との組合せ、`-num_threads 4`、分割する batch と分割しない batch。
   - 確かめ方：subject の数が `prelim_hitlist_size` を超える入力を作り、NCBI BLAST+ 2.17.0（`/home/kawato/micromamba/bin/blastn`、比較だけ）とバイト比較する。
   - 結論は supported / unsupported / inconclusive。指摘があれば、メイン（Opus）が NCBI のソースで判断し、修正を設計して実装する。そのあとゲートと監査をやり直す。
3. **ゲート**（エンジンを変えたときだけ。sonnet が実行と集計）。スクリプトは `~/.cache/losat-web-gui-target/s07p-resume/s07pp_gates2.sh`。実行の前に `date -u +%Y%m%dT%H%M%SZ > ~/.cache/losat-web-gui-target/s07p-resume/ts2.txt` で run の名前を決める。V-PERF の段階はアプリ側の lock（`vperf_lock.sh`）を取る。NCBI の参照の注釈の検査は `~/.cache/losat-web-gui-target/s07p-resume/verify_refs.py <変えた .rs>`。
4. **ゲート記録（sonnet が下書き、メインが確かめる）。** `docs/evidence/losat_web_e2f/README.md` を `docs/evidence/losat_web_e2c/README.md` と同じ構成で書く：段階、ブランチ、コミット、権威の記録、実行記録、判定、変更の内容、完了条件と結果、独立監査、性能の計測、残件と扱い。`docs/evidence/losat_web_e2f/evidence.sha256` も作る（リポジトリのルートで `sha256sum --check` が通るもの）。`docs/web/verification_cells.tsv` には S07++ の升目を足す（batch の sweep、`split_check.py`、`check_inputs.py`、fixture）。
5. **README の表と計画。**
   - README の表：S07++ を「完了」にし、この S07++b の行を完了にする。
   - 計画の状態を更新し、次の指示書 [S07+++](session_s07ppp_e2g_blastn_inventory.md) を実測に合わせて直す（下の「S07+++ への引き継ぎ」）。
   - コミットして push し、`main` への PR を作って merge する（S07+ は PR #107）。
6. 最終回答は README の規則 8 に従う。

### S07+++ への引き継ぎ（次の指示書に書くこと）

- **棚卸しの第 1 段（やり直し）。** NCBI の参照の注釈を機械的に突き合わせる script を、`docs/evidence/losat_web_e2g/` に置く。前回の一時的な結果は失われた。前回は、BLASTN から届く 54 の Rust ファイルで 2180 の注釈を数え、669 の NCBI 関数に対応した。
- **分類（第 2 段）。** 7 つの範囲ごとに 1 つずつ、sonnet の agent に serial に回す：
  - A：app と引数（`blastn_app.cpp`、`blast_args.cpp`、`CFastaReader`）；
  - B：API（`CLocalBlast`、`prelim_stage`、`setup_factory`、`blast_setup_cxx`、`seqsrc_multiseq`、`traceback_stage`、分割）；
  - C：統計と DUST（`blast_setup.c`、`blast_parameters.c`、`blast_stat.c`、`blast_filter.c`、`dust_filter.cpp`、`symdust.cpp`）；
  - D：lookup・scan・ungapped（`blast_nalookup.c`、`blast_nascan.c`、`na_ungapped.c`、`blast_extend.c`）；
  - E1：gapped（`blast_engine.c`、`blast_gapalign.c`、`greedy_align.c`）；
  - E2：traceback と hit の保存（`blast_traceback.c`、`blast_hits.c`、`hspfilter_collector.c`、`blast_hspstream.c`、`blast_itree.c`）；
  - F：整形（`blast_seqalign.cpp`、`blast_format.cpp`、`showalign.cpp`、`tabular.cpp`、`create_defline.cpp`）。

  表の列と状態は S07+++ の指示書のとおり。行ごとに影響の大きさ（high / medium / low / none）も付けさせる。
- **一括の transpile の設計と実装はメイン（Opus）が行う。** 前回の調べで、確かめる候補が出ている：
  - ほかの qsort を Rust の安定な並べ替えにした箇所。オラクルの glibc は 2.39 で、qsort は安定な merge sort。
  - `run.rs` の初期の hit の並べ替え。`sort_unstable_by` を使っている。
