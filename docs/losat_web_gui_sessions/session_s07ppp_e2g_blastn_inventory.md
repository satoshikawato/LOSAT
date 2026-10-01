# Session S07+++ — E2g：BLASTN の経路の棚卸しと一括の移植

## INSTRUCTION PROMPT

LOSAT の段階 E2g を実行する。BLASTN（LOSATN）の NCBI との一致を、独立監査の指摘を 1 つずつ直すのでなく、NCBI の実行経路の棚卸しと一括の transpile で仕上げるセッションである（計画 DW-12）。先に [セッション README](README.md) の共通規則を読み、特に規則 4（`AGENTS.md`、`verify-ncbi-parity-and-speed`）に従う。完了条件の正本は、総合計画書 §7 の S07+++ の行である。

背景：S07+ では 16 回の独立監査のたびに、NCBI の経路の未移植（小さい lookup table の word の延長、subject の曖昧な文字の `CRandom`、query の batch と分割など）と、C の細部（`Int4` の変換、`size_t` の折り返し、`get_data` の位置のずれ、得点 0 の HSP の除去）が 1 つずつ見つかった（`docs/evidence/losat_web_e2c/AUTHORITY.md` §H〜§R）。

1. **棚卸し。** アプリが出す BLASTN のオプション（`describe` と S12 の指示書の範囲）で、NCBI の `blastn`（`-subject` の bl2seq）の実行経路に現れる関数を洗い出す：`blastn_app.cpp` の引数の処理と batch、`CLocalBlast`、`CBlastPrelimSearch`（query の分割を含む）、`algo/blast/core` の setup・lookup・scan・ungapped・gapped・traceback・hit の保存・統計、`blast_seqalign.cpp`、`blast_format.cpp` と `align_format` の outfmt 0/6/7。関数ごとに、NCBI のファイル・行、LOSAT の対応箇所、状態（忠実な移植、差のある移植、未移植、明示的な拒否、承認済みの例外）を表にする（`docs/evidence/losat_web_e2g/INVENTORY.tsv`）。LOSAT のコードの `NCBI reference` の注釈を機械的に突き合わせて、分類を始める（その script も証拠に置く）。option の値ごとに分岐が変わる関数（例：`BlastChooseNaExtend`、`s_SmallNaChooseScanSubject`）は、分岐ごとに行を分ける。
2. **一括の transpile。** 「未移植」と「差のある移植」を、NCBI の関数ごとに、簡略化せずに移植する（直上に NCBI のファイル・行と断片）。LOSAT が速度のために NCBI と違う実装にしている箇所は、出力が同じなら新しく移植する部分にも同じ方式を使ってよい（DW-12）。そのときは表に書く。C の細部（整数の型の幅と折り返し、浮動小数点から整数への変換、ポインタの位置）も NCBI と同じにする。
3. **明示的な拒否の見直し。** S07+ と S07++ の明示的な拒否（`AUTHORITY.md` §G など）のうち、transpile で NCBI と同じにできるものは拒否をなくす。残すものは表に理由を書く（NCBI が落ちる、NCBI の C++ の層の丸ごとの移植が要る、など）。
4. **試験。** S07+ と S07++ のすべての検査（`check_inputs.py`、`scoring_sweep.py`、`word_size_sweep.py`、`slice_sweep.py`、`title_sweep.py`、`batch_sweep.py`）、既存の BLASTN のゲート（Gate A）、outfmt 0 の fixture に退行なし。移植した関数ごとに、その分岐に入る入力を NCBI と比べる。V-PERF の非退行。独立監査は、棚卸しの表を基準に、経路の網羅と移植の忠実さを確かめる。

記録は `docs/evidence/losat_web_e2g/`（`README.md`、`INVENTORY.tsv`、`evidence.sha256`）。変更前の成果物は S07++ の成果物である（ゲートの実行ファイルのハッシュは `docs/evidence/losat_web_e2f/run-20260930T163727Z/artifacts.sha256`、native は `~/.cache/losat-web-gui-target/s07pg-native/release/LOSAT`、SHA-256 `a4ff5abb…80f5`）。

## S07++ からの引き継ぎ（S07++b、2026-10-01 の実測）

S07++（E2f）は完了した（[ゲート記録](../evidence/losat_web_e2f/README.md)、`AUTHORITY.md` §A〜§E）。独立監査は第 1 回が unsupported（予備の hit list の大きさ `prelim_hitlist_size` を traceback の後に適用していた）、第 2 回が supported（325 件の比較で差 0、`docs/evidence/losat_web_e2f/audit_round2/`）。これを前提に、次の点を守る。

### 進め方（保守者の指示、2026-10-01）

- **オーケストレータは Opus（メインのセッション）。** 設計を伴う作業だけを行う：NCBI の移植の設計と実装、監査の指摘の判断、分類の結果の確認と統合、agent の成果の確認。機械的な作業は Agent ツールに `model: "sonnet"` で回す：ゲートや sweep の実行と集計、NCBI との比較の実行、記録・`evidence.sha256`・`verification_cells.tsv` の行の下書き、決まった観点でのファイルの走査と分類。
- **並行は同時に 4 つまで。** 互いに独立な読み取り専用の作業（第 2 段の分類の 7 範囲、独立監査の観点ごとの確認、NCBI との比較の実行）だけを並行させる。**エンジンのソース（特に `run.rs`）を変える作業は serial に 1 つずつ行う**（同じファイルを触る・ビルドの出力先が競合する）。ビルドの出力先は agent ごとに別の `--target-dir` にする。ゲートの実行は 1 つだけ（V-PERF は lock を取る）。アプリ側の S09（`/mnt/c/Users/genom/GitHub/LOSAT-web-gui-app` に中断した作業が未コミットで残っている）はこのセッションでは触らない。
- **使用量の上限で作業が失われないようにする。** 2026-09-30 に並行した agent が同時に上限に達して成果を失った。agent には、結果を**途中で、ファイルに書きながら**進めさせる（分類は範囲ごとの TSV を `docs/evidence/losat_web_e2g/` か scratch に行ごとに追記する。比較の実行は case ごとの結果を TSV に書く）。agent の指示は自己完結にして、成果の保存先を指定する。agent が終わるたびにメインが読んで確かめ、コミットする。上限に達したら、並行数を 2 に減らして続ける。
- V-PERF で閾値を超えた case の切り分けの再計測は `perf_cases.py run … --repeat 5`。V-PERF の段階はアプリ側の lock（`~/.cache/losat-web-gui-target/s07p-resume/vperf_lock.sh`）を取る。ゲートの script の雛形は `~/.cache/losat-web-gui-target/s07p-resume/s07pp_gates2.sh`（実行の前に `date -u +%Y%m%dT%H%M%SZ > …/s07p-resume/ts2.txt` で run の名前を決める。無ければ `run-20260930T163727Z/` の記録から作り直す）。NCBI の参照の注釈の検査は `…/s07p-resume/verify_refs.py <変えた .rs>`。
- 独立監査の第 2 回の script（`docs/evidence/losat_web_e2f/audit_round2/{gen,cases,driver}.py`）は、棚卸しの後の監査でそのまま再利用できる。

### 棚卸しの第 1 段（やり直し）

NCBI の参照の注釈（`NCBI reference: …` と直後の断片）を機械的に突き合わせる script を `docs/evidence/losat_web_e2g/` に置く（前回の一時的な結果は失われた）。前回は、BLASTN から届く 54 の Rust ファイルで 2180 の注釈を数え、669 の NCBI 関数に対応した。この数を再現できることを、script の最初の検査にする。

### 分類（第 2 段）

7 つの範囲ごとに 1 つずつ、sonnet の agent に回す（独立な読み取り専用の作業なので、4 つまで並行してよい）。各 agent は行を範囲ごとの TSV に追記しながら進め、メインが読んで統合する。

| 範囲 | NCBI のファイル |
|---|---|
| A：app と引数 | `blastn_app.cpp`、`blast_args.cpp`、`CFastaReader` |
| B：API | `CLocalBlast`、`prelim_stage`、`setup_factory`、`blast_setup_cxx`、`seqsrc_multiseq`、`traceback_stage`、分割 |
| C：統計と DUST | `blast_setup.c`、`blast_parameters.c`、`blast_stat.c`、`blast_filter.c`、`dust_filter.cpp`、`symdust.cpp` |
| D：lookup・scan・ungapped | `blast_nalookup.c`、`blast_nascan.c`、`na_ungapped.c`、`blast_extend.c` |
| E1：gapped | `blast_engine.c`、`blast_gapalign.c`、`greedy_align.c` |
| E2：traceback と hit の保存 | `blast_traceback.c`、`blast_hits.c`、`hspfilter_collector.c`、`blast_hspstream.c`、`blast_itree.c` |
| F：整形 | `blast_seqalign.cpp`、`blast_format.cpp`、`showalign.cpp`、`tabular.cpp`、`create_defline.cpp` |

表の列と状態は上の 1.（`INVENTORY.tsv`）のとおり。これに加えて、行ごとに**影響の大きさ**（high / medium / low / none）を付けさせる（出力が変わりうる入力の広さ。例：既定のオプションで出る差は high）。

### 一括の transpile の確かめる候補

設計と実装はメイン（Opus）が行う。前回の調べで、次が確かめる候補になっている。

- ほかの `qsort` を Rust の安定な並べ替えにした箇所。オラクルの glibc は 2.39 で、`qsort` は安定な merge sort。
- `run.rs` の初期の hit の並べ替え。`sort_unstable_by` を使っている。
- 予備の hit list（`74223fee3`）の単体試験の不足。S07++ の独立監査の第 2 回の指摘（重大度 低）：`hsp.rs` の `test_hitlist_update_keeps_best` は大きさ 1・異なる e-value の場合だけで、subject の番号による同点の順、e-value・得点が等しい場合、1e-180 の規則、2 つ以上の list の heap、collector、`merge_prelim_hit_list`、最後の hit list の大きさの上限を固定していない。試験だけの変更でもソースの行が動いて実行ファイルが変わるので、このセッションの最初のエンジンの変更と同時に足す。
- 残る明示的な拒否（S07++ の `AUTHORITY.md` §C）：得点の表が無い得点で、最初の batch が無効な query だけのとき。NCBI は、その batch の結果と警告を書いた後に、次の batch で誤りを出す。棚卸しの表で、拒否のままにするか移植するかを決める。
- 証拠の範囲：S07++ の比較は EDL933 と Sakai だけで、`BATCH_SIZE` を与えた NCBI との比較は行っていない（LOSAT は `BATCH_SIZE` を拒否する）。棚卸しの後の独立監査では、別の生物の配列（ウイルス、IUPAC の多い配列）も入力に含める。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S08 — TBLASTX outfmt 0/7](session_s08_e2b_tblastx_outfmt0_7.md)。S08+ は同じ棚卸しの方式で行う（S08+ の指示書）。
