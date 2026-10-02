# Session S08 — E2b：TBLASTX outfmt 0 と 7

## INSTRUCTION PROMPT

LOSAT の段階 E2b を実行する。TBLASTX（TLOSATX）に outfmt 0 と 7 を移植する。先に [セッション README](README.md) の共通規則を読み、特に規則 4 に従う。完了条件の正本は、総合計画書 §7 の S08 の行である。1 セッションで終わらなければ、S06・S07 と同じく「権威と fixture」（S08a）と「実装とゲート」（S08b）に分け、README の表を直す。

現状：TBLASTX は outfmt 6 だけを受け付ける（`LOSAT/src/blastinput/value_parsers.rs:313-320`）。最終の結果は平らな `Vec<Hit>` で、subject frame・整列文字列を持たない（`LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs` と `LOSAT/src/common.rs`）。承認済みの例外がある：局所 `-subject` の検索では、既定以外の `-db_gencode` を subject の翻訳・検索・表示に使う（ルートの `AGENTS.md` の「Mandatory Compliance Requirements」の 3）。

1. NCBI の経路を記録する：`c++/src/app/blast/tblastx_app.cpp` から `CBlastFormat` と `c++/src/objtools/align_format/`（`showalign.cpp`、`showdefline.cpp`、表形式の注釈は `tabular.cpp`）へ。翻訳検索に特有の表示（`Frame = +1/-2` のような両側の frame、ギャップの無いアラインメント、positives、翻訳した配列の行、核酸座標）と、outfmt 7 の query ごとの見出しと末尾を、ファイルと行の番号付きで表にする。
2. fixture を選ぶ（`LOSAT/tests/tblastx_v010_parity_manifest.tsv` の入力を基に、6 frame、複数 HSP、複数の query と subject、ヒット無し、遺伝暗号 1 と既定以外のコード）。NCBI BLAST+ 2.17.0 の outfmt 0 と 7 を oracle として固定する。既定以外の subject の遺伝暗号の場合は、承認済みの例外のとおりに差を分類する（ほかの差は認めない）。
3. 最終の HSP 一覧から `PairwiseHit` を作るのに必要な情報（subject frame、翻訳した配列の区間）を、NCBI と同じ時点で保持する。計算の順序は変えない。`PairwiseHit` を `hits` と outfmt 0 の formatter の両方に渡し、観測者を outfmt 0 と 7 につなぐ。
4. outfmt 0 と 7 を移植し、`tblastx_outfmt` が受け付けるようにする。S07 で分かった NCBI の振る舞いと、使える部品は `docs/evidence/losat_web_e2a/AUTHORITY.md` §G にある（特に §G.3 の説明の一覧表の規則。TBLASTX は sum statistics の `N` の列を持つ。BLASTP と TBLASTN もこの規則に切り替えるかを、ここで決める）。Rust の移植箇所の直上に NCBI の参照を書く。S07 の部品（座標の桁数の関数、`eBar` 以外の中線を含む整列の writer、`docs/evidence/losat_web_e2a/run_oracle.py` と `LOSAT/tests/outfmt0_manifest.tsv`。TBLASTX の fixture はこの manifest に足す）を使う。subject の見出しと説明の一覧の title は、NCBI の `CDeflineGenerator::GenerateDefline` で作る（TBLASTX の subject は核酸なので、S07+ の `LOSAT/src/report/defline.rs` の `ncbi_nucleotide_title` と、HTML の文字参照の拒否を使う。`docs/evidence/losat_web_e2c/AUTHORITY.md` §N）。S07+ で分かった入力の読み方の拒否（§E・§M・§N）のうち、TBLASTX にも当てはまるものは S08+ で扱う。
S07+ で分かった、TBLASTX にも当てはまる NCBI の振る舞い（`docs/evidence/losat_web_e2c/AUTHORITY.md`）：NCBI は Seq-align を作るとき、subject ごとの HSP を e-value（同じなら `ScoreCompareHSPs`）で並べ直す（`blast_seqalign.cpp:1571-1577`）。TBLASTX は context（frame）ごとに Karlin block を持つので、同じ得点の HSP の順序が frame で変わりうる。outfmt 0/7 の HSP の順序をこの規則で確かめる。検索されなかった query の batch には outfmt 7 の「# N hits found」が無い（`tabular.cpp:1264-1284`）。NCBI と同じ終了コードと文言の失敗は `LOSAT/src/cli.rs` の `NativeError` で返す（S07+ で BLASTX と共有にした）。

5. 試験：fixture すべてで NCBI とバイト一致（例外の分類を除く）。既存の outfmt 6 の回帰ゲート（`LOSAT/tests/audit_tblastx_v010.py` と Gate A のハッシュ）が変わらない。`docs/web/verification_cells.tsv` の TBLASTX の outfmt 0/7 の升目を埋め、TBLASTX の全升目（0/6/7 × スレッド 1/2/4）の V-ABI が通ること。`cargo fmt --check`・`clippy`・`cargo test --all-features`。
6. アダプタ：TBLASTX の形式に 0 と 7 を足す（`web/adapter/src/run.rs` の `Program::formats`、`tests/v_abi.js` と `tools/v_abi_cases.py` の `FORMATS`）。これで全 program が既定の `-outfmt 0` を受け付けるので、`run.rs` の `parse` が挿入している `-outfmt 6` をやめ、`validate` の文言を CLI と同じにする（S05 の独立監査の指摘。`docs/web/abi_v2.md` §4 を直す）。
7. 独立監査を受ける。

完了条件は計画 §7 の S08 の行による。記録は `docs/evidence/losat_web_e2b/README.md`。

## S07+++b からの引き継ぎ（2026-10-02 の実測）

S07+++b（E2g）は BLASTN の棚卸し（`docs/evidence/losat_web_e2g/INVENTORY.tsv`、1015 行）の一括の transpile を終え、独立監査は第 4 回で supported になった（[ゲート記録](../evidence/losat_web_e2g/README.md)）。毎晩の WASI の TIMEOUT（`threaded/tblastx/p11_avclpv_psclpv`）は runner の差で、TBLASTX の速度の退行ではない（手元の 4 コアで HEAD 相当 2164 秒、緑の run の成果物 2809 秒）。期限を 7200 秒にした（`beea0bddf`）ので、S08 の最初の作業にはしない。S08 で使えること、守ること：

- **進め方：** オーケストレータは Opus（設計を伴う作業：移植の設計と実装、監査の指摘の判断、agent の成果の確認）。機械的な作業（NCBI との比較の実行、sweep、ゲート、記録の下書き、決まった観点の走査）は Agent ツールに `model: "sonnet"` で回す。sonnet の agent は互いに独立なら 4 つを超えて並行してよい。エンジンのソースの変更、ゲート、V-PERF は 1 つずつ。agent には結果を途中でファイルに書かせ、作業ディレクトリと `--target-dir` を分ける。
- **棚卸しの方式（DW-12）：** S08 以降も、NCBI の経路の関数を棚卸しして（`docs/evidence/losat_web_e2g/` の `inventory_refs.py`・`build_inventory.py`・`stage2/` の方式）、未移植と差のある移植を一括で移す。
- **NCBI の不具合の扱い（DW-15、`docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md`）：** NCBI が落ちる・debug でだけ防ぐ前提が崩れる・buffer の外を読む入力は、近い入力の NCBI の出力と一致を示せる妥当な結果を LOSAT が出せるなら承認済みの例外（保守者に諮る）。NCBI が決まった結果を出すものは、見かけが誤りでも再現する。妥当な結果を確かめられない失敗は明示的な拒否。NCBI が受け付ける入力は移植する。そのような挙動が見つかったら、まとめて保守者の判断を仰ぐ（保守者の指示）。TBLASTX の outfmt 0 の題も BLASTN と同じ `ncbi_nucleotide_title` で、句読点だけの題は文字列の終わりで止める（例外 2）。
- **項目ごとの検査：** `~/.cache/losat-web-gui-target/s07p-resume/e2g_item_check.sh <項目> <変えた .rs>` が、pure-Rust の境界の検査（CI の `rust` のジョブと同じ `LOSAT/tests/check_pure_rust_runtime_boundary.py`）、`verify_refs.py`、fmt、clippy、`cargo test --all-features`（`opt-level=1`）、release のビルド、`ci_fast_regressions.py` をまとめて行う。`extern "C"` の import はこの検査が拒否する。
- **NCBI の凍結 fixture：** `LOSAT/tests/blastn_regression_fixtures.py` は case ごとの環境変数（列 `env`）、NCBI だけの環境（列 `oracle_env`：NCBI が失敗する承認済みの例外の期待値を、近い設定の NCBI の出力で凍結する）、stderr を stdout に合わせた case（case の名前の末尾 `.merged`）を持てる。TBLASTX でも同じ方式の fixture を作り、変更前の実行ファイルでも `check` して、case がその分岐を区別することを記録する。
- **警告と報告の順（BLASTN、`e37099f44`・`f4057718a`）：** NCBI は警告を cerr に出し、cerr は cout に tie されているので、警告の前に書いた分の stdout が出る（`2>&1` で観測できる）。outfmt 0 の prolog は最初の query の batch を読む前に書いて flush し、書き込みの失敗はそこで止まる。LOSAT は `report/query_warnings.rs` の `QueryWarnings`（各 query の報告の前で、報告を flush してから警告を書く）と、`run.rs` の `write_pairwise_prologs` で同じにした。TBLASTX の outfmt 0/7 を作るときも、stdout と stderr を 1 つにした NCBI との比較を入れる。NCBI の flush の位置は gdb で write の system call を数えると分かる（`docs/evidence/losat_web_e2g/order/flushmap.py`）。
- **outfmt 0 の書き込みの失敗：** 「BLAST failed to write output」、終了コード 6（`/dev/full`）。閉じたパイプは承認済みの例外 5（LOSAT は同じ文言と終了コード 6、NCBI は SIGPIPE）。TBLASTX の outfmt 0 も同じ扱いにする。
- **`CTOOLKIT_COMPATIBLE`：** 共有の outfmt 0 の説明の一覧の見出しは「(bits)」になる（`report/pairwise.rs`）。TBLASTX の outfmt 0 を作ったら、変数を与えた NCBI との比較に TBLASTX を足す（S07+++b では TBLASTX に outfmt 0 が無く、比べられなかった 4 件）。
- **`%#8.3g`：** Lambda・K・H の書き方は C と同じになった（全 program 共有）。
- **NCBI の application の層：** `LOSAT/src/blastinput/ncbi_environment.rs` の `check_ncbi_application_settings(program)` が、`DIAG_*`・`NCBI_CONFIG_*` などと、出力を変える `<program>.ini`・`.ncbirc` を拒否する。今は blastn だけが呼ぶ（`main.rs`）。S08+ でほかの program の入口にも足す。
- **整数でない環境変数の値：** 明示的に拒否する（計画 TD-15）。
- **保守者の確認を待つ細部（BLASTN、ゲート記録の「残件」）：** 負の `CHUNK_SIZE` の組で batch を分けるがもう一度は分けない場合（拒否のまま）、reward 32768 以上（拒否のまま）、outfmt 6/7 の閉じたパイプの時間に依る終了コード 0。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S08+ — BLASTP・TBLASTN・TBLASTX の既定以外のオプション](session_s08p_e2e_protein_options.md)。
