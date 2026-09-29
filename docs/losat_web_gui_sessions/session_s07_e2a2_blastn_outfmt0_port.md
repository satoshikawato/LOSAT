# Session S07 — E2a-2：BLASTN outfmt 0（実装とゲート）

## INSTRUCTION PROMPT

LOSAT の段階 E2a-2 を実行する。S06 で固定した権威と fixture に従って、BLASTN（LOSATN）に outfmt 0 を移植する。先に [セッション README](README.md) の共通規則を読み、特に規則 4（`AGENTS.md`、`verify-ncbi-parity-and-speed`）に従う。完了条件の正本は、総合計画書 §7 の S07 の行である。入力は `docs/evidence/losat_web_e2a/`（S06 の成果物）である。

1. S06 の経路の表の順に移植する。Rust の移植箇所の直上に、NCBI のファイル・行と C/C++ の断片を書く。`LOSAT/src/report/pairwise.rs` の既存の部品を使い、同じ働きの writer を二重に作らない。
2. BLASTN の最終の HSP 一覧（並べ終えたもの）から `PairwiseHit` を作る。整列文字列は、traceback の編集操作（`gap_info`）から作る。この一覧を `run_local` の `hits` と、outfmt 0 の formatter の両方に渡す（計画 §4.4）。観測者を outfmt 0 の節にもつなぐ。
3. `parse_blastn_output_format` が 0 を受け付けるようにする。未対応の組合せ（S06 の一覧）は明示的に拒否する。
4. 試験：S06 の fixture すべてで、NCBI とバイト一致すること。最初に違ったところは、query / subject の ID、鎖、座標、score、E 値、HSP の順序、表示の文字列に分けて記録し、該当するものをまとめて直す。既存の BLASTN の outfmt 6/7 の回帰ゲートと Gate A のハッシュが変わらないこと。`docs/web/verification_cells.tsv` の BLASTN の outfmt 0 の升目を埋め、BLASTN の全升目（0/6/7 × スレッド 1/2/4）の V-ABI が通ること。`cargo fmt --check`・`clippy`・`cargo test --all-features`。
5. 独立監査を受ける。

完了条件は計画 §7 の S07 の行による。記録は `docs/evidence/losat_web_e2a/README.md`（S06 の記録の後に S07 の節を足す）。

### S05 の独立監査の残り（アダプタ。最初に行う）

S07 は reactor を作り直して full の V-ABI を再実行するので、S05 の独立監査の指摘のうち reactor のバイトを変えるものをここで直す（記録は `docs/evidence/losat_web_e1d/README.md` の「独立監査」）。

1. 計画 TD-11：`web/adapter/tools/build_reactors.py` で、checkout のパスと `CARGO_HOME` を `--remap-path-prefix` で置き換え、`tools/check_build_identity.py` はこの置き換えだけを rustflags の差として認める。別の場所に展開した同じ commit から同じバイトの reactor ができることを確かめ、`docs/web/abi_v2.md` §2 を直す。
2. `validate`：解析に渡す binary の名前を CLI と同じ `LOSAT` にする。`-num_threads` を `-out`・`-outfmt` と同じく拒否する（host が付ける）。`tests/v_abi.js` の `checkSurface` に、必須の引数が無い場合と `-num_threads` の場合を足す。
3. `register`：記録の照合を関数に分け、ID・長さ・記録の数が食い違う場合に失敗することの単体試験を足す（S05 の指示書 4）。
4. `RangeRecorder`：対応の無い `hsp_end` と重複した番号を、release のビルドでもエラーにする（今は `debug_assert` だけ）。
5. `tools/run_v_abi_parallel.py`：部分を実行する前に、その部分の古い結果を消す。`summary.json` に 3 つの成果物の SHA-256 を記録する。

### S06 の結果（移植の順とゲート）

権威は `docs/evidence/losat_web_e2a/AUTHORITY.md`、fixture は `LOSAT/tests/outfmt0_manifest.tsv`（34 件。期待する出力は `LOSAT/tests/fixtures/outfmt0/`）である。LOSAT は `LOSAT/` から、manifest と同じ引数で実行する（`docs/evidence/losat_web_e2a/run_oracle.py` の `search_argv`）。`(Input: …)` の行は `-subject` の綴りをそのまま含むので、作業ディレクトリと綴りを変えない。

移植の順（`AUTHORITY.md` §F。各項目の NCBI の根拠は §A〜§D）：

1. **BLASTN の panic を直す**（§D.2）。`blast_get_start_for_gapped_alignment_nucl` の後ろ向きの走査を subject の先頭で、前向きの走査を subject の長さで止める（NCBI は subject の両端の番兵で止まる。`blast_gapalign.c:3340-3349`、`seqsrc_multiseq.cpp:353-358`、`blast_encoding.c:121`）。S06 で試した差分は `docs/evidence/losat_web_e2a/run-20260929T000046Z/probe-gapped-start-bound.diff`。独立したエンジンの変更とし、回帰試験は `multi.ws7_e1000.blastn` の入力（`-task blastn -word_size 7`）で outfmt 6 が NCBI と一致すること。BLASTN の既存のゲート（`compare_blastn_parity.py`、Gate A）が変わらないこと。
2. **座標の桁数の関数**（§D.1）。`showalign.cpp:1366-1367` の規則（両方の行の 0 始まりの start / stop の最大値の桁数）の関数を 1 つ作り、BLASTP の `write_alignment_with_sequences` と TBLASTN の `write_tblastn_alignment` から使う。BLASTX の `write_blastx_alignment` は変えない（DW-10）。`width.blastp` と `width.tblastn` が一致すること。`capture_outputs.py` で、ハッシュが変わる凍結出力をすべて列挙し、それぞれが NCBI の出力と一致することを示す（一致しなければ止める）。
3. `write_alignment_with_sequences` の引数化：`eBar` の中線、0 始まりの桁数、subject の向き（±1）、空の行の規則、マスクの前の文字での比較。
4. `write_hsp_info` の ` Strand=` の行。program は両方の task で `"blastn"` とする（`"megablast"` は Positives を出してしまう）。
5. query の末尾の統計の Gumbel の列を省けるようにする。
6. 末尾（epilog）の `Matrix: blastn matrix R P` と、`double` の gap extension（0 のときの式）。
7. `write_translated_pairwise_intro` に文献の引数（megablast は eMegaBlast）。
8. `-max_target_seqs` が与えられたかを区別する（§D.4。BLASTX の `blastx/args.rs:915-922` の方法）。
9. BLASTN の driver：最終の HSP 一覧から `PairwiseHit`、説明とアラインメントの数での切り詰め、query と subject のマスク、minus 鎖の反転、無効な query の末尾、観測者の `hsp_begin` / `hsp_end` / `subject_begin` / `subject_end`。
10. stderr の警告（§D.3。`Examining 5 or more matches is recommended` と無効な query の警告。どの形式でも出す）。BLASTX の関数は `[blastx]` を直接書くので呼ばない。

ゲート：BLASTN の 32 件で stdout がバイト一致し、`stderr_sha256` のある 3 件は stderr も一致する。`width.blastp` と `width.tblastn` が一致する。BLASTN の outfmt 6/7 の回帰ゲートと Gate A・Stage G のハッシュは、2 の変更（NCBI と一致することを示したもの）を除いて変わらない。V-ABI の BLASTN の升目（0/6/7 × スレッド 1/2/4）が通る（`web/adapter/src/run.rs` の `Program::formats`、`web/adapter/tests/v_abi.js` と `web/adapter/tools/v_abi_cases.py` の `FORMATS` に BLASTN の 0 を足す）。

範囲：fixture の得点は認証済みの既定の範囲（blastn 2/−3・5/2、megablast 1/−2・0/0）と、outfmt 6 が一致した `multi.r2p3g00.megablast` に限る。それ以外の得点の組合せは S07+ で扱う（§D.5）。`AUTHORITY.md` §F の「拒否を続けるオプション」は拒否のままにする。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S07+ — BLASTN の得点のオプション](session_s07p_e2c_blastn_scoring_options.md)。BLASTN で作った部品のうち TBLASTX に使えるものを S08 の指示書に書き足す。
