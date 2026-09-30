# Session S06 — E2a-1：BLASTN outfmt 0（権威と fixture）

## INSTRUCTION PROMPT

LOSAT の段階 E2a-1 を実行する。BLASTN（LOSATN）に outfmt 0 を移植する前に、何を NCBI と一致させるかを固定するセッションである。コードは変えない。先に [セッション README](README.md) の共通規則を読み、特に規則 4 に従う。完了条件の正本は、総合計画書 §7 の S06 の行である。

現状：BLASTN は outfmt 6/7 だけを受け付け、0 を拒否する（`LOSAT/src/algorithm/blastn/hsp.rs:30-43`）。CLI の既定値は "0" なので、`-outfmt` を省くとパースで失敗する（`LOSAT/tests/cli_v2.rs`）。BLASTP・TBLASTN・BLASTX の outfmt 0 は `LOSAT/src/report/pairwise.rs` にある。

1. NCBI の経路を記録する：`c++/src/app/blast/blastn_app.cpp` から `CBlastFormat`（`c++/src/algo/blast/format/blast_format.cpp` の `PrintProlog`、`PrintOneResultSet`、`PrintEpilog`）、`c++/src/objtools/align_format/` の `showdefline.cpp`（説明の一覧）と `showalign.cpp`（アラインメント）へ。核酸に特有の表示（`Strand=Plus/Minus`、identities と gaps、positives が無いこと、小文字のマスクの表示、長い defline の折り返し、megablast と blastn の task の違い、ヒットが無い場合、query ごとの見出しと末尾の統計）を、ファイルと行の番号付きで表にする。
2. `LOSAT/src/report/pairwise.rs` の既存の関数のうち、BLASTN に使えるものと核酸に固有に書くべきものを分けて記録する（重複した writer を新しく作らない）。
3. fixture を選ぶ。LOSAT の BLASTN の対応範囲（task megablast / blastn、局所 `-subject`）の中で、plus / minus 鎖、1 subject に複数 HSP、複数の query と subject、ヒット無し、`-lcase_masking`、`-dust`、ギャップ、長い defline、`-max_target_seqs` の境界を含める。入力は `LOSAT/tests/blastn_parity_manifest.tsv` の既存の入力をなるべく使う。
4. NCBI BLAST+ 2.17.0 の `blastn -outfmt 0` を比較用の oracle として実行し、出力を SHA-256 とともに固定する。実行ファイルの場所と版、コマンド、入力の SHA-256 を記録する（既存の比較スクリプト `LOSAT/tests/run_comparison.sh` と `LOSAT/tests/compare_blastn_parity.py` が使う設定に合わせる）。
5. S07 のゲート（固定した fixture すべてでバイト一致、既存の 6/7 に退行なし）と、対応しないオプション（明示的に拒否するもの）の一覧を書く。

成果物は `docs/evidence/losat_web_e2a/AUTHORITY.md`（経路の表）、fixture の manifest（TSV）、固定した NCBI の出力と `evidence.sha256`。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S07 — BLASTN outfmt 0：実装とゲート](session_s07_e2a2_blastn_outfmt0_port.md)。経路の表から、S07 で移植する関数の一覧と順番を S07 の指示書に書き足す。
