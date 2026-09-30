# Session S07+++ — E2g：BLASTN の経路の棚卸しと一括の移植

## INSTRUCTION PROMPT

LOSAT の段階 E2g を実行する。BLASTN（LOSATN）の NCBI との一致を、独立監査の指摘を 1 つずつ直すのでなく、NCBI の実行経路の棚卸しと一括の transpile で仕上げるセッションである（計画 DW-12）。先に [セッション README](README.md) の共通規則を読み、特に規則 4（`AGENTS.md`、`verify-ncbi-parity-and-speed`）に従う。完了条件の正本は、総合計画書 §7 の S07+++ の行である。

背景：S07+ では 16 回の独立監査のたびに、NCBI の経路の未移植（小さい lookup table の word の延長、subject の曖昧な文字の `CRandom`、query の batch と分割など）と、C の細部（`Int4` の変換、`size_t` の折り返し、`get_data` の位置のずれ、得点 0 の HSP の除去）が 1 つずつ見つかった（`docs/evidence/losat_web_e2c/AUTHORITY.md` §H〜§R）。

1. **棚卸し。** アプリが出す BLASTN のオプション（`describe` と S12 の指示書の範囲）で、NCBI の `blastn`（`-subject` の bl2seq）の実行経路に現れる関数を洗い出す：`blastn_app.cpp` の引数の処理と batch、`CLocalBlast`、`CBlastPrelimSearch`（query の分割を含む）、`algo/blast/core` の setup・lookup・scan・ungapped・gapped・traceback・hit の保存・統計、`blast_seqalign.cpp`、`blast_format.cpp` と `align_format` の outfmt 0/6/7。関数ごとに、NCBI のファイル・行、LOSAT の対応箇所、状態（忠実な移植、差のある移植、未移植、明示的な拒否、承認済みの例外）を表にする（`docs/evidence/losat_web_e2g/INVENTORY.tsv`）。LOSAT のコードの `NCBI reference` の注釈を機械的に突き合わせて、分類を始める（その script も証拠に置く）。option の値ごとに分岐が変わる関数（例：`BlastChooseNaExtend`、`s_SmallNaChooseScanSubject`）は、分岐ごとに行を分ける。
2. **一括の transpile。** 「未移植」と「差のある移植」を、NCBI の関数ごとに、簡略化せずに移植する（直上に NCBI のファイル・行と断片）。LOSAT が速度のために NCBI と違う実装にしている箇所は、出力が同じなら新しく移植する部分にも同じ方式を使ってよい（DW-12）。そのときは表に書く。C の細部（整数の型の幅と折り返し、浮動小数点から整数への変換、ポインタの位置）も NCBI と同じにする。
3. **明示的な拒否の見直し。** S07+ と S07++ の明示的な拒否（`AUTHORITY.md` §G など）のうち、transpile で NCBI と同じにできるものは拒否をなくす。残すものは表に理由を書く（NCBI が落ちる、NCBI の C++ の層の丸ごとの移植が要る、など）。
4. **試験。** S07+ と S07++ のすべての検査（`check_inputs.py`、`scoring_sweep.py`、`word_size_sweep.py`、`slice_sweep.py`、`title_sweep.py`、`batch_sweep.py`）、既存の BLASTN のゲート（Gate A）、outfmt 0 の fixture に退行なし。移植した関数ごとに、その分岐に入る入力を NCBI と比べる。V-PERF の非退行。独立監査は、棚卸しの表を基準に、経路の網羅と移植の忠実さを確かめる。

記録は `docs/evidence/losat_web_e2g/`（`README.md`、`INVENTORY.tsv`、`evidence.sha256`）。変更前の成果物は S07++ の成果物である。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S08 — TBLASTX outfmt 0/7](session_s08_e2b_tblastx_outfmt0_7.md)。S08+ は同じ棚卸しの方式で行う（S08+ の指示書）。
