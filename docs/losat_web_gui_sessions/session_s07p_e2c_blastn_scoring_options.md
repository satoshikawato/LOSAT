# Session S07+ — E2c：BLASTN の得点のオプション

## INSTRUCTION PROMPT

LOSAT の段階 E2c を実行する。BLASTN（LOSATN）の得点のオプション（`-reward`、`-penalty`、`-gapopen`、`-gapextend`）のうち、認証済みの既定の組合せ以外を、NCBI と同じにするか、明示的に拒否するセッションである。先に [セッション README](README.md) の共通規則を読み、特に規則 4（`AGENTS.md`、`verify-ncbi-parity-and-speed`）に従う。完了条件の正本は、総合計画書 §7 の S07+ の行である。

現状（S06 の `docs/evidence/losat_web_e2a/AUTHORITY.md` §D.5）：認証済みの BLASTN の fixture は、各 task の既定の得点（blastn：reward 2、penalty −3、gap 5/2。megablast：reward 1、penalty −2、gap 0/0）だけを使う。`docs/evidence/losat_web_e2a/scoring_sweep.py` で、10 の reward / penalty の組 × 4 の gap の設定 × 2 つの task を outfmt 6 で比べた結果（`docs/evidence/losat_web_e2a/run-20260929T000046Z/scoring-sweep-edf825d36.tsv`）：

- NCBI が拒否する 34 の組合せを、LOSAT は受け付けて実行する。NCBI の文言は `Gap existence and extension values … are not supported for substitution scores …`（`c++/src/algo/blast/core/blast_stat.c:3916-3919`）と `Greedy extension must be used if gap existence and extension options are zero`（`c++/src/algo/blast/core/blast_options.c:1701-1709`）。
- NCBI が受け付ける 17 の組合せで、HSP か得点が違う。どれも reward 1 で gap が 0 でない。
- LOSAT には NCBI の Karlin の表がある（`LOSAT/src/core/blast_stat/lookup_tables.rs`）ので、原因は表の外にある。

これらのオプションは、ABI v2 の `describe` を通してアプリの検索画面（S12）に出る。

1. `scoring_sweep.py` を現在の commit で実行し、S06 の結果と比べて現状を確かめる。
2. **NCBI の検査を移植する。** NCBI の `blastn` がこれらの組合せを拒否する経路（オプションの検査、`Blast_KarlinBlkNuclGappedCalc` の組合せの表、task ごとの伸長の方法の検査）をソースで追い、同じ条件で、NCBI と同じ文言で拒否する。NCBI の stderr の形と終了コードは oracle で確かめる（`AGENTS.md` の fail-fast の規則）。
3. **違いの原因を調べる。** reward 1 と 0 でない gap の組合せで、NCBI が得点と統計の値をどう決めるか（表の選び方、得点の倍率、X-drop の換算など）をソースで追い、最初に値が食い違う箇所（lambda、K、H、生の得点、X-drop、`eff_searchsp`）を記録してから直す。直す箇所の直上に NCBI のファイル・行と断片を書く。
4. 直した組合せは、fixture として `LOSAT/tests/outfmt0_manifest.tsv` に足し、`docs/evidence/losat_web_e2a/run_oracle.py` で固定する（既存の行のハッシュが同じになることを、先に照合の実行で確かめる）。outfmt 6/7 の比較は `scoring_sweep.py` と、必要なら `LOSAT/tests/blastn_parity_manifest.tsv` の行で行う。
5. 直せなかった組合せは、LOSAT が対応していないことを示す文言で明示的に拒否する。黙って NCBI と違う結果を出さない。拒否がアダプタの `validate` にも出ることを確かめる。
6. 試験：`scoring_sweep.py` の全組合せが、NCBI と同じ拒否か、outfmt 6 のバイト一致か、明示的な拒否のどれかになる。BLASTN の既存のゲート（`compare_blastn_parity.py`、Gate A）と S07 の outfmt 0 の fixture が変わらない。エンジンの変更の V-PERF の非退行。`cargo fmt --check`・`clippy`・`cargo test --all-features`。独立監査。

記録は `docs/evidence/losat_web_e2c/`（`README.md`、`evidence.sha256`、変更前と変更後の sweep の結果）。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S08 — TBLASTX outfmt 0/7](session_s08_e2b_tblastx_outfmt0_7.md)。対応する得点の組合せ（直したものと拒否するもの）を [S12 の指示書](session_s12_w3_search_ui.md) に書き足す。
