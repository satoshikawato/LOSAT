# Session S09+ — R2：TBLASTN の subject の前処理（条件付き）

## INSTRUCTION PROMPT

LOSAT Web の段階 R2 を TBLASTN について実行する。先に [セッション README](README.md) の共通規則を読み、それに従う。完了条件の正本は、総合計画書 §7 の S09+ の行である（「連続実行の出力が CLI と一致。独立監査」）。設計は計画 §4.6 の R2、DW-8、§2.3 の PreparedSubject の行にある。エンジン側のセッション（`LOSAT/`・`web/adapter/` を変える）で、ブランチは `feature/losat-web-gui`、作業場所は worktree `/mnt/c/Users/genom/GitHub/LOSAT-web-gui` である。

**入口の条件**：S09（W1）で、ブラウザの warm 実行のうち subject だけで決まる前処理の割合が、TBLASTN で 20% 以上になった（DW-8 の条件。[W1 のゲート記録](../evidence/losat_web_w1/README.md)の「DW-8」）。ほかの program（BLASTP・BLASTN・TBLASTX）は 20% 未満で、R2 の対象ではない。

### S09 の実測（W1 のゲート記録から）

- 測った検索：V-PERF の TBLASTN の fixture（`docs/evidence/tlosan_stage_g/benchmark/first_AvCLPV_protein.faa`、1 本の 1046 aa、`-task tblastn`、subject `LOSAT/tests/fasta/AvCLPV.fasta`、416,069 nt）。
- ブラウザ（Chromium、serial、warm、1 回の暖機と 3 回の計測の中央値）：検索の段階は 92 ms、同じ subject に query の先頭 30 残基だけを当てた検索は 52 ms で、比（subject だけの前処理の上限）は 0.56（4 スレッドでは 0.71）。
- ネイティブの callgrind（命令数）：subject だけの前処理は全体の 40.9%（下限 36.1%、上限は約 50%）。内訳は 6 フレームの翻訳 `tblastx::translation::generate_frames`（36.1%）と、subject の ncbi2na への変換 `resolve_local_subject_ncbi2na`（4.8%）。呼ぶ場所は `LOSAT/src/algorithm/tblastn/search_seed.rs:428`（変換）と `:448`（翻訳）、1 回の検索で subject ごとに 1 回で、検索をまたいだキャッシュは無い。
- 翻訳の命令数の 76% は、codon ごとに曖昧な塩基の組合せをたどる `GeneticCode::get`（`LOSAT/src/utils/genetic_code.rs:188-231`、codon あたり約 159 命令）である。64 通りの表引きと曖昧な塩基の時だけの元の関数にすれば、この割合は約 16% になる見込み（推定。測っていない）。

### 作業

1. **翻訳の表引き化（出力を変えない）**：`GeneticCode::get` を、遺伝暗号ごとに曖昧でない 64 codon の表と、曖昧な塩基を含む codon だけ今の関数を使う形にする。翻訳の結果は今の関数とバイト単位で同じでなければならない（全 27 の遺伝暗号 × 全 15^3 の IUPAC codon で、前後の関数の結果を比べる単体試験を置く）。TBLASTX と BLASTX も同じ関数を使うので、変える前と後の全 program の出力の SHA-256（`docs/evidence/losat_web_e1a/capture_outputs.py`）、TLOSAN の Stage G のゲート、Gate A、V-PERF の非退行を確かめる。
2. **測り直す**：W1 の `web/app/tests/e2e/measure.spec.ts`（`LOSAT_WEB_MEASURE=dw8`）と callgrind で、TBLASTN の subject だけの前処理の割合を測り直す。**20% 未満になれば、DW-8 の条件を満たさないので、キャッシュは入れない**（作業 3〜5 を行わず、結果を記録して終える）。
3. **キャッシュ（20% 以上のままのとき）**：アダプタの登録済みのハンドル（R1 の `register`）に、TBLASTN の subject の前処理の結果（ncbi2na に変換した塩基と、`-db_gencode` ごとの 6 フレームの翻訳）を持たせ、`run` の検索がそれを使うようにする。query に依存する構造（lookup など）は持たない（計画 §4.6）。NCBI の拠り所は、subject 集合から `CLocalDbAdapter` を一度作り、query の batch ごとに `CLocalBlast` を作り直す構成（計画 §4.9）。ハンドルの解放で捨てる。メモリの増え方（subject の長さの約 3 倍）を記録する。
4. **ゲート**（計画 §4.6）：同じ subject に対して query・条件・遺伝暗号を変えて連続で実行したとき、どの出力も、毎回新しく実行した CLI の出力と一致すること。V-ABI（Subject を保持した連続実行を含む）、V-BR の R1 の試験（`web/app/tests/e2e/engine.spec.ts`）に TBLASTN の連続実行を足す（アプリの変更は、アプリ側のセッションに頼むか、S09 の試験の形のまま case を足す）。
5. 独立監査（README の「レビュー」）、V-PERF の非退行、ゲート記録 `docs/evidence/losat_web_r2_tblastn/README.md`。

### 保守者の判断

作業 1 と 2 の順（表引き化で割合が 20% を下回れば、キャッシュを入れない）は、W1 のゲート記録で保守者に諮った推奨の案である。保守者が別の案（作業 1 をせずにキャッシュを入れる、など）を選んだときは、それに従う。

## 終了・引き継ぎ

README の規則 8 に従う。
