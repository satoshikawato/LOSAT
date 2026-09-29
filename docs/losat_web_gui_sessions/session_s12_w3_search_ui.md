# Session S12 — W3：検索画面

## INSTRUCTION PROMPT

LOSAT Web の段階 W3 を実行する。先に [セッション README](README.md) の共通規則を読み、それに従う。完了条件の正本は、総合計画書 §7 の S12 の行である。設計は計画 §5.3〜§5.5、設計書 §3.1〜§3.2 と §3.4（画面の骨格）にある。

始める前に、program の表示名（BLASTN 系か、LOSATN 系か、併記か）を保守者に確認する（計画 §10）。決まらなければ、今の `web/app/src/domain/programs.ts` の表示名のまま進め、ゲート記録に「判断待ち」と書く。

1. program のタブと、Query / Subject の入力（貼り付け、複数ファイル、ドラッグ＆ドロップ）。大きなファイルは textarea に流し込まず、ファイル名・配列数・総長・先頭の表示・レコード一覧・警告を見せる。
2. レコード一覧と、明示的な除外。重複 ID は内部の番号で区別する。不正なレコードの判定は、エンジンの結果を表示する。配列種別の警告はアプリの推定として出し、program は変えない。
3. 領域の指定（座標の入力とプレビューでの選択）。その役割のレコードが 1 つのときだけ指定でき、S11 で移植した `-query_loc` / `-subject_loc` の argv になる（DW-9）。
4. パラメーターのフォーム。`ProgramDescriptor` に、表示する引数・節・ラベルを足し、既定値・選択肢・help は `describe` から取る（初期値は BLAST+ CLI の既定値）。既定値と同じ値は argv に書かない。遺伝暗号の選択肢は、program ごとの許可リスト（`describe`）から出す。
5. Combined / Separate の切り替え（計画 §5.2）。Separate は、ファイルごとの RunSnapshot を同じグループ ID で積み、「グループを取消」を用意する。
6. スレッドの Auto / 手動、実行中の段階と経過時間、診断情報（RunRecord の経路・スレッド数・切り替えの理由）。query ごとの途中経過は出さない（DW-5）。
7. モバイルでの縦の配置。Wake Lock の選択、実行中にタブを閉じる前の注意、スリープやバックグラウンドからの復帰時の状態の確認（設計書 §9.4、`REQ-22`）。
8. BLASTN の得点と入力（S07+、`docs/evidence/losat_web_e2c/AUTHORITY.md`）。フォームは次を守る：
   - `-reward`・`-penalty`・`-gapopen`・`-gapextend`・`-word_size` は省略できるオプションで、省いたときは task の既定値（blastn：11、2、−3、5/2。megablast：28、1、−2、0/0）になる。task を切り替えたときに、既定値を argv に書き込まない（既定値と同じ値を明示しても NCBI と同じ結果になるが、argv は短く保つ）。
   - reward は 1 以上（0 は LOSAT が拒否する）、penalty は 0 以下、gap は任意の整数。NCBI が拒否する組合せ（penalty 0、gap extend 0 の gap open、blastn の task の gap 0/0、NCBI の Karlin の表に無い組）は、`validate` が NCBI の文言（`BLAST query/options error: …` または `BLAST engine error: Error: …`）で返すので、その文言をそのまま見せる。表の組と、表を超える gap（ungapped の block を写す）は、NCBI とバイト一致で実行される（`scoring_sweep.py` の 880 の組合せ）。
   - word size は 4 以上（100 を超えると `validate` が NCBI の文言で拒否する）、e-value は 0 より大きい。LOSAT は、reward − penalty が 3000 を超える得点と、megablast の 32767 を超える gap を、LOSAT の文言で拒否する（`AUTHORITY.md` §G）。
   - `register` のハンドルは、登録した program の `run` でだけ使える（BLASTN の入力の検査を他の program の登録で避けられないように）。program を切り替えたら登録し直す。
   - `register` は、NCBI が違う読み方をする BLASTN の入力（空の定義行・先頭の空白・tab などの制御文字・非 ASCII、残基の無いレコード、IUPAC の文字以外の残基。`U` は `T` として受け付ける）を 「not supported by LOSAT's BLASTN」を含む文言で拒否する（TD-12）。レコード一覧の警告として、その文言を見せる。
   - query の batch に依存する場合（表を超える gap で組成の違う query が最初の batch に収まらない、など）は、`run` が LOSAT の文言で失敗する。失敗として見せ、結果を部分的に出さない。
9. E2E：代表的な研究作業（Subject を保持したまま Query を変えて繰り返す、キューに複数積む、実行中に次のジョブを編集する）と、境界条件（空の入力、不正なレコード、除外の後の再実行、取消）を Playwright で試す。

完了条件は計画 §7 の S12 の行による。画面の記録を画面レビューに見せ、指摘と対応を `docs/evidence/losat_web_w3/README.md` に記録する。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S13 — 結果画面](session_s13_w4_results_ui.md)。
