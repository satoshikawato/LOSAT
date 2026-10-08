# Session S14 — W5：抽出と候補

## INSTRUCTION PROMPT

LOSAT Web の段階 W5 を実行する。先に [セッション README](README.md) の共通規則を読み、それに従う。完了条件の正本は、総合計画書 §7 の S14 の行である。設計は計画 §5.7、設計書 §11.3〜§11.4 にある。

1. 座標の変換を `web/app/src/domain/` の 1 か所に置く（1 始まりで両端を含む座標、strand、単位 nt / aa）。outfmt 6 の座標（start > end は minus 鎖）から、元のレコードの区間を得る。
2. 配列の回収：ヒット区間、左右を別々に指定する flank、Subject の全長、選んだ候補の multi-FASTA。既定は元のレコードの向き。端で切り詰めたときは、要求した範囲と実際の範囲を両方記録する。核酸の長さとタンパク質の残基数は単位を示す。原配列は Data worker が元の File から読み、大文字小文字を保つ。
3. 同じ Subject の複数の HSP は、別々の配列にするか、間を含む 1 区間にするかを選べる。重なりを勝手にまとめない。CDS の復元やスプライシングの推定はしない。
4. ギャップ付きアラインメントの書き出しは、HSP レコードの整列文字列から作り、原配列の抽出とは区別する。逆向きのヒットは、検索で得た向きのまま示す。
5. 候補トレイ：複数の実行から候補を集め、元の結果（run、`q_idx`、`rank`）へ戻れる。メモ、並べ替え、削除、由来の一覧。候補の登録・抽出・書き出しは、完了した実行だけで行える（`REQ-10`。アプリ層で拒否する）。
6. 試験：領域を指定した実行、端、逆向き、翻訳（BLASTX・TBLASTN・TBLASTX）の各場合で、抽出した配列が元のファイルの配列と一致すること。重複 ID のレコードを取り違えないこと。

完了条件は計画 §7 の S14 の行による。記録は `docs/evidence/losat_web_w5/README.md`。

## 開始の条件と最初の作業

- 場所と手順：Linux の clone の worktree `/home/kawato/losat-work/.worktrees/web-gui-app`（ブランチ `feature/losat-web-gui-app`）。場所は clone の `CLAUDE.local.md`、手順は skill `losat-campaign`・`losat-worktree`・`losat-gates`・`losat-ship-pr` に従う（README の規則 1・5 の `/mnt/c` のパスは `CLAUDE.local.md` の表で読み替える）。
- 開始の条件：エンジン側が W4（このブランチ）を `feature/losat-web-gui` に merge し、SF（E2h、`CFastaReader` の移植）がそのブランチに入っていること（DW-23 (6)：抽出は SF の読み方に乗る）。満たさないときは始めず、保守者に伝える。
- 最初に `git merge origin/feature/losat-web-gui` を行い（衝突の解消は独立したコミット）、reactor とネイティブの CLI をその木から作り直す（`$BUILD_ROOT/s14-reactors`・`$BUILD_ROOT/native`）。FakeEngine の `describe.json` は `LOSAT_WEB_REACTORS=<dir> npx vitest run tests/unit/engine-runtime.test.ts -u` で書き直す。
- S13+（E2j、subject ごとの集約の値と HSP の鎖）が先に merge されていたら、その指示書の 6.（アプリへの申し送り）を最初に行う：列定義表の「採用・エンジン待ち」の列を NCBI の表形式の並びで出し、鎖の欄で BLASTN の 1 文字の HSP の向きを示す。

## S13（W4）から引き継ぐこと

[W4 のゲート記録](../evidence/losat_web_w4/README.md)の要点：

- **結果画面の構成**：`src/application/results.ts` の `ResultsBrowser` が、選択（run、query、subject、HSP）・表示用のフィルター（ViewState）・並べ替えを持つ。選択の中心は HSP の ID（`runId`・`qIdx`・`rank`、`HspId`）で、表の行番号ではない。候補トレイ（5.）の「元の結果へ戻る」は、この ID で `ResultsBrowser.open(runId)` と `selectQuery`・`selectSubject`・`selectHsp` を呼べばよい。
- **データの読み方**：`RunStore.readHitTable(runId)` は HSP レコードを列の形（`src/domain/hsp-table.ts` の `HspTable`、typed array を Data worker から transfer）で返し、整列文字列（`query_aligned`・`subject_aligned`）を含まない。ギャップ付きアラインメントの書き出し（4.）は `readHits(runId)` のレコードの整列文字列から作る。outfmt 6 の行と outfmt 0 の節は `readOutputRange(runId, format, start, end)` で、HSP の `out6`・`out0`・`out0_subject` の範囲を読む。
- **座標と向き**：表示用の向きは `src/domain/result-index.ts` の `orientation`（`forward`・`reverse`・`unknown`）、frame は `hsp-table.ts` の `frame`（0 は無し）。単位は `src/domain/programs.ts` の `residueUnit`。1.（座標の変換を 1 か所に置く）は、これらと重ねずに domain にまとめ、結果画面もそれを使うように寄せる。BLASTN の 1 文字の HSP は start = end で向きが座標から分からない（S13 は「not in the record」と outfmt 0 の `Strand=` を示す。S13+ が鎖の欄を足すまで、抽出では向きを決めずに理由を示す）。
- **E2E の補助**：検索の入力を作って Run を完了させる補助は `tests/e2e/support/search.ts`（`paste`・`openFiles`・`settled`・`program`・`submit`・`waitStatus`・`result`）。結果画面の操作は `tests/e2e/results.spec.ts` を見る。キューの「Open results」（`run-N-open`）で完了した Run の結果を開ける。
- **FakeEngine**：検索の形をした出力（200 query まで、3 subject まで、3 つ目の subject は outfmt 0 に無い、逆向きの HSP、BLASTN の 1 文字の HSP）を書く。値は `FAKE`。FakeEngine の outfmt 6 は先頭に印の行がある（行は範囲で読む）。
- **検証バッジ**：`build/verification.ts` が `docs/web/verification_cells.tsv` と認証の記録から表を生成する（手で書かない）。認証の記録を足せば、表は次のビルドで変わる。
- **Run の入力**：`Coordinator.runInputs` は解放しない（Run の削除は S13 で入れなかった。W4 の判断 5）。候補トレイの削除は Run の削除ではない。
- **画面の記録と画面レビュー**：`tests/e2e/screens.spec.ts` が検索画面と結果画面の状態を撮る。W4 の記録は `$BUILD_ROOT/s13-gate-screens/`（SHA-256 は W4 の実行記録の `screens.sha256`）。W5 の画面レビューは、検索画面と結果画面が変わらないことも確かめる。
- **ゲートの script**：W4 の [`run_gate.sh`](../evidence/losat_web_w4/run_gate.sh) を元にできる。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S15 — 出力と再現性](session_s15_w6_export_session.md)。
