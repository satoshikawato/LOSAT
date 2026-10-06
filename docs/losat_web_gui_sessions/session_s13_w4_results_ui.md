# Session S13 — W4：結果画面

## INSTRUCTION PROMPT

LOSAT Web の段階 W4 を実行する。先に [セッション README](README.md) の共通規則を読み、それに従う。完了条件の正本は、総合計画書 §7 の S13 の行である。設計は計画 §4.4〜§4.5、§5.7、§6.1、設計書 §10〜§11 にある。

1. 列定義表を作り、保守者の確認を受ける。項目は、名称・単位・対象の program・値の出どころ（outfmt 6 の列、outfmt 0 の見出し、エンジンのレコード）・集約の単位・採否。Total score や Query cover のように NCBI Web にしか無い集約列を採用する場合は、NCBI の `align_format`（`c++/src/objtools/align_format/showdefline.cpp`）からエンジンへ移植するセッションを、このセッションの後に挿入する（README の表に行を足す）。TS で計算してはいけない。
2. Run / Query の選択、Subject 一覧、HSP の一覧。表の値は outfmt 6 の該当行（HSP レコードの `out6` の範囲）から取り、並べ替えにはエンジンの原値を使う。選択の中心は HSP の ID（run、`q_idx`、`rank`）にする。
3. 詳細：選んだ HSP の subject 見出しと、outfmt 0 の該当節（`out0` の範囲）を原文のまま表示する。`out0` が null の HSP（outfmt 0 に現れないもの）は、その理由を示す。
4. ドットプロット：選んだ query と subject の組を Canvas で描く。ズーム、パン、向きの識別、HSP の選択。軸の単位（nt / aa）と frame を示す。BLASTN の 1 文字の HSP は start と end が同じで、向きが座標から分からない（`docs/evidence/losat_web_e2c/AUTHORITY.md` §P、`docs/web/abi_v2.md` §8）。HSP レコードに鎖の欄を足す（ABI の追加。エンジンの HSP の query の frame）か、outfmt 0 の節の `Strand=` から読むかを決める。
5. Run の詳細：RunSnapshot と RunRecord、CLI コマンド。検証バッジ（計画 §6.1）の判定表を、`docs/web/verification_cells.tsv` と既存の認証記録から生成するスクリプトを作る。手で書かない。
6. 表示用のフィルター（ViewState）。再検索とは別の操作にし、表示 0 件、ヒット無し、失敗、結果数の上限に達した可能性、を区別して示す。上限に達しただけで、取れなかったヒットがあると断定しない。
7. 試験：HSP レコードと outfmt 6 の行・outfmt 0 の節の対応（個数と座標、null の場合を含む）を、対応済みの全 program の fixture で確かめる。対応済みの全 program の E2E（SX がまだなら BLASTX を除き、そのことを記録する）。

完了条件は計画 §7 の S13 の行による。画面の記録を画面レビューに見せ、`docs/evidence/losat_web_w4/README.md` に記録する。

## S12（W3）から引き継ぐこと

[W3 のゲート記録](../evidence/losat_web_w3/README.md)の要点：

- **画面の構成**：`src/ui/AppView.vue` は `coordinator`・`draft`（`SearchDraft`、検索画面の下書き）・`attention`（Wake Lock、離れる前の注意、復帰時の確認）を受け取る。検索画面は `v-show` で保たれ、結果タブを見ている間も次のジョブの編集が残る。結果タブは S01 の最小の `ResultsPanel.vue`（run の選択、形式、CLI コマンド、原文、書き出し）のままで、`result-output` の `data-shown="<run 番号>:<形式>"` が読み終えた表示を示す（E2E はこれを待つ）。W4 はこれを置き換える。
- **run の識別**：RunSnapshot に `group`（`groupId`・`position`・`size`。Separate の検索）が増えた。Run の一覧は番号と program のラベル、入力名、`Options:`（argv の入力より後の語）を示す。結果画面の Run の選択も、同じ情報で区別できるようにする。
- **E2E の補助**：`tests/e2e/search.spec.ts` の `paste`・`openFiles`・`settled`・`showRecords`・`result` は、入力を作って run を完了させるのに使える。エンジンのビルドだけの試験（長い検索、エンジンの文言）は `BUILD_HAS_ENGINE` で分けている。長い検索は BLASTP の `NZ_CP006932.faa` どうし。
- **FakeEngine の `describe`**：`src/infra/fake/describe.json` はエンジンの `describe` の写しで、reactor があるときの単体試験（`tests/unit/engine-runtime.test.ts`）が一致を確かめる。この試験の file snapshot（`toMatchFileSnapshot`）なので、エンジンの option を変えたら `LOSAT_WEB_REACTORS=<dir> npx vitest run tests/unit/engine-runtime.test.ts -u` で書き直す。
- **キューから結果へ**：W3 の画面レビュー（1 回目の L7）は、完了した Run から結果を開く手段がキューに無いことを指摘した（今は結果タブを開いて Run を選ぶ）。W4 で、Run の一覧から結果を開けるようにする。2 回目の L-e（実行中の Run が完了した Run の下にあり、狭い画面では遠い。待ちと取消の印が同じ灰色）も、Run の一覧を変えるときに扱う。
- **Run の入力**：`Coordinator` は Run の入力のバイト（`runInputs`）を解放しない。Run を消す操作を入れるなら、そこで解放する（W3 のゲート記録の「残件」）。
- **多数のレコード**：W3 の計測では、10 万レコード（34 MB）の入力も索引・検査・積むことができた（索引は 3〜10 秒）。結果画面の Query の選択と一覧は、それだけの query の Run でも仮想化して扱えるようにする。
- **画面の記録と画面レビュー**：`tests/e2e/screens.spec.ts`（`LOSAT_WEB_SCREENS=<directory>`、エンジン入りのビルド）が、同じ表示サイズ 1280×900 と 390×844 で、3 ブラウザの全ページの PNG を撮る。W4 は結果画面の状態をここに足す。W3 の記録は `/home/kawato/.cache/losat-web-gui-target/app-s12-gate-screens/`（リポジトリに入れない。SHA-256 は W3 の実行記録の `screens.sha256`）。W4 の画面レビューは、検索画面が変わらないことも確かめる。
- **ゲートの script**：W3 の [`run_gate.sh`](../evidence/losat_web_w3/run_gate.sh)（V-ABI quick、check、単体、両方のビルドの E2E、繰り返し、計測、画面の記録）を元にできる。

## 終了・引き継ぎ

README の規則 8 に従う。集約列の移植を挿入した場合は、その指示書を書く。挿入しない場合、次は [S14 — 抽出と候補](session_s14_w5_extraction_candidates.md)。
