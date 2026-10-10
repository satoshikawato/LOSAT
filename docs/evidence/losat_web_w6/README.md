# LOSAT Web W6（Session S15）ゲート記録

- 段階：W6 出力と再現性（[総合計画書](../../losat_web_gui_plan.md) §7 の S15、[指示書](../../losat_web_gui_sessions/session_s15_w6_export_session.md)、[セッションファイルの形式](../../web/session_file.md)、[対応表](../../web/ncbi_ui_mapping.md)の §2「1 ページの結果（classic）」）。完了条件は計画 §7 の S15 の行：セッションを再読込しても再計算しない、元配列が無いとできない操作を明示する、勝手につなぎ直さない
- ブランチ：`feature/losat-web-gui-app`（アプリ側。worktree `$WORK_ROOT/.worktrees/web-gui-app`、Linux の clone）。起点は `306b458d`（W5 のゲート記録と S15 の指示書の引き継ぎ）。エンジン側のブランチ `feature/losat-web-gui` は `141615ec` で、S15 の間に動いたもの（SFe）はこのブランチに取り込んでいない。reactor は `$BUILD_ROOT/s14-reactors`、ネイティブの CLI は `$BUILD_ROOT/s14-gate-native/LOSAT`（sha256 `0236291b…`）で、どちらも `141615ec` の木から作った（W5 のものをそのまま使った。実行記録の `environment.txt` に wasm のハッシュ）。エンジン（`LOSAT/`、`web/adapter/`）、計画、README の表、依存ライブラリ（`package.json`）は変えていない。`docs/web/` は `session_file.md`（新規）と `ncbi_ui_mapping.md`（1 ページの結果、ドットプロットの SVG など）だけを変えた
- 実行記録：[`run-20261010T154425Z/`](run-20261010T154425Z/)（ゲート 1、`da3a3565`。結果画面の E2E の繰り返しで 1 件落ちて止まった）、[`run-20261010T161642Z/`](run-20261010T161642Z/)（ゲート 2、`f25c0d69`。画面の記録の WebKit の電話 1 件を除いて通った）、[`run-20261010T170931Z/`](run-20261010T170931Z/)（ゲート 3、`08e4557b`。すべての段階が通った）。どれも作成後は書き換えない。再現は [`run_gate.sh`](run_gate.sh)、NCBI との比較は [`check_commands.py`](check_commands.py)。ファイルのハッシュは [`evidence.sha256`](evidence.sha256)。画面の記録（PNG）はリポジトリに入れない（下の「残件と注意」）
- 判定：**完了条件を満たした**。最後の木（`08e4557b`）のゲート 3 はすべての段階が通り（計測もこの木で測り直した）、画面レビュー 3 回目はその記録で合格（High・Medium・Low 無し、参考 4 件）。コードレビューは 2 回（1 回目：High 無し、Medium 2 件（M1・M2）と Low 5 件（L1〜L5）をすべて直した。2 回目：High・Medium 無し、Low 3 件（L1〜L3）を直した）。画面レビュー 1・2 回目は不合格（どちらも M1）で、直して撮り直した。ゲート 1・2 は試験の側の 1 件ずつで止まり、ゲート 2 は計測の query の場面の記録も空だった（計測の spec を `08e4557b` で直した。下の「ゲートの実行」）
- 保守者の判断待ち：判断 1・2・12・13 は保守者の判断（2026-10-10。1・2 は始めに質問して答えを得た。12・13 は保守者の依頼）で、判断待ちではない。この 4 件の DW を計画に書くのはエンジン側（アプリ側は計画と README の表を変えない。下の「合流」）。それ以外の判断 3〜11・14〜17 は推奨案で進めた（Owner-delegated、2026-09-29 の常設の指示、2026-10-07 に再掲）。3〜11・14 は 2026-10-10、15〜17 は 2026-10-11（JST）に取った。保守者が予約した判断は含まない
- `feature/losat-web-gui` への merge：このセッションでは行わない。エンジン側が行う（下の「合流」）
- このセッションは取りまとめ役（Opus）と作業単位 WP-0・A・B・C・D・E・G・H の agent、修正 4 回（fix 1〜4）で行った。各 agent は一時的な worktree と作業用の branch で作業し、取りまとめ役が `feature/losat-web-gui-app` に merge した。コードレビューと画面レビューは別の agent（`losat-reviewer`）。ゲートの script の実行は取りまとめ役が背景で行った

## コミット

`306b458d..8fb986dc` の 60 コミット（merge 8 を含む）と、この記録のコミット 2 つ（古い順）。`7f6015e9`・`be77ba07` は agent の作業用 branch にアプリ側の先頭を取り込んだ merge。

| コミット | 内容 |
|---|---|
| `a3fefe00` | アプリ：書き出しの Writer の契約（設計書 §12.1）：`ExportSink`、ブラウザの Blob の sink、`ExportWriter`、`writeFile` |
| `7625890b` | アプリ：W5 の残件（L-A の試験、L-B の loop、L-a のトレイ ID の折れ）と、reader の変更の通知（保守者の判断 2） |
| `270e9845` | アプリ：互換出力、抽出、整列行を Writer で区切って書く（指示書 1、W5 コードレビュー L1） |
| `f4d99ccb` | アプリ：ドットプロットの SVG（Download SVG）：表示のとおり、Writer で書き、エスケープし、script も参照も無い |
| `2ecf06f0` | アプリ：`Downloader.save` を外した（すべての書き出しが `open` と Writer を通る） |
| `bdad94b4` | アプリ：LOSAT Web 独自の HSP の file：CSV、JSON、静的な HTML のレポート（指示書 2） |
| `5eb6d491` | アプリ：設定ファイル、Edit Search、Run の再現コマンド（指示書 3・5） |
| `4ab5c781` | merge：WP-A（既存の書き出しを Writer の契約に） |
| `36a1e5bf` | アプリ：検索フォームの設定の行、Edit Search、Run details の「Reproduce this run」 |
| `267aebeb` | アプリ：レポートは行末で始まる文字列を `<pre>` の中でそのまま保つ |
| `97e4fb59` | アプリ：Outputs の「LOSAT Web formats」：範囲ごとの CSV・JSON・レポート |
| `d01df9a3` | 文書：`check_commands.py`。固定した argv の再現コマンドを LOSAT と NCBI BLAST+ 2.17.0 で動かして比べる（指示書 5、設計書 §12.3） |
| `6ede10c9` | merge：WP-B（CSV・JSON・レポート） |
| `71b0267c` | アプリ：セッションファイル（指示書 4）：完了した Run の保存、検索せずに開く、一致するときだけ原 FASTA をつなぎ直す |
| `8492e1e6` | アプリ：セッションファイルの形式の文書、その E2E、約 1 MiB のブロックでの gzip の出力 |
| `14526fbe` | アプリ：Edit Search が埋められない理由を言う。比較の固定ケースは LOSAT が検索できるオプションを使う |
| `f670f5b7` | merge：WP-D（セッションファイル） |
| `50feb0aa` | merge：WP-C（設定ファイル、Edit Search、再現のパネル、Run の入力 FASTA） |
| `d937042e` | アプリ（試験だけ）：セッションのエンジン入りの E2E：実検索を保存し、検索せずに開く |
| `a8b10de4` | アプリ：セッションから読んだ Run の再現のパネルは、入力の出所を言う |
| `fc7ac1d7` | アプリ（試験だけ）：W6 の画面の記録 29〜40（Outputs、再現、設定、セッション、つなぎ直し） |
| `5cfb5d6b` | アプリ（試験だけ）：計測に Run の file とセッション file を足す（W5 の規模） |
| `736a121f` | アプリ（試験だけ）：W6 のゲートの script（W5 の段階に、W6 の単体試験・E2E・記録と NCBI の比較コマンドを足したもの） |
| `c972b648` | アプリ：読み込んだ Run の HSP レコードを型付き配列に入れる前に JSON として検査し、読み込みを atomic に（コードレビュー M1、L5） |
| `e73bb74a` | アプリ：読み込んだ Run の検証バッジは、このサイトのエンジンが書いた出力にだけ付ける（コードレビュー M2） |
| `7f0f7bdf` | アプリ：読み込んだ Run の入力 FASTA は、組み直した入力の SHA-256 を検査し、記録された file の出所を言う（コードレビュー L1、L2） |
| `b9ef6f97` | アプリ：つないだ入力は、選んだ file の順序に依らずつなぎ直せる。拒否した試みと置き換えた試みの保持物を解放（コードレビュー L3） |
| `189c9890` | アプリ：Run を保存できなくなったら読み込みをすぐ止める（コードレビュー L4） |
| `a9e2dd00` | アプリ：読み込んだ Run のバッジは、保存した LOSAT Web と出力を書いたエンジンの build を分けて示す（コードレビュー M2） |
| `fa81f4e5` | merge：修正 1 回目（コードレビュー 1 回目の M1、M2、L1〜L5） |
| `c62324b3` | アプリ：読み込んだ Run の HSP レコードは null の score を許す（アダプタは有限でない値を null で書く） |
| `05f7ad62` | アプリ（試験だけ）：計測は最長のフレームの間隔の終わりとタイマーの遅れを記録する。読み込んだ Run の入力の文が記録どおりであることを確かめる |
| `7f6015e9` | merge：修正 1 回目をゲートの作業用 branch に取り込み |
| `4495f852` | アプリ（試験だけ）：修正 1 回目の後のエンジン入りのセッションの試験（バッジの行、別の順序でのつなぎ直し） |
| `001b8db7` | アプリ：HSP レコードの batch を、1 件ごとの OPFS の読みではなく少数の範囲の読みで読む |
| `af8fb81f` | アプリ：レポートは自分の部分と outfmt 0 の窓だけを待ち、整列を除くことができる |
| `551803b4` | アプリ：画面レビュー 1 回目の修正（reader の通知はそのプログラムで終わる（M1）、Run details の保存ボタンとコマンドの注、部分的に適用した設定の警告、句読点、関係の文言、radio、オプションの数、Session の候補とメモの数） |
| `0415cb55` | アプリ：大きな書き出しの間、ページが描画を続ける：Writer の block の再利用、JSON の group、Outputs の text の contain |
| `948a75e3` | アプリ：Run details の保存ボタンは、原 FASTA をつなぐとすぐ追随する（パネルが生きた Run を読む）。落としていた E2E を足した |
| `034b9e9a` | アプリ：10 万 query のセッションを、ページが描画を続ける区切りで保存し、開く |
| `d09638cd` | アプリ（試験だけ）：計測のタイマーの遅れは、タイマーを置いた時点から取る |
| `157488f3` | アプリ：フレームが描かれなかったときはフレームを待つ pause、書き出しは時計でも pause |
| `9868b008` | アプリ：結果の最初のタブは NCBI の classic の 1 ページ（保守者の依頼、判断 12）：Graphic Summary、その直下の Descriptions、選んだ subject の Alignments を NCBI と同じように切る |
| `44440cd3` | アプリ（試験だけ）：Ranges の確認は節ではなくブロックを scroll する |
| `ec29d263` | merge：修正 2 回目（書き出しとセッションの大きな規模での速さ） |
| `be77ba07` | merge：修正 2 回目を WP-H の作業用 branch に取り込み |
| `7559817a` | アプリ：1 ページの Descriptions と Graphic Summary の「Show all」は、結果のタブを離れて戻っても、別の Run か query になるまで残る |
| `da55ab67` | アプリ（試験だけ）：結果の計測は、切られているときだけ 200 subject の query の Descriptions を全部出す |
| `4fbace1d` | アプリ：対応表に、1 ページの計測（W5 との比較）と「Show all」が残ることを記録 |
| `df35f774` | アプリ：切った Descriptions の見出しは「100 of 260 listed」、「select all listed」は表示中の行だけ、印は表示中の subject に従い、デスクトップは約 22 行（画面レビュー 2 回目 M1・I1・I2、コードレビュー 2 回目 L2） |
| `84c79848` | アプリ：読み込んだ Run の HSP レコードの rank は query ごとに 0〜n−1、候補の rank は Run のレコードの範囲内（コードレビュー 2 回目 L1） |
| `fbb48a9d` | アプリ：入力の関係の文言（除いた record が 1 つのとき）（画面レビュー 2 回目 I5） |
| `a00e4a7e` | アプリ：電話の Filter Results は閉じた disclosure、設定した数を表示（画面レビュー 2 回目 L1） |
| `da3a3565` | アプリ（試験だけ）：Descriptions の見出しと印、デスクトップの高さ、電話の disclosure の E2E。行はマウスで押す（画面レビュー 2 回目 L2 は Playwright の click だけで起きた）。**ゲート 1 の木** |
| `d3342a0e` | アプリ（試験だけ）：1 ページの試験は、並べ替えやグラフの矢印キーの後に 2 px の動きを許す（ゲート 1 の Firefox の 873 対 874）。**ゲート 2 の木のアプリ側** |
| `f25c0d69` | 文書：ゲートの最初の実行（`da3a3565`）の記録。**ゲート 2 の木** |
| `621f23aa` | アプリ（試験だけ）：画面の記録の spec は、Run の結果を開くとき「Open results」へ click を点ではなく event で送る（判断 17。アプリは変えていない） |
| `ea4a35bf` | 文書：ゲートの 2 回目の実行（`f25c0d69`）の記録（画面の記録の WebKit の電話で止まった。計測の query の場面は記録が空） |
| `08e4557b` | アプリ（試験だけ）：計測は Descriptions の見出しを「N shown」と「L of N listed」の両方で読み、query の場面に記録された失敗で試験が落ちる（ゲート 2 の計測の記録が空だった）。**最後のアプリの木、ゲート 3 の木** |
| `8fb986dc` | 文書：ゲートの 3 回目の実行（`08e4557b`）の記録（すべての段階が通った） |
| このゲート記録を含むコミット | 文書：W6 のゲート記録と `evidence.sha256` |
| S16 への引き継ぎのコミット | 文書：S16 の指示書に「S15（W6）から引き継ぐこと」 |

最後のアプリの木は `08e4557b`（アプリの本体は修正 4 回目の `a00e4a7e` から変わっていない。後は試験と記録だけ）。エンジン（`LOSAT/`、`web/adapter/`）、計画、README の表、セッションの README、`package.json` は、アプリ側のブランチでは変えていない（README 規則 1）。

## 完了条件と結果

指示書 1.〜6.（と計画 §7 の S15 の行）の各項。試験の名前は `web/app/tests/` のもの。

| 条件 | 結果 | 証拠 |
|---|---|---|
| 1. 互換出力（outfmt 0/6/7）を 1 バイトも変えずに出し、範囲は Run の結果全体だけ。書き出しは全体をメモリに組み立てず、ブロックを順に書く | 通過。保存したテキストをそのまま Writer に流す。抽出と整列行の書き出しも Writer を通る（W5 コードレビュー L1。抽出は 1 回の読みを 4,194,300 残基以下にし、全体を読むレコードは `checkRecord` で長さと残基の数を照らす。`Downloader.save` は無くなった）。エンジン入りの E2E は outfmt 0/6/7 が保存前・保存後・ネイティブの CLI で一致することを見る | `export-writer.test.ts`（5 件。ブロック、UTF-8 の byte 数、失敗したら何も保存しない）、`extraction-steps.test.ts`（5 件）、`extraction-engine.test.ts`、`exports.spec.ts` のエンジン入り「the files of a BLASTN search hold its stored outputs as written」、`session.spec.ts` のエンジン入り（下の 6.） |
| 2. CSV・JSON（HSP レコードから作る。ViewState のフィルターや選択を反映）、静的な HTML のレポート。LOSAT Web の形式と明記し、NCBI と一致するというラベルを付けない | 通過。範囲は Run 全体・フィルター後・印を付けた subject の 3 つ。CSV は RFC 4180（CRLF、BOM 無し、ID はそのまま）、JSON は HSP レコードの値を書き直したもの（値は同じ）、レポートは外部の script も参照も計測も無く（CSP `default-src 'none'; style-src 'unsafe-inline'; base-uri 'none'; form-action 'none'`）、入力のテキストはすべてエスケープする。画面に「LOSAT Web formats」と書き、NCBI の形式ではないと言う。10 万 query（149,828 HSP）で JSON 1.8〜2.9 s、レポートは整列を含めても約 2 s（下の「実測」） | `hsp-export.test.ts`（14 件）、`result-export.test.ts`（21 件）、`exports.spec.ts` の FakeEngine「CSV, JSON and the report of each scope say what they are and keep hostile IDs as text」、判断 6・14 |
| 3. 設定ファイル：program、入力名を除いた argv、スレッドだけを書き出し、読み込むとフォームに入る | 通過。入力は書かず、読み込んでも検索は始まらない。未知の field、新しい schema、BLASTX、壊れた file は拒否し、一部しか入らないときは入らなかったものを警告する。Edit Search は設定と Job Title を入れる（入力は入れない） | `settings-file.test.ts`（12 件）、`run-files.test.ts`（15 件）、`repro.spec.ts`「settings: saved from the form, the form changed, loaded back; broken and newer files are refused」、同「Run details: the commands from the run's argv…; Edit Search fills the form without searching」 |
| 4. セッションファイル：manifest（schema の版、app / engine の版、RunSnapshot、レコード表、入力の SHA-256）と、長さ付きのブロックの container を `CompressionStream` で gzip。読み込んでも検索しない。不一致・欠け・長すぎる値・不正な参照を拒否。入力の名前や内容を HTML・script・URL として解釈しない。OPFS のパスとトークンは書かない。flank・全長の抽出は原 FASTA を選び直し、SHA-256 が一致したときだけつなぐ | 通過。container 1・schema 1（[形式](../../web/session_file.md)）。候補とメモを含めるか（保守者の判断 1）は「Include candidates and notes」で、既定は入れる。拒否は atomic（何も残さない）で、HSP レコードは型付き配列に入れる前に JSON として検査する（コードレビュー M1）。どの位置で切れた file も拒否する。つなぎ直しは、選び直した file をレコードごとに照合し、組み直した入力の SHA-256 が manifest と同じときだけ | `session-file.test.ts`（24 件）、`session-data.test.ts`（10 件）、`session-save-load.test.ts`（12 件。「writes no storage path, token, run, revision or source ID of the session that saved it」）、`session.spec.ts`（FakeEngine 4 件：「a damaged or foreign file is refused with its reason and changes nothing」「names, IDs, titles and notes that look like HTML, scripts or URLs are shown as text, and run nothing」を含む。エンジン入り 1 件） |
| 5. 再現の手順：LOSAT のコマンドと NCBI との比較用コマンドを argv から作る。既定以外の subject の遺伝暗号の TBLASTX / TBLASTN は承認済みの例外と注記。除外・つないだ Run は「この実行に使った入力 FASTA」を書き出せる | 通過。NCBI のコマンドは、すべてのオプションが NCBI 2.17.0 の `-help` で確かめた表にあるときだけ出し、NCBI の CLI が拒む遺伝暗号（gencode 32）には理由を出す。比較の検査は下の「コマンドの比較」。入力 FASTA は検索したバイトのまま、argv の名前で保存する | `reproduce.test.ts`（9 件）、`run-files.test.ts`、`repro.spec.ts`「a run that NCBI BLAST+ cannot run as LOSAT did gets the reason instead of an NCBI command」、`check_commands.py`（下） |
| 6a. セッションを再読込しても再計算しない | 通過。エンジン入りの E2E は、開く間にエンジンと thread の worker が起動しないこと、Run の状態が queued・running にならないことを見て、出力が保存前と同じ（outfmt 0/6/7、CSV、JSON の HSP、レポート）で、ネイティブの CLI のものとも一致することを確かめる | `session.spec.ts` のエンジン入り「real searches saved and opened in a new page: nothing is searched, the outputs are those before saving and the native CLI's, and the originals attach only when chosen and matching」、`session-save-load.test.ts`「saves two runs with candidates and notes, and opens them in a fresh app without searching」 |
| 6b. 元配列が無いとできない操作を明示する | 通過。読み込んだだけの Run は、配列の抽出と「この実行の入力 FASTA」の保存が「Needs the original … FASTA」で無効になり、整列行の書き出しは保存した出力だけで足りるので使える | `session-save-load.test.ts`「refuses extraction until the original is chosen again and matches; the alignment export needs none」、`session.spec.ts`「Run details attaches the original FASTA only when it matches; then sequences are extracted as before saving」、画面の記録 34・35・37 |
| 6c. 勝手につなぎ直さない | 通過。つなぎ直しは利用者が選んだ file だけ。セッションから読んだ Run にしか付かない。選んだ file が違えば（レコードの SHA-256、数、順序）、何が違うかを名指して拒否する | `session-save-load.test.ts`「attaches only to a run loaded from a session file, and only by an explicit choice」、`session-file.test.ts`「names the first record of chosen files that differs from the saved input, or the count」「matches chosen files to the recorded sources by their records, whatever order they were chosen in」 |
| 6d. 壊れた session file を拒否する | 通過。gzip でない、壊れた、切れた、新しい schema、欠けたブロック、manifest の上限超過、ブロックの余り、存在しない Run・HSP・record を指す候補、ranks が 0〜n−1 でない HSP レコード | `session-save-load.test.ts`（「refuses a file that is not gzip, damaged gzip data, and a cut file」「refuses a newer schema, a missing block and a manifest over the limits, and shows hostile text only as data」「refuses HSP records and candidates that the runs cannot have」ほか）、`session-file.test.ts`（「refuses a container cut at any byte: at every block boundary and inside every block」ほか） |
| 保守者の依頼（判断 12・13）：結果は NCBI の classic のように、タブを開かずに Graphic Summary、直下の Descriptions、Alignments が 1 ページに出る。NCBI のように切る | 通過。最初のタブ「Descriptions」が 1 ページ。Graphic Summary と Descriptions は先頭の 100 subject（「Show all N」）、Alignments は選んだ subject だけ。デスクトップでは「Open results」の直後にグラフの棒が見える。10 万 query の Run を開くのは 276〜306 ms（W5 は 281 ms） | `results.spec.ts`「the one page (NCBI classic): the Graphic Summary, the Descriptions and the Alignments without a tab click; a Description chosen brings its alignments into view」、画面の記録 07・09・14・18〜21・24・25・33（画面レビュー 2 回目の Check 2）、`logs/wp-h/report.md` の計測 |
| ドットプロットの SVG（W4b の判断 24） | 通過。表示のとおり（拡大・移動・フィルター後の HSP）、選択の halo・hover・凡例は無く、script も参照も無い | `plot-svg.test.ts`（19 件）、`results.spec.ts`「the dot plot: Download SVG saves the plot as shown, as a plain SVG with the canvas's titles and without the halo (S15)」 |
| W5 からの残件（L-A、L-B、L-a、I2） | 通過。トレイの欄の幅は loop、電話のトレイの Subject ID は `_ . \| :` の後でだけ折れる、プログラムを替えて reader の kind が変わると、その役割を索引し直して除外を消し、検索画面がそう言う（保守者の判断 2） | `search.spec.ts`「a program of another reader kind reads the sources again, clears their exclusions and says so (Owner decision 2)」「the reader-change notice never outlives its program (screen review M1)」、画面の記録 04・05・27・28 |
| 速さ（指示書「S15 の書き出しも、この規模で同じ操作の速さを保つ」） | 通過（WebKit の書き出しの開始の 128〜155 ms を除く。下の「実測」と「残件と注意」）。10 万 query（149,828 HSP）の書き出しとセッションの保存・読み込みは 3 ブラウザで 0.5〜3.5 s、最長のフレームの間隔は 100 ms 以下 | `logs/fix2/report.md`（タスクフォルダ）、`run-20261010T161642Z/measure/` |
| 3 ブラウザ・両方のビルドで既存と新しい E2E が通る | ゲート 1・2 とも FakeEngine のビルド 146 件（Chromium 53・Firefox 47・WebKit 46）、エンジン入りのビルド 170 件（61・55・54。V-BR を含む）。W5 の最後の木の 110 件、140 件から、`exports.spec.ts`・`repro.spec.ts`・`session.spec.ts`、1 ページの結果の試験、Descriptions の見出しと電話の disclosure の試験などが増えた。ゲート 3 は 同じ件数で通った（FakeEngine 146 件、エンジン入り 170 件、繰り返し 90 件） | `run-20261010T154425Z/`、`run-20261010T161642Z/`、[`run-20261010T170931Z/`](run-20261010T170931Z/) の `npm-e2e-fake-engine.log`・`npm-e2e-engine.log` |
| HSP の対応の試験が通る（変えずに） | 通過。`hsp-correspondence.test.ts` は W4 から変えていない。reactor 付きの単体試験は 38 ファイル、716 件通過・1 件 skip（W5 は 550 件）。V-ABI quick：15 の検索 × 経路とスレッド = 60 の実行がネイティブの CLI と一致、NCBI の凍結ハッシュ 16/16 一致 | 各実行の `unit-cases.log`・`v-abi.log` |

### コマンドの比較（`check_commands.py`）

NCBI BLAST+ 2.17.0（`~/micromamba/bin`）と LOSAT のネイティブの CLI に、再現のパネルが作る固定ケースのコマンドを与えて出力を byte で比べた。ゲートでは oracle の lock（`oracle.lock`）の下で動かす。42 の比較のうち、33 が一致、6 が承認済みの例外（既定以外の subject の遺伝暗号 `-db_gencode` の TBLASTX / TBLASTN。PD-TLOSAN-LOCAL-GENCODE-32 と同じ扱い）、3 は NCBI が拒むケース（gencode 32）。ゲート 1・2 とも同じ。結果は各実行の `check-commands.tsv`・`.json`・`.log`。

### ゲートの実行

**ゲート 1（`da3a3565`）**：2026-10-11 00:44〜01:15（JST）、約 31 分で止まった。V-ABI quick（60 実行が一致、16/16）、`npm run check`（38 ファイルのうち 36 が通り 2 が skip、656 件通過・61 件 skip）、reactor 付きの単体試験（38 ファイル、716 件通過・1 件 skip）、`check_commands.py`（42 の比較）、FakeEngine のビルドの E2E（146 件）、エンジン入りのビルドの E2E（170 件）が通った。結果画面の E2E の 2 回の繰り返し（90 件）の 1 件が Firefox で落ちた：1 ページの試験（`results.spec.ts:1097`）で、並べ替えの後の `scrollY` が 873 対 874 だった。ページの動きの 1 px の揺れで、試験が 2 px を許さなかっただけ。`d3342a0e` で 2 px の許容にし、直した試験を Chromium・Firefox・WebKit で 10 回繰り返して通した。ここで script が止まったので、計測と画面の記録はこの実行に無い。

**ゲート 2（`f25c0d69`、アプリ側は `d3342a0e`）**：2026-10-11 01:16〜01:54（JST）、約 37 分。V-ABI quick（60 実行が一致）、`npm run check`（656 件・61 skip）、単体試験（716 件）、`check_commands.py`（42 の比較）、FakeEngine の E2E（146 件）、エンジン入りの E2E（170 件）、結果画面の E2E の繰り返し（90 件、3 ブラウザとも 30 件）、計測（3 ブラウザ、各 3 件が通った扱い。ただし query の場面は下のとおり記録が空）が通った。画面の記録は 17 件が通り、1 件が落ちた：WebKit の電話で、Run 1 の「Open results」の click が Run 2 の行に当たった（Run 3 の Alignments の節が Queue の上で遅れて読まれて内容を押し下げたため。WebKit は scroll anchoring を持たない）。単独で 3/3 通った。画面の記録の spec が click を点ではなく event で送るようにして直した（判断 17。アプリは変えていない）。**計測の query の場面は、この実行では記録が取れていない**：10,000 と 100,000 query の場面が 3 ブラウザとも `object null is not iterable` で記録に `error` を残した（試験は通った扱い）。原因は `results-measure.spec.ts:815` が Descriptions の件数を `/^([\d,]+) shown$/` で読むのに、修正 4 回目（`df35f774`）が切った一覧の見出しを「100 of 260 listed」にしたこと。`08e4557b` で直した：見出しを「N shown」と「L of N listed」の両方の形で読み、query の場面に `error` が残れば試験が落ちる（ドットプロットの大きな組の見込みどおりの失敗は今までどおり記録だけ）。直した spec は Chromium の 10,000 query で通り、Run の file とセッションの項も記録した。そのため 10 万 query の書き出しとセッションの数値は、この実行の記録ではなく修正 2 回目の計測（下の「実測」）による。1,500・3,000 写しの場面（ドットプロット、トレイ、書き出し、セッション）は記録が取れている。

**ゲート 3（`08e4557b`）**：`ea4a35bf` で一度始めたが、V-ABI の途中で止めた（ゲート 2 の計測の記録が空だったことに気づいたため。途中の記録はリポジトリに入れていない）。計測の spec を `08e4557b` で直して最初からやり直した。2026-10-11 02:09〜03:04（JST）、約 55 分（Firefox の計測が 6.7 分と長い）。**すべての段階が通った**（[`run-20261010T170931Z/`](run-20261010T170931Z/)、`gate run passed`）：V-ABI quick（15 の検索 × 経路とスレッド = 60 の実行がネイティブの CLI と一致、NCBI の凍結ハッシュ 16/16 一致）、`npm run check`（38 ファイルのうち 36 が通り 2 が skip、656 件通過・61 件 skip。reactor の要る分は skip）、reactor 付きの単体試験（38 ファイル、716 件通過・1 件 skip）、`check_commands.py`（42 の比較：一致 33、承認済みの例外 6、NCBI が拒む 3）、FakeEngine のビルドの E2E（146 件）、エンジン入りのビルドの E2E（170 件）、結果画面の E2E の 2 回の繰り返し（90 件）、計測（3 ブラウザ、各 3 件通過。10,000 と 100,000 query の場面は 3 ブラウザとも記録が取れ、`error` は下の見込みどおりの組の失敗だけ）、画面の記録（3 ブラウザ × 2 サイズ、258 枚、状態 01〜40 と 30b、18 件通過。SHA-256 は `screens.sha256`）。落ちた試験、flaky、やり直しは無い。

どの実行の E2E の log にも `FAIL [opfs] storage full …` と `FAIL [run output] storage full …` が出る（FakeEngine に 3 行、エンジン入りに 6 行）。ブラウザの中の契約の確認が自分の記録として印字するもので、それを印字する Playwright の試験は通っている（落ちた試験は無い）。W4b・W5 と同じ。計測の見込みどおりの失敗は W4b・W5 と同じ：4000 写しの自己検索は 3 ブラウザでエンジンのメモリ不足（`memory allocation of 49861 bytes failed`）、3000 写しは WebKit で保存の上限（512 MB）。

## 作業との対応

| 作業 | 実装 | 試験 |
|---|---|---|
| Writer の契約（指示書 1・2、設計書 §12.1、W5 コードレビュー L1） | `application/export-writer.ts`（`ExportSink`、`ExportWriter`、`writeFile`。1 つの block を再利用して `encodeInto`、sink は block を write が終わる前に写す）、ブラウザの Blob の sink と gzip（`infra/browser/compression.ts`、256 KiB の区切りで task を挟む）、`infra/browser/next-task.ts`（約 40 ms ごとに timer、フレームが 50 ms 描かれなければフレームを待つ）。抽出は 1 回の読みを 4,194,300 残基以下、小さな FASTA レコードは約 1 MiB にまとめ、整列は 1,000 レコードずつ、互換出力は 8 MiB の範囲で読む（判断 5） | `export-writer.test.ts`、`extraction-steps.test.ts`、`session-data.test.ts`（`browserCompression`）、`extraction-engine.test.ts`、`data-extraction.test.ts` |
| CSV・JSON・レポート（判断 6、14） | `domain/hsp-export.ts`、`domain/hsp-report.ts`、`application/result-export.ts`。範囲は Run 全体・フィルター後・印を付けた subject。JSON は HSP の行を 64 KiB 以下を読み飛ばし 8 MiB 以下の少数の範囲読みで読む。レポートは「Include the alignments in the report」（既定は入れる） | `hsp-export.test.ts`、`result-export.test.ts`、`exports.spec.ts` |
| 設定ファイルと再現（判断 8） | `domain/settings-file.ts`、`domain/reproduce.ts`、`application/run-files.ts`、`ui/SettingsFileControl.vue`、`ui/ReproducePanel.vue` | `settings-file.test.ts`、`reproduce.test.ts`、`run-files.test.ts`、`repro.spec.ts`、`check_commands.py` |
| セッション（判断 1・7・9・10・11） | `domain/session-file.ts`（container と manifest の検査）、`application/session.ts`（保存、読み込み、つなぎ直し）、`ui/SessionPanel.vue`、`ui/LoadedRunOrigin.vue`、`docs/web/session_file.md`。読み込みは block を Data worker に背圧付きで流す（`stagedBytes`、1 MiB の message、同時に 16 以下） | `session-file.test.ts`、`session-data.test.ts`、`session-save-load.test.ts`、`session.spec.ts` |
| ドットプロットの SVG（判断 4） | `domain/plot-svg.ts`、`domain/plot-layout.ts`（canvas と共有）、`application/plot-export.ts` | `plot-svg.test.ts`、`results.spec.ts` |
| reader の変更の通知（保守者の判断 2、判断 3） | `DraftSource.notice`、`{role}-source-{n}-notice` | `search.spec.ts`、`draft.test.ts`、`run-files.test.ts` |
| 1 ページの結果（判断 12・13・15・16） | `ui/ClassicResults.vue`、`ui/shownWhole.ts`、`SubjectTable.vue`、`GraphicSummary.vue`、`ResultsBrowser.setListLimit`、対応表 §2 | `results.spec.ts`、`candidates.test.ts`、`results.test.ts` |
| 速さ（判断 14） | 範囲読み、レポートの部分、`pre.output { overflow-anchor: none; contain: content }`、区切りのあるセッションの検査 | `results-measure.spec.ts`、`logs/fix2/report.md` |
| 計測と画面の記録 | `results-measure.spec.ts`（Run の file とセッションの項）、`screens.spec.ts`（状態 29〜40、30b）、[`run_gate.sh`](run_gate.sh) | ゲートの実行 |

## 作ったもの

- `src/domain/`：`session-file.ts`、`settings-file.ts`、`reproduce.ts`、`hsp-export.ts`、`hsp-report.ts`、`plot-svg.ts`、`plot-layout.ts`。
- `src/application/`：`export-writer.ts`、`result-export.ts`、`plot-export.ts`、`session.ts`、`run-files.ts`。`candidates.ts`・`results.ts`・`draft.ts`・`coordinator.ts` は書き出し・つなぎ直し・1 ページの選択に合わせて変えた。
- `src/ports/`・`src/infra/`：`ports/compression.ts`、`infra/browser/compression.ts`、`infra/browser/next-task.ts`。Data worker（`infra/data/data-service.ts`）は HSP レコードの検査（`checkHspRecords`）、範囲読み、入力の記述と解放、`stagedBytes` を足した。`ports/download.ts` から `Downloader.save` を外した。
- `src/ui/`：`ClassicResults.vue`、`ExportFormats.vue`、`ReproducePanel.vue`、`SessionPanel.vue`、`SettingsFileControl.vue`、`LoadedRunOrigin.vue`、`shownWhole.ts`。`SubjectTable.vue`、`GraphicSummary.vue`、`ResultFilters.vue`、`CandidatesPanel.vue`、`OutputsView.vue`、`RunDetails.vue`、`styles.css` を変えた。
- 文書：`docs/web/session_file.md`（新規）、`docs/web/ncbi_ui_mapping.md` の追記、[`run_gate.sh`](run_gate.sh)、[`check_commands.py`](check_commands.py)、この記録、S16 の指示書への引き継ぎ。

## 判断

保守者の判断は 1・2・12・13（2026-10-10）、それ以外は推奨案で進めたもの（Owner-delegated、2026-09-29 の常設の指示、2026-10-07 に再掲。3〜11・14 は 2026-10-10、15〜17 は 2026-10-11 JST）。出典は取りまとめ役の `DECISIONS.md`（タスクフォルダ `/home/kawato/losat-baselines/s15-w6-20261010/`）と各担当の記録（`logs/wp-*/progress.md`、`logs/wp-h/report.md`、`logs/fix*/`）。

1. **保守者**：セッションファイルは、候補トレイの候補とメモを「Include candidates and notes」の checkbox（既定で入れる）で含める。含めるのは、Run がセッションに入っている候補だけ。読み込むとトレイとメモが戻る。配列の抽出には原 FASTA の選び直しと照合が要る（計画 §5.8、§10）。
2. **保守者**（W5 の判断 12）：プログラムを替えて役割の reader の kind（1 核酸 / 2 タンパク質）が変わると、その役割を索引し直して除外を消し、検索画面がそう言う（「…the exclusions were cleared because …」）。
3. WP-0：reader の変更の通知は source の文字列（`DraftSource.notice`、`{role}-source-{n}-notice`）で、次の選択の変更か削除で消える。除外が無い source は何も言わない。トレイの Subject ID は `_ . | :` の後でだけ折れる（`<wbr>`、`overflow-wrap: break-word`）。
4. WP-E：ドットプロットの SVG は表示のとおり（拡大・移動・フィルター後の HSP）、選択の halo・hover・凡例は無い（`<desc>` が色の意味を言う）。file 名は `losat-run{N}-dotplot-q{Q}-s{S}.svg`。線は入れ子の `<svg overflow="hidden">` で切り（`url(` も `href` も無い）、XML に入らない文字は置き換える。
5. WP-A：抽出は 1 回の `readResidues` で 4,194,300 残基（60 × 69,905）以下、全体を分けて読んだレコードは `checkRecord` で長さと残基の数を照らす、小さな FASTA レコードは約 1 MiB にまとめる、整列は 1,000 HSP レコードずつ、互換出力は 8 MiB の範囲で読む。`Downloader.save` は無い。最大で、ブラウザの保管場所の Blob にレコード全体の写しが 1 つと、heap に約 4.2 MB の部分が約 3 つ。
6. WP-B：CSV は RFC 4180・CRLF・UTF-8（BOM 無し）・見出しの行・注記無し（画面が形式を言い、表計算が `=` で始まる ID を式と読むことを警告する。ID はそのまま）。JSON の数値は HSP レコードの値を書き直したもの（text は違い得るが値は同じ。outfmt 6 の文字列は整形した値）。JSON は HSP レコードを 1,000 ずつ読む。範囲は開始時に決め、書き出しは同時に 1 つ。レポートの CSP は `default-src 'none'; style-src 'unsafe-inline'; base-uri 'none'; form-action 'none'`、link は無い。file 名に範囲は入れない（JSON とレポートの中に書く）。
7. WP-D：セッションの形式 v1（[`session_file.md`](../../web/session_file.md)）：gzip の `LOSAT-WEB-SESSION 1`、`<name> <len>\n<bytes>\n` のブロック（manifest、Run ごとに out0・out6・out7・hits・diagnostics、含めるなら candidates）、`LOSAT-WEB-SESSION-END`。未知の manifest の field は拒否（新しい field は新しい schema）。保存は読み込みの検査を通す（保存した file は必ず開く）。上限：manifest と candidates 64 MiB、Run 1,000、argv 1,000 語 × 10,000 文字、名前 1,000、ID・題 10,000、メモ 100,000、候補 1,000,000、行 64 byte。参照は 1 始まりの Run の位置。読み込んだ Run は新しい Run ID と次の番号を持つ（`RunView.fromSession`）。復元した候補は追加され、選択され、`addedAt` を保つ。つなぎ直しは記録された source の数だけの file が要る（修正 1 回目から順序は問わず、レコードで照合、名前は違ってもよい）。組み直した Run の入力の SHA-256 が manifest と同じでなければならない。つなぐまで「Download FASTA」は要件を示して無効。
8. WP-C：設定ファイルのスレッドは 1〜16（それ以上はブラウザの Auto として読み、「Not applied」に載せる）。BLASTX の設定ファイルは拒否（SX まで）。schema 1 は未知の field を拒否。NCBI のコマンドは、すべてのオプションが NCBI 2.17.0 の `-help` と照合した program ごとの表にあるときだけ出し、そうでないとき（と NCBI の CLI が拒む遺伝暗号、例 TBLASTN の 32）は理由を言う。既定以外の subject の遺伝暗号には、コマンドと `approvedExceptions` の文を付ける。Edit Search は設定と Job Title（Run に無ければ消す）を入れ、入力は入れず、検索を始めない。Run の入力 FASTA は argv の名前で保存する。
9. 統合：読み込んだ Run の `RunFiles.inputBytes` は `Session.attachedInput`（つないだ原 file から組み直した Run の入力）を使い、無ければ何をするかを言って拒否する。
10. 修正 1 回目（コードレビュー 1 回目）：読み込んだ Run の HSP レコードは Data worker で生の JSON として検査してから使う（`RunStore.checkHspRecords`）。`index` は 0〜count−1 を 1 回ずつ（順序は問わない）、score は有限の数か null（アダプタは有限でない値を null で書く）、追加の field は許す。読み込みは Run とトレイの項目を作ってから `addSessionRuns`（atomic）。読み込んだ Run のバッジは、このサイトと同じエンジンの build なら通常のバッジに「Loaded from a session file (saved by LOSAT Web V, build B)…」を加え、別の build なら level `outside`（「Written by another engine build」）。Run の入力 FASTA は組み直して SHA-256 を manifest と比べる。つないだ入力は順序を問わずつなぎ直し、拒否した試みと置き換えた試みは解放する（`DatasetStore.releaseSources`）。`stagedBytes` は Run が失敗すると拒否するので、読み込みはすぐ止まる。復元する候補は 1,000 ずつ、整列行無しで読む。
11. 修正 3 回目（画面レビュー 1 回目）：reader の通知は、後のどのプログラムの変更でも消す（履歴として残さない）。読み込んだ Run の入力 FASTA の保存ボタンは、原 FASTA をつなぐまで「Needs the original <role> FASTA (choose it above)」で無効（パネルは生きた Run を読む）。コマンドの注は、入力が選んだ file そのままでないとき、Run の保存した入力 FASTA を使うよう言う。部分的に適用した設定は警告する。除外の関係の文言は 1 つ。Session のパネルは含める候補とメモの数を言う。
12. **保守者の依頼**：結果画面は NCBI の classic（Traditional）の結果のように、Graphic Summary、直下の Descriptions、選んだ subject の Alignments を、タブを切り替えずに 1 ページに出す。Dot Plot・Run details・Outputs はタブのまま。S13b のタブの構成のこの部分を置き換える。対応表 §2 に記録した。
13. **保守者の依頼**：1 ページの結果は NCBI のように描くものを切る。Graphic Summary は選んだ query の先頭 100 subject（「Show all N」）、Descriptions は先頭 100 subject（「Show more / all」）、Alignments は選んだ subject だけ。大きな Run の計測は最初の画面の退行の確認で、全部を描く理由にしない。
14. 修正 2 回目：JSON は一括して HSP の行を少数の範囲読み（64 KiB 以下を読み飛ばし、1 回 8 MiB 以下）で読む。レポートに「Include the alignments in the report」（既定は入れる。外すと outfmt 0 を読まず、そう書く）。Writer は 1 つの block を再利用して書く（`encodeInto`、`texts()`）。書き出しは約 20 ms ごとに pause する（message、40 ms の timer、フレームの待ち）。`pre.output` に `overflow-anchor: none; contain: content`。セッションの検査は 20,000 レコードずつ（manifest の byte は同じ）。既知：WebKit は書き出しごとの最初に 128〜155 ms フレームを描かない（WebKit 自身の仕事）。
15. WP-H：結果の最初のタブは NCBI の名前「Descriptions」のまま 1 ページを持つ。Description の行を選ぶと Alignments の見出しが focus とともに見える所に来る。グラフの棒、「Show in results」、ドットプロットの「Show alignment」は HSP の Range に行く。ページが動くのは選択のときだけ（矢印キー、並べ替え、絞り込み、query の変更、印では動かない）。「select all」は表示中の subject だけ。100 を超えた subject を選んでも一覧は広がらない。グラフの band は 1 行 8 px・棒 4 px・30 行。「Show all」は Run か query が替わるまで残る。
16. 修正 4 回目：切った一覧は「100 of N listed」「select all listed」。印は表示中の行に従う（`ResultsBrowser.setListLimit`。subject が一覧から出たら外す）。デスクトップの一覧は 22 行（電話は 10）。電話の Filter Results は閉じた disclosure で「(N set)」を出す。読み込んだ HSP レコードの rank は query ごとに 0〜n−1、候補の rank は Run の範囲内。1 つの record の除外は単数の文言。時間の試験は 5,000 と 40,000 候補を比べる（上限 24）。画面レビュー 2 回目の L2（Firefox の横 scroll）は Playwright の click だけで起きたので、画面の spec は行をマウスで押す。
17. ゲート 2（WebKit の電話の画面の記録）：`screens.spec.ts` の `openResults` は「Open results」への click を点ではなく event（`dispatchEvent('click')`）で送る。アプリは変えない。Alignments の節が Queue の上で遅れて読まれて内容を押し下げるのは WebKit に scroll anchoring が無いため。残件として記録する（節の高さを読む間確保する）。

## 実測

計測は `tests/e2e/results-measure.spec.ts`（`--workers=1`、機械は他の処理と共有）。1 回の暖機と複数回の中央値。W5 と同じ作り方で、書き出しの項（Run の file、セッション）が W6 の追加。

**書き出しとセッション、10 万 query（149,828 HSP）**（`logs/fix2/report.md`。前は WP-G の記録 `logs/wp-g/measure/`（Firefox・WebKit、`5cfb5d6b`）と `logs/wp-g/measure-merged/`（Chromium、`4495f852`）、後は `logs/fix2/measure-after/`（`157488f3`）。ms は click から画面が要約か message を出す（outfmt 6 は download）まで。前 / 後）：

| ブラウザ | file | 前 / 後 bytes | 前 / 後 ms | 最長のフレームの間隔 ms |
|---|---|---|---|---|
| Chromium | CSV | 16,242,367 / 16,242,367 | 545 / 502 | 43 / 25 |
| Chromium | JSON | 127,249,939 / 127,249,939 | 137,895 / 1,808 | 3,033 / 100 |
| Chromium | レポート | 176,468,917 / 176,469,035 | 2,971 / 2,224 | 300 / 17 |
| Chromium | outfmt 6 | 8,966,787 / 8,966,787 | 44 / 21 | 17 / 17 |
| Chromium | セッション保存 | 25,874,344 / 25,874,288 | 2,534 / 2,822 | 100 / 17 |
| Chromium | セッションを開く | | 2,661 / 1,999 | 133 / 67 |
| Firefox | CSV | 16,242,367 / 16,242,367 | 453 / 689 | 50 / 83 |
| Firefox | JSON | 127,249,835 / 127,249,939 | 79,190 / 2,882 | 83 / 34 |
| Firefox | レポート | 176,468,917 / 176,469,035 | 25,211 / 2,537 | 534 / 34 |
| Firefox | outfmt 6 | 8,966,787 / 8,966,787 | 118 / 60 | 50 / 17 |
| Firefox | セッション保存 | 26,522,323 / 26,522,888 | 2,417 / 2,752 | 450 / 50 |
| Firefox | セッションを開く | | 4,636 / 3,066 | 51 / 50 |
| WebKit | CSV | 16,242,367 / 16,242,367 | 998 / 725 | 176 / 155 |
| WebKit | JSON | 127,249,835 / 127,249,939 | 3,205 / 1,901 | 137 / 132 |
| WebKit | レポート | 176,468,917 / 176,469,035 | 2,560 / 1,869 | 167 / 128 |
| WebKit | outfmt 6 | 8,966,787 / 8,966,787 | 23 / 24 | 16 / 16 |
| WebKit | セッション保存 | 27,485,522 / 27,485,522 | 2,605 / 3,528 | 517 / 52 |
| WebKit | セッションを開く | | 2,773 / 1,830 | 486 / 85 |

- 目標に対して：CSV 0.5〜0.7 s、JSON 1.8〜2.9 s（前は 80〜138 s）、レポート 1.9〜2.5 s（前は 2.6〜25 s）、セッションの保存 2.8〜3.5 s、開く 1.8〜3.1 s。最長のフレームの間隔は 100 ms 以下：Chromium（JSON の 100 ms は約 107 ms の main thread の task で、probe では garbage collection の停止）と Firefox（34〜83 ms）で満たし、WebKit はセッション（52・85 ms、前は 517・486）で満たした。**満たしていない**：WebKit は書き出しごとの最初に 128〜155 ms フレームを描かない（CSV 155、JSON 132、レポート 128）。その間ページの timer は 20〜41 ms で動き、フレームを待つ pause でも縮まず、書き出しのボタンを無効にするだけで 47〜63 ms、Blob の作成では 0 ms で、1 万 query や 1 組の場面でも同じ（WP-G で 107〜138 ms）。WebKit 自身の仕事で、ページの script の外（この環境に WebKit の profiler が無く、これ以上は追っていない）。
- JSON の遅さの原因は、`readHspRecords` が 1 件ごとに OPFS を読んでいたこと（Chromium の 148 s のうち 134 s）と、main thread の garbage collection（5.65 s の停止）。レポートは Firefox で await が約 200 万回（1 行・見出し・節ごと）。直し方は判断 14 と `logs/fix2/report.md`。セッションファイルのバイトは同じ（manifest の JSON は `4495f852` の `checkManifest` の出力と一致、セッションの単体試験は変えずに通る）。
- ゲート 2 の 1,500・3,000 写しの場面（`run-20261010T161642Z/measure/`）：Chromium の 3,000 写し（5,993 HSP）は JSON 363,539,407 byte・2,244 ms、レポート 720,989,181 byte・2,858 ms、セッションの保存 87,496,876 byte・10,901 ms・最長の間隔 17 ms、開く 9,181 ms・50 ms。Firefox の 3,000 写しは保存 86,703,145 byte・10,236 ms、開く 25,531 ms・最長の間隔 166 ms（WP-G の記録の 14,805 ms・133 ms より遅い。単一の実行で、追っていない）。WebKit の 1,500 写しは保存 18,448,644 byte・2,506 ms、開く 2,027 ms。書き出しの最長のフレームの間隔は Chromium・Firefox で 17〜26 ms、WebKit の 1,500 写しで 125〜150 ms（上の WebKit の開始の間隔）。
- **WP-H（1 ページの結果）**（`logs/wp-h/report.md`。Chromium、中央値、`da55ab67` の木 C。W5 の括弧）：10 万 query（149,828 HSP）の Run を開くのは 276〜306 ms（W5 は 281 ms）、最初の HSP の詳細まで 278〜301 ms（261 ms。グラフと Descriptions を Alignments と一緒に mount するので +17〜+40 ms）。選択・絞り込み・並べ替えは 1 フレーム（18 ms 以下）で、1 組のドットプロットのクリックは 29 ms（27 ms）。5,993 HSP の 1 組を開くのは 625〜678 ms（634 ms）。トレイの追加・タブを開く・scroll・並べ替え・抽出 100・全部を外すは 11/12/10/13/14/45/11 ms（W5 は 12/12/11/13/14/45/12）。
- **計測の木**：10 万 query の書き出しとセッション（上の表）は `157488f3`（修正 2 回目）、1 ページの Run を開く時間は `da55ab67`。修正 3・4 回目の後の木では、query の場面を測り直せていない（ゲート 2 では spec が読み取れず記録が空）。ゲート 3（`08e4557b`、最後の木）で 3 ブラウザとも測り直した（[`run-20261010T170931Z/measure/`](run-20261010T170931Z/measure/)。下の表）。

**ゲート 3 の計測、10 万 query（149,828 HSP）**（中央値。ms は ゲート 3 / 修正 2 回目、最長のフレームの間隔も同じ。W5 は S15 の指示書の値）：

| ブラウザ | Run を開く / 最初の HSP の詳細 | CSV | JSON | レポート | outfmt 6 | セッション保存 | セッションを開く |
|---|---|---|---|---|---|---|---|
| Chromium | 258 / 259（W5 281 / 261、WP-H 276〜306 / 278〜301） | 409 / 502（間隔 28 / 25） | 1,794 / 1,808（117 / 100） | 1,374 / 2,224（17 / 17） | 22 / 21 | 2,302 / 2,822（17 / 17） | 1,680 / 1,999（17 / 67） |
| Firefox | 554 / 473（W5 487 / 494） | 372 / 689（67 / 83） | 1,751 / 2,882（17 / 34） | 1,458 / 2,537（17 / 34） | 33 / 60 | 1,867 / 2,752（34 / 50） | 1,944 / 3,066（34 / 50） |
| WebKit | 410 / 461（W5 444 / 450） | 710 / 725（134 / 155） | 1,759 / 1,901（120 / 132） | 1,419 / 1,869（121 / 128） | 20 / 24 | 2,064 / 3,528（48 / 52） | 1,443 / 1,830（48 / 85） |

- 書き出しとセッションは、修正 2 回目と同じか速い。バイト数も同じ（CSV 16,242,367、JSON 127,249,939、レポート 176,469,035、outfmt 6 8,966,787。セッションは保存の時刻の分だけ違う）。Run を開くのは Chromium と WebKit で W5 より速く、Firefox は W5 の 487 ms に対して 554 ms（1 ページに Graphic Summary と Descriptions を一緒に mount する分。最初の HSP の詳細は 473 ms で W5 の 494 ms より速い）。
- 100 ms を超えたフレームの間隔：Chromium の JSON の 117 ms（修正 2 回目は 100 ms。garbage collection の停止）、WebKit の書き出しの開始の 120〜134 ms（既知。修正 2 回目の 128〜155 ms より短い）。
- 1 組の場面：1,500・3,000 写しは Chromium と Firefox で通り、WebKit は 1,500 写しだけ（3,000 は保存の上限、4,000 は 3 ブラウザでエンジンのメモリ不足。見込みどおり）。3,000 写し（5,993 HSP）のドットプロットを出す・拡大・縮小は Chromium 21 / 9 / 4 ms、Firefox 37 / 27 / 27 ms。3,000 写しのセッション（約 87 MB）は Chromium で保存 8,218 ms（間隔 17 ms）・開く 7,686 ms（33 ms）、Firefox で保存 14,863 ms（**間隔 633 ms**）・開く 16,265 ms（100 ms）。Firefox の 633 ms はページの task の遅れが 54 ms で、ページの script の外（ゲート 2 の同じ場面は 50 ms。1 回ずつの値で、追っていない。下の「残件と注意」）。
- トレイ（10 万 query の 200 候補）：Chromium の印・追加・タブを開く・scroll・並べ替え・抽出 100・全部を外すは 11 / 11 / 10 / 13 / 14 / 45 / 12 ms で、WP-H・W5 と同じ。WebKit の抽出は 91 ms。数値の全体はタスクフォルダの `drafts/gate3-measure.md`。

## 独立レビュー（コード）

agent `losat-reviewer`（code review の役）。記録はタスクフォルダの `reviews/code-review.md`、`code-review-2.md`。

**1 回目**（`306b458d..50feb0aa` の `web/app`・`docs/web`・`docs/evidence/losat_web_w6`、73 ファイル、+10,170/−245）：High 無し、Medium 2 件、Low 5 件。修正 1 回目（`c972b648..a9e2dd00`、`c62324b3`）ですべて直した。

| 指摘 | 内容 | 対応 |
|---|---|---|
| M1 | HSP レコードが型付き配列に入れた後でしか検査されない（`"s_idx": null`、`4294967296`、frame 259 などが通る）。候補があると `addSessionRuns` の後に throw し、Run が coordinator に残って保管場所が空になる | `c972b648`：Data worker で生の JSON として検査し、候補を `addSessionRuns` の前に作る（atomic）。null の score は許す（`c62324b3`、アダプタの書き方）（判断 10） |
| M2 | 読み込んだ Run の検証バッジが、出力を書いた build ではなくこのサイトのものを示す | `e73bb74a`・`a9e2dd00`：build が違う・不明なら `outside`、同じなら通常のバッジに「Loaded from a session file…」（判断 10） |
| L1 | 読み込んだ Run の入力 FASTA を全部メモリで組み直し、SHA-256 を再検査しない | `7f0f7bdf`：SHA-256 を manifest と比べる |
| L2 | 読み込んだ Run の「input FASTA」の文が誤り | `7f0f7bdf`：source と除外を file から取る |
| L3 | つないだ入力の再接続が file の順序に依存する | `b9ef6f97`：レコードによる照合（順序を問わない）と、拒否・置き換えた試みの解放 |
| L4 | Run を保存できなくなっても読み込みが続く | `189c9890`：`stagedBytes` が失敗を返して止める |
| L5 | 復元する候補を 1 度に読む | `c972b648`：1,000 ずつ、整列行無し |

**2 回目**（`50feb0aa..4fbace1d`、56 ファイル、+3,406/−732）：High・Medium 無し。修正 1 回目の M1・M2・L1〜L5 は直っていること、修正 2 回目が出力のバイトを変えないこと（`readHitLines` を乱数の 60 Run × 10 batch で、Writer の `encodeInto` を乱数 3,000 件で確かめ、不一致 0）、WP-H が 1 つの選択を保つことを確かめた。Low 3 件を修正 4 回目で直した：

| 指摘 | 内容 | 対応 |
|---|---|---|
| L1 | HSP レコードの `rank` が 0 以上としか検査されず、型付き配列が wrap して 2 つの HSP が同じ識別になり得る（作り込んだ file だけ） | `84c79848`：query ごとに rank が 0〜n−1、候補の rank は Run の範囲内 |
| L2 | 並べ替えや絞り込みで先頭 100 の外へ出た印が、見えないまま数えられ、追加・書き出しに入る | `df35f774`：印は表示中の subject に従う（判断 16） |
| L3 | トレイの線形時間の試験が負荷で落ち得る（比の上限 9 に対し 10.34 を測った） | `df35f774`：5,000 と 40,000 候補を比べ、上限 24（判断 16） |

情報として残したもの（直さない）：null の score は型付き配列で 0 になり、最良の E value として並べ替え・絞り込みに使われる（検索した Run も同じ。アダプタが null を書くのは有限でない値だけ）。貼り付けた入力は session に `query.fa` / `subject.fa` として記録されるので、読み込んだ Run は「query.fa has the same bytes as the file query.fa」と言い、「保存した入力 FASTA を使う」注を加えない（`run-files.ts:195-205` のコメントに書いてある）。隠れた tab では timer が約 1 s に絞られ、書き出しとセッションの保存が遅くなる（待ちには上限がある。修正 2 回目の前より悪くはない）。

確かめて指摘の無かったもの（2 回目から）：生の HSP レコードの検査（型、index、レコード表と byte の範囲、`hitCount` を `Uint8Array(count)` の前に比べること）、atomic な読み込み（遅い失敗でも coordinator・トレイ・選択・保管場所が変わらない）、修正 2 回目のバイト（JSON の group、CSV の `texts`、レポートの 64 Ki の部分、`RangeReader.cached`）、検査の評価順序（最初の拒否も同じ）、WP-H の選択と focus とページの動き、エンジン入りのセッションの E2E が検索しないことを worker とキューの状態で証明すること、ゲートの script（W5 の段階 + `check_commands.py` を oracle lock の下で。記録は `losat_web_w6/run-<UTC>/`、比較の作業は repository の外）。

## 画面レビュー

agent `losat-reviewer`（screen review の役）が、画面の記録を前の記録、対応表、前のレビューと並べて見た（記録はタスクフォルダの `reviews/screen-review.md`、`screen-review-2.md`。行ごとの差の script と切り出しは `tmp/review-screens*/`）。

**1 回目：不合格**（High 無し）。`4495f852` の記録 258 枚（3 ブラウザ × デスクトップ 1280 px・電話 390 px、状態 01〜28 と 29〜40・30b。SHA-256 は 258 枚とも `logs/wp-g/screens-2.sha256` と一致）。W5 の最後の記録（`s14-review-screens`、`608f6f75`）と比べて、状態 01〜28 は S15 が変えるものの他は変わらず、ドットプロットの本体は Chromium でピクセル単位で同一だった。W5 の L-a（電話のトレイの Subject ID）は 3 ブラウザで直っていた。

| 指摘 | 内容 | 対応 |
|---|---|---|
| M1 | 04-tblastn / 04-tblastx / 05：reader の変更の通知が、変えたプログラムを過ぎても残り、今の reader と矛盾する（「Read again as protein for BLASTP…」が TBLASTX の注意の上に出る） | 修正 3 回目（`551803b4`）：どの再読みも通知を置き換え、後のプログラムの変更で消す（判断 11）。2 回目の画面レビューで 3 ブラウザとも直ったことを確認 |
| L1 | 34・37：読み込んだ Run の「Save …fasta」が有効に見えるのに、同じ tab が使えないと言う | `551803b4`・`948a75e3`：つなぐまで無効、要件を示す |
| L2 | 30・34・37：除外した Run のコマンドの注が、保存した file を使うよう言わない | `551803b4`：注を足した |
| L3 | 31：部分的に適用した設定が、ヒントと同じ控えめな style | `551803b4`：notice の style |
| 参考 I1〜I5 | ブラウザの gzip のエラーで句読点が重なる、関係の文言の違い、電話の radio、オプションの数え方、Session の候補とメモの数 | すべて `551803b4` で直した |

**2 回目：不合格**（High 無し）。`4fbace1d` の記録（`s15-wph-screens`、3 ブラウザ × 2 サイズ、状態 01〜40・30b）。1 回目の指摘は 3 ブラウザですべて直っていた。1 ページの結果は Owner の依頼に答える：デスクトップでは「Open results」の直後に Graphic Summary の棒が見え（18）、直下に Descriptions、その下に Alignments。電話では積み重なって横 scroll が無い。ドットプロット（08・15・17・22・23・26）はタブの行を除いてピクセル単位で同一。

| 指摘 | 内容 | 対応 |
|---|---|---|
| M1 | 09・20：Descriptions の見出しが「260 shown」なのに一覧は 100 行（「select all」は表示中の 100 だけを印す） | `df35f774`：「100 of 260 listed」、「select all listed」（判断 16）。ゲートの記録で撮り直す |
| 低 L1 | 電話：「Open results」の後の最初の画面に 1 ページの部分が入らない | `a00e4a7e`：Filter Results を閉じた disclosure に |
| 低 L2 | Firefox の電話：行を選ぶと一覧が横に scroll している | 利用者の操作では起きない（Playwright の `click()` だけ。scrollLeft は `click()` の後 158、マウス・Tab・focus では 0）。画面の spec は行をマウスで押す（`da3a3565`） |
| 参考 I1 | 入れ子の scroll の箱（Descriptions は 10 行） | `df35f774`：デスクトップは 22 行 |
| 参考 I2 | 20：Job Title や通知があると Descriptions が最初の窓の下に入る | コメントを直した（`GraphicSummary.vue`）。受け入れる |
| 参考 I3 | グラフの行が 8 px で、タッチの標的が小さい | 判断 15 で受け入れた（Descriptions の 28 px の行とキーボードがある） |
| 参考 I4 | 電話のタブの行が 3 + 1 に折れる | **直さない**（見た目だけ） |
| 参考 I5 | 1 つの record の除外の文言が複数形 | `fbb48a9d` |

**3 回目**：**合格**（High・Medium・Low 無し、参考 4 件。撮り直しは要らない）。記録は `$BUILD_ROOT/s15-gate-screens/`（ゲート 3 の木。258 枚）、基準は `s15-wph-screens`（`4fbace1d`）と W5 の最後の記録。258 枚とも実行記録の `screens.sha256` と一致し、WebKit の電話を含む 18 件が通った記録。2 回目の M1 は 3 ブラウザ・両サイズで直った：切った一覧（状態 09・20 の Run 2）の見出しは「100 of 260 listed」（コントラスト約 5.2:1）、「select all listed」、foot は「The first 100 of 260 are listed.」と「Show all 260」。切らない一覧は「3 shown」「select all」のまま（07・19・21・24・25）。L1：電話の「Filter Results」は条件が無い間は閉じた「▸」の disclosure で、電話の結果のページは 221 px 短くなり、状態 18 の最初の画面にタブの行、Graphic Summary の見出しと色の凡例が入る。L2：09・21・25 の一覧は 3 ブラウザとも「#」の列から始まる（2 回目の Chromium の 09 の 29 px の横 scroll も 0）。I1：デスクトップの一覧は 22 行で scroll（09・20 で +336 px）、電話は 10 行。I2：`GraphicSummary.vue` の注を直した。I5：34・37 は「the record left out of it is not in it」。I3・I4 は直していない（下の「残件と注意」）。それ以外は 2 回目の記録と比べて、修正 4 回目が意図して変えたもの（と実行ごとに変わる値：時刻、エンジンのメモリ、build のハッシュ、保存の時刻、Queue の時間、保存の行）の他は変わらない。電話の 11〜13 が 16〜17 px 高いのは、11 で利用者が絞り込みのパネルを開き、12・13 でも開いたままのため。W4 の規則（電話の幅 390 px、デスクトップの表が 950 px の列に収まる、新しい toggle は約 30 px）と文言の規則（「Run LOSAT」、NCBI の branding 無し、LOSAT Web の形式の注）は守られている。WebKit の電話の 43 状態はすべてあり、結果の状態はどれも意図した Run を出す（ゲート 2 の誤った click は再び起きていない）。参考：(I1) 電話 18 の最初の画面から、グラフの最初の棒はまだ約 190 px 下（前は約 410 px）。「Open results」でタブの行まで scroll すれば入る。(I2) 状態 30 の「the records left out of the chosen subject are not in it」は 6 つのうち 1 つを除いた場合も複数形（`src/domain/reproduce.ts:231` の数の分からない分岐で、意図して複数形）。(I3) 閉じた電話の絞り込みのパネルの「(N set)」を写した画面の記録は無い（E2E が確かめる）。(I4) 2 回目の I3・I4 は記録のとおり変わらない。記録は タスクフォルダの `reviews/screen-review-3.md`、比較の script と切り出しは `tmp/review-screens-3/`

## 見つけて直したもの

- コードレビュー 1 回目の M1・M2・L1〜L5、2 回目の L1〜L3（上）。
- 画面レビュー 1 回目の M1・L1〜L3・I1〜I5、2 回目の M1・L1・I1・I5（上）。
- 修正 3 回目で見つかった：読み込んだ Run の Run details の保存ボタンが、原 FASTA をつなぐとすぐ追随しない（パネルが古い Run を読んでいた）。落としていた E2E を足した（`948a75e3`）。
- 修正 2 回目で見つかった：JSON の 80〜138 s、レポートの Firefox 25 s、セッションの WebKit の 517 ms のフレームの間隔（上の「実測」）。
- ゲート 1 で見つかった：1 ページの試験の 1 px の判定（`d3342a0e`）。
- ゲート 2 で見つかった：WebKit の電話の画面の記録の click が別の行に当たる（判断 17）。計測の spec が修正 4 回目の見出しの文言を読めず、query の場面が記録されない（上のゲート 2。`08e4557b` で直し、記録が空のまま試験が通らないようにした。ゲート記録の下書きを書いた agent が `measure/results-*.json` の `error` から見つけた）。
- WP-A：新しい書き出しはどれも `Downloader.save` を使わず、Writer の `open` を通る。

## 合流

`feature/losat-web-gui` への merge は、このセッションでは行わない。エンジン側が行う。merge のときにエンジン側が行うこと：

1. `origin/feature/losat-web-gui-app` を `feature/losat-web-gui` に merge する（アプリ側の起点は W5 のゲート記録の `306b458d`。エンジン側が `141615ec` から進んでいれば通常の merge）。エンジン側のパス（`LOSAT/`、`web/adapter/`）、計画、README の表、`package.json` は、このブランチで変えていない。
2. README の表の S15 の行を「完了」にして、この記録へのリンクを置く（計画 §7 の S15 の行と、セッションの README の S15 の行も、W5 のときと同じ扱いで）。
3. 計画に、保守者の判断 1・2・12・13 の DW を書く（アプリ側は計画を変えない）。1：セッションに候補とメモを含める（既定で含める）。2：reader の kind が変わると除外を消して言う。12：結果は classic の 1 ページ。13：NCBI のように切る。
4. `docs/web/abi_v2.md` の §4（`losat_web2_scan_begin` の行）と §9（*scan* の kind 0 の項）の「アプリが kind 1・2 に切り替えるまで残す」はまだ kind 0 を言っている。S14 の残件で、アプリは S14 で切り替えた（S15 でも変えていない）。エンジン側がその文を書き直す。
5. 判断 3〜11・14〜17 は推奨案で進めた（Owner-delegated）と記す（W4・W4b・W5 の判断と同じ扱い）。

## エンジン側への申し送り

- **reactor とネイティブの CLI**：`141615ec` の木から作ったもの（W5 と同じ）。アプリの試験（`extraction-engine.test.ts`、`hsp-correspondence.test.ts`、`engine-runtime.test.ts`、`session.spec.ts` のエンジン入り）は、エンジンの reader（SF）、整列文字列、HSP レコード（ABI の stream 1）に依存する。SFe など reader・scan・`register`・HSP レコードを変える変更の後は、アプリ側で reactor を作り直してこれらを実行する。HSP レコードの field を足すと、セッションの読み込みは追加の field を許す（修正 1 回目）が、manifest の field は許さない。
- **エンジンの build の名前**：読み込んだ Run の検証バッジは `record.engineBuild` をこのサイトの build 名（`composition.ts` の `engineBuildName`）と比べる。エンジンを替えた release で古いセッションを開くと level は `outside` になる。これは意図した動き。
- **E2j（S13+）の後のアプリ**：S15 は E2j が merge される前に終わったので、E2j の指示書の 6.（アプリへの申し送り）は行っていない。E2j が merge された後のアプリ側の最初の作業が行う（W4b・W5 の記録の「エンジン側への申し送り」の内容のとおり）。
- **BLASTX**：座標・単位・frame の規則と抽出は合成のレコードの単体試験だけ。設定ファイルは BLASTX を拒否する。SX が merge されたら、BLASTX の検索画面、実検索の抽出の試験、設定ファイルと再現コマンドの BLASTX を足す。
- **`approvedExceptions`**：再現のパネルが付ける承認済みの例外の文は、`LOSAT/` の PD-TLOSAN-LOCAL-GENCODE-32 と同じ内容（既定以外の subject の遺伝暗号）。エンジン側の承認の範囲が変わったら、`reproduce.ts` の表を合わせる。
- **エンジンのメモリ不足**（W4・W4b・W5 の続き）：4000 写しの自己検索が 3 ブラウザで `memory allocation of 49861 bytes failed`。W6 の計測でも再現した。

## S16 への申し送り

[S16 の指示書](../../losat_web_gui_sessions/session_s16_w7_delivery.md)の「S15（W6）から引き継ぐこと」に書いた（最後の木と件数、ゲートの script と段階と時間、書き出しとセッションの研究データの扱い、Service Worker が入れてはいけないもの、残す計測、画面の記録の場所、残件、E2j）。

## 残件と注意

- **WebKit の書き出しの開始の間隔**：WebKit は書き出しごとの最初に 128〜155 ms フレームを描かない（CSV・JSON・レポート）。WebKit 自身の仕事で、ページの script の外（判断 14）。S17 の V-MOB の実機で見る。
- **WebKit の Alignments の節の遅い読み**（判断 17）：WebKit は scroll anchoring を持たないので、Run 3 の Alignments の節が Queue の上で遅れて読まれると、下の内容が押し下げられる。画面の記録の spec は event で click を送って避けている。アプリは節の高さを読む間確保しない。直すなら placeholder の高さ（レイアウトの変更）。
- **Firefox の大きなセッションの保存のフレームの間隔**：ゲート 3 の Firefox の 3,000 写しのセッション（約 87 MB）の保存で、最長のフレームの間隔が 633 ms だった（ページの task の遅れは 54 ms でページの script の外。ゲート 2 の同じ場面は 50 ms）。1 回ずつの値で、追っていない。Chromium の JSON の 117 ms（garbage collection の停止）も 100 ms の目標を超える。
- **計測の spec の失敗の扱い**：`08e4557b` から query の場面に記録された失敗は試験を落とす。1 組の場面は見込みどおりの失敗（メモリ不足、保存の上限）があるので、今も記録だけで、結果の JSON の `error` を見る必要がある。
- **コードレビュー 2 回目の情報**：null の score が最良の E value として並ぶ、貼り付けた入力の `query.fa` の注、隠れた tab の遅さ（上）。
- **W5 コードレビュー L3**：Data worker の行の対応（`message-line.ts`）はアダプタの先読みの状態を持たない。混ざった行末で、除外の按が隣のレコードを名指し得る。直していない。
- **画面レビュー 2 回目の参考 I3・I4**：グラフの行が 8 px（タッチの標的が小さい）、電話のタブの行が 3 + 1 に折れる。直していない。
- **索引の時点の拒否**（W5 の判断 4）：NCBI の reader が拒む入力（gap の行、`Near line N`、Seq-id の最初の行）は索引の時点で失敗し、除外の按は出ない。設定ファイルやセッションを読み込むときも同じで、読み込みで索引が失敗することがある。
- **BLASTX の検索は SX まで動かせない**：抽出と出力は合成のレコードの単体試験だけ（実検索の試験は無い）。
- **画面レビュー 3 回目の参考**：電話で「Open results」の後の最初の画面に、グラフの最初の棒がまだ約 190 px 届かない（タブの行まで scroll すれば入る）。状態 30 の再現のコマンドの注は、除いた record の数が分からない分岐で複数形のまま。閉じた電話の絞り込みのパネル（「(N set)」）の画面の記録は無い（E2E が確かめる）。
- **E2j の申し送りは未実施**（上）。
- **W4b・W5 から残る残件**：エンジンのメモリ不足（4000 写し）、WebKit の保存の上限（OPFS の無い WebKit は結果をメモリに置き、512 MB を超える Run は「Not enough temporary storage」で失敗する。V-MOB の実機で確かめる）、`runInputs` の解放（Run の削除が入るまで解放しない）、iOS Safari での確認（V-MOB、S17）、TBLASTN の電話の「Axes not to scale」と小格子の密さ、`results_columns.md` の列の表（E2j の後に対応表と合わせる）、撮らなかった NCBI の画面。
- **画面の記録（PNG）はリポジトリに入れない**：ゲート 3 の画面の記録は `$BUILD_ROOT/s15-gate-screens/`（`/home/kawato/.cache/losat-work/s15-gate-screens/`。258 枚、SHA-256 は実行記録の `screens.sha256`）。途中の記録は `s15-wpg-screens`（`4495f852`、258 枚、`logs/wp-g/screens-2.sha256`）、`s15-wph-screens`（`4fbace1d`）、`s15-fix3-screens`、`s15-fix4-screens`、`s15-g2fail-screens`。状態の一覧：01〜28 は W5 と同じ並び、29〜40 と 30b が W6（Outputs の LOSAT Web formats、再現のコマンド、承認済みの例外（30b）、設定ファイル、Session のパネルと保存、読み込んだ Run、原 FASTA のつなぎ直しと拒否など）。ゲート 2 の画面の記録の落ちた 1 件は `s15-g2fail-screens` に残る。
- **ゲートの script の順序**：`run_gate.sh` の `all` は、結果画面の E2E の繰り返しで落ちると、計測と画面の記録に進まない（ゲート 1）。
- 背景の作業の記録（担当の報告、レビュー、判断の一覧）はタスクフォルダ `/home/kawato/losat-baselines/s15-w6-20261010/`（リポジトリには入れていない）。
