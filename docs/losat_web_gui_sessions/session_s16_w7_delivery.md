# Session S16 — W7：配信の仕上げ

## INSTRUCTION PROMPT

LOSAT Web の段階 W7 を実行する。先に [セッション README](README.md) の共通規則を読み、それに従う。完了条件の正本は、総合計画書 §7 の S16 の行である。設計は計画 §5.9〜§5.10、設計書 §13〜§14 にある。Cloudflare へのデプロイと GA4 の設定は保守者が行う（README の規則 7）。

1. release ごとの資産：資産を `/r/<release-id>/…` に置き、小さな起動用の `index.html` が、`localStorage` に固定した release ID の資産を読み込む（読み書きは try/catch で囲む）。資産の一覧とハッシュを ReleaseManifest として生成する。
2. Service Worker：固定した release の資産（2 つの wasm、Worker、例題を含む）を取得してハッシュを確かめ、それから「Offline ready」を出す。研究データは Service Worker のキャッシュに入れない。
3. 明示更新：新しい release を知らせる → 利用者が選ぶ → 資産を取得して確かめる → 実行中・未保存の作業への影響を案内する → 固定 ID を切り替えて再読み込みする。取得が不完全なら今の release のままにする。`skipWaiting` を無条件には使わない。旧版の資産は、それを使っているタブがある間は消さない。
4. GA4：計画 §5.10 の (a)（別 origin の計測文書を iframe で埋め込む）と (b)（自分の Cloudflare Worker から Measurement Protocol へ転送する）を、COEP の下で対象の全ブラウザで試し、結果とともに保守者に推奨を示す。決まった方式を実装する。初期のイベントはアプリの訪問と release ID だけで、同意の前は何も読み込まず、何も送らない。計測が塞がれていても、未同意でも、オフラインでも、検索・結果・抽出・保存が動くこと。
5. V-PRIV（計画 §5.10）と V-OFF（計画 §6.2）を Playwright で作り、CI に入れる。V-PRIV には、配信する Wasm と JavaScript にビルドした機械のパス（checkout のパス、`CARGO_HOME`、ユーザー名）が含まれないことの検査を含める（計画 TD-11）。
6. 導入支援：「例を試す」、5 program の小さな例題、最短のチュートリアル、任意のブラウザ自己試験（期待値の SHA-256 と比べる、V-BR の小さな版）。
7. Cloudflare のプレビューの設定（`web/app/public/_headers` をそのまま使う）を用意し、保守者がデプロイした後に、プレビューで `crossOriginIsolated === true` と CSP を確かめる。304（条件付きの要求への応答）にも COOP / COEP / CORP / CSP が付くことと、WebKit / Safari で取消の後の検索が動くことも確かめる（WebKit は、304 に COEP が無いと、作り直す Engine worker の読み込みを拒否する。[W3 のゲート記録](../evidence/losat_web_w3/README.md)の判断 12。開発とプレビューのサーバーは `vite.config.ts` の plugin で付けている）。ハッシュ付きの資産（`/assets/*`、`/r/<release-id>/…`）に immutable のキャッシュを付けることも検討する。
   旧版のタブのデータを新しい版が消さないこと（`tmp/` と Web Lock の名前 `losat-web:tmp:` を変えないこと。設計書 §7.2、S10 の判断 8）を E2E で確かめる。

完了条件は計画 §7 の S16 の行による。画面の記録を画面レビューに見せ、`docs/evidence/losat_web_w7/README.md` に記録する。

## S15（W6）から引き継ぐこと

[W6 のゲート記録](../evidence/losat_web_w6/README.md)の要点。上の 1.〜7. の条件はこの節で緩めない。

- **最後の木と件数**：最後のアプリの木は `08e4557b`（ゲート 3 の木。アプリの本体は `a00e4a7e` から変わっていない）。`npm run check` 656 件（61 件 skip）、reactor 付きの単体 716 件（1 件 skip。HSP の対応の試験は W4 から変えていない）、FakeEngine のビルドの E2E 146 件（Chromium 53・Firefox 47・WebKit 46）、エンジン入りのビルドの E2E 170 件（61・55・54。V-BR を含む）、結果画面の E2E の繰り返し 90 件、画面の記録 258 枚（3 ブラウザ × デスクトップ・電話、状態 01〜40 と 30b）、NCBI BLAST+ 2.17.0 との比較コマンド 42 件（33 一致、承認済みの例外 6、NCBI が拒む 3）。ゲート 3（`08e4557b`、[`run-20261010T170931Z/`](../evidence/losat_web_w6/run-20261010T170931Z/)）はすべての段階が通った（約 55 分）。画面レビュー 3 回目はその記録で合格（High・Medium・Low 無し、参考 4 件）。
- **ゲートの script**：[`run_gate.sh`](../evidence/losat_web_w6/run_gate.sh)（W5 のコピーに NCBI の比較の段階を足したもの）と [`check_commands.py`](../evidence/losat_web_w6/check_commands.py)。段階：V-ABI quick、`npm ci`、`npm run check`、reactor 付きの単体試験、`check_commands.py`（oracle の lock の下）、FakeEngine のビルドの E2E、エンジン入りのビルドの E2E、結果画面の E2E の繰り返し 2 回、計測、画面の記録。`LOSAT_WEB_GATE_STEPS` は `all`・`after-review`・`measure`。`all` には `LOSAT_WEB_NCBI_BIN`（NCBI BLAST+ 2.17.0 の bin。`~/micromamba/bin`）も要る。全段階で約 55 分（ゲート 3：02:09〜03:04 JST。計測の Firefox が 6.7 分）。繰り返しで落ちると計測と画面の記録に進まない。S16 のゲートは記録の場所を `losat_web_w7` に替えたコピーを作る。reactor と native は `141615ec` の木のもの（エンジンが進んでいれば作り直す）。
- **研究データの扱い（書き出しとセッション）**：書き出し（互換出力、CSV、JSON、レポート、SVG、設定ファイル、セッションファイル、抽出した FASTA、Run の入力 FASTA）は、すべてページの中で作り、利用者が download として保存する。どこにも送らず、URL にも載せない。CSV・JSON・レポート・セッションの整列行は研究データ（配列を含む）。セッションファイルは OPFS のパス、作業のトークン、Run・group・revision・source の ID を書かず、入力の FASTA も入れない（入力の SHA-256・レコード表・除外だけ）。候補とメモは既定で入る（保守者の判断。チェックを外せる）。file 名は入力の名前を含まない（`losat-session-YYYYMMDD-HHMMSS.losat-session.gz`、`losat-run{N}-{program}-hsps.csv` など）。GA4 などの計測のイベントに、これらの file の名前・内容・大きさ・件数を載せない。V-PRIV はセッションファイルの内容にパスやトークンが無いことを検査に含められる（`session-save-load.test.ts`「writes no storage path, token, run, revision or source ID of the session that saved it」が単体の検査）。
- **Service Worker が入れてはいけないもの**：セッションファイルは利用者が file として選び（`<input type="file">`）、書き出しは Blob の download で、どちらもネットワークの応答ではない。Service Worker のキャッシュには release の資産（2 つの wasm、Worker、JavaScript、CSS、例題）だけを入れ、実行時の取得をキャッシュしない。利用者の file、書き出した file、Blob URL、セッションの内容を `fetch` の応答として扱わない。ドットプロットの SVG とレポートの HTML は download の file で、ページの URL ではない（レポートは自分の CSP の `<meta>` を持つ。サイトの CSP と `_headers` は別）。例題（S16 の「例を試す」）は release の資産であって利用者のデータではない。
- **セッションの版と release**：container 1・schema 1（形式は [`docs/web/session_file.md`](../web/session_file.md)）。manifest の field を足す・変えるときは新しい schema にする（未知の field は拒否する）。読み込んだ Run のバッジは、記録された `engineBuild` がこのサイトのエンジンの build 名（`composition.ts` の `engineBuildName`）と同じときだけ通常のバッジで、別の build や不明は `outside`（「Written by another engine build」）にする。release ID と build 名の付け方を変えるときは、この比較が壊れないようにする。旧版のタブのデータ（`tmp/` と Web Lock の名前 `losat-web:tmp:`）は W6 で変えていない。
- **残す計測**：10 万 query（149,828 HSP）の書き出しとセッション（Chromium・Firefox・WebKit、`157488f3`、`logs/fix2/report.md`）：CSV 0.5〜0.7 s、JSON 1.8〜2.9 s、レポート 1.9〜2.5 s、セッションの保存 2.8〜3.5 s・開く 1.8〜3.1 s、最長のフレームの間隔は 100 ms 以下（WebKit の書き出しの開始の 128〜155 ms を除く）。1 ページの結果で 10 万 query の Run を開くのは 276〜306 ms（W5 は 281 ms）。3,000 写し（5,993 HSP）のセッションは約 87 MB で、ゲート 3 では Chromium の保存 8.2 s・開く 7.7 s、Firefox の保存 14.9 s・開く 16.3 s（ゲート 2 の Firefox の開くは 25.5 s。1 回ずつの値で揺れが大きい）。ゲート 2 の query の場面は記録が空だった（計測の spec が修正 4 回目の「100 of 260 listed」を読めなかった）。`08e4557b` から、query の場面に記録された失敗は試験を落とす。計測の JSON の `error` も見ること。ゲート 3（`08e4557b`）で 3 ブラウザとも測り直し、書き出しとセッションは修正 2 回目と同じか速かった（10 万 query で JSON 1.75〜1.79 s、セッションの保存 1.9〜2.3 s・開く 1.4〜1.9 s。Run を開くのは Chromium 258・Firefox 554・WebKit 410 ms）。Firefox の 3,000 写しのセッションの保存に 1 回 633 ms のフレームの間隔があった（ページの script の外）。S16 の release の計測は、今の木から取り直す。
- **残件（S16 に関わるもの）**：
  - WebKit は書き出しごとの最初に 128〜155 ms フレームを描かない（WebKit 自身の仕事）。
  - WebKit は scroll anchoring が無く、Alignments の節が遅れて読まれると下の内容が動く（判断 17。画面の記録の spec は event で click を送って避けている）。節の高さを読む間確保するのは未実施。
  - 隠れた tab では timer が約 1 s に絞られ、書き出しとセッションの保存が遅くなる。
  - 貼り付けた入力は、セッションに `query.fa` / `subject.fa` として記録される。
  - W5 コードレビュー L3（行の対応が先読みの状態を持たない）、画面レビュー 2 回目の I3・I4（グラフの行 8 px、電話のタブの行）は直していない。
  - BLASTX の検索は SX まで動かせない。抽出と出力は合成のレコードの単体試験だけで、設定ファイルは BLASTX を拒否する。
  - `abi_v2.md` の kind 0 の文はエンジン側が書き直す（S14 の残件。合流の項）。
  - エンジンのメモリ不足（4000 写し）、WebKit の保存の上限（512 MB）、`runInputs` の解放、iOS Safari（V-MOB、S17）。
- **画面の記録と画面レビュー**：S16 の画面レビューは W6 の最後の記録と比べる。`$BUILD_ROOT/s15-gate-screens/`（`/home/kawato/.cache/losat-work/s15-gate-screens/`、ゲート 3 の木 `08e4557b`、258 枚）。SHA-256 は [`run-20261010T170931Z/screens.sha256`](../evidence/losat_web_w6/run-20261010T170931Z/screens.sha256)。状態：01〜28 は W5 と同じ並び（結果画面は 1 ページの Descriptions）、29〜40 と 30b が W6（Outputs、再現のコマンド、設定ファイル、Session のパネル、読み込んだ Run、原 FASTA のつなぎ直し）。
- **E2j（S13+）**：E2j は S15 の時点で merge されておらず、E2j の指示書の 6.（アプリへの申し送り）は S15 でも行っていない。E2j が merge された後のアプリ側の最初の作業が行う（W4b・W5 の記録の「エンジン側への申し送り」のとおり。Filter Results の Query Coverage、BLASTN の 1 文字の HSP の鎖、結果の列の置き場所）。鎖が分かれば、トレイの「Strand not decided」と抽出の「unknown strand」の扱いが変わる。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S17 — 公開判定](session_s17_g_release_decision.md)。
