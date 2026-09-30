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

## 終了・引き継ぎ

README の規則 8 に従う。集約列の移植を挿入した場合は、その指示書を書く。挿入しない場合、次は [S14 — 抽出と候補](session_s14_w5_extraction_candidates.md)。
