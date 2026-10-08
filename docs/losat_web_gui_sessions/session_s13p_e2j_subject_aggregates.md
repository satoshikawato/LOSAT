# Session S13+ — E2j：subject ごとの集約の値と HSP の鎖（エンジン側）

## INSTRUCTION PROMPT

LOSAT Web の段階 E2j を実行する（エンジン側。Linux の clone の worktree `/home/kawato/losat-work/.worktrees/web-gui`、ブランチ `feature/losat-web-gui`。場所は clone の `CLAUDE.local.md`、手順は skill `losat-campaign`・`losat-worktree`・`losat-gates`・`losat-oracle-runs`・`losat-ship-pr` に従う）。先に [セッション README](README.md) の共通規則を読み、それに従う（規則 1・4 の `/mnt/c` のパスは `CLAUDE.local.md` の表で読み替える）。完了条件の正本は総合計画書 §7 の S13+ の行（S13 の合流で足す）。この段階は S13（W4）が結果画面の列定義表（[`docs/web/results_columns.md`](../web/results_columns.md)）で採用を決め、TS では計算しない値をエンジンから出すためのものである（計画 §5.7、[W4 のゲート記録](../evidence/losat_web_w4/README.md)の判断 1・2）。順番はエンジン側の SF（E2h）の後、S09+（R2）の前を推奨する（保守者の確認待ち）。

1. **棚卸し（DW-12）**：アプリが出すオプションの範囲で、NCBI が subject ごとの集約の値を作って書く経路の関数を、BLASTN・BLASTP・TBLASTN・TBLASTX ごとに棚卸しし（忠実な移植・差のある移植・未移植・明示的な拒否）、未移植を簡略化せずに transpile する。少なくとも次を含む（固定 commit `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`）：
   - `CBlastFormat::x_ConfigCShowBlastDefline`（`c++/src/algo/blast/format/blast_format.cpp:497-527`）：`-sorthits` を与えると `eShowPercentIdent`・`eShowTotalScore`・`eShowQueryCoverage` を立てる。
   - `CShowBlastDefline` の説明表（`c++/src/objtools/align_format/showdefline.cpp`。テキストの経路 `:660-980` 付近、`x_GetScoreInfo` `:1483-1515`、`x_GetScoreInfoForTable` `:1518-1560`）。
   - `CAlignFormatUtil::GetSeqAlignSetCalcParams`（`align_format_util.cpp:4247-4316`。`highest_bits`、`total_bits`、bit score が最も高い HSP の E 値と identity（WB-1175）、`percent_coverage`（四捨五入、100 を超えれば 99））、`GetSeqAlignCoverageParams`、`GetAlignmentLength`、`GetPercentIdentity`、`GetScoreString`（`total_bit_score_buf`）。
   - `-sorthits` / `-sorthsps` の引数と並べ替え（`blast_args.cpp`、`CAlignFormatUtil` の `SortHitBy*` / `SortHspBy*`）。outfmt 6 の `qcovs`・`qcovhsp`（`tabular.cpp:707-735`）は、LOSAT が custom の列を出す program で比較の補助にできる。
2. **CLI と比較の oracle**：BLAST+ の CLI は `-sorthits` のとき outfmt 0 の説明表に Total Score・Query Coverage・Percent Ident を出すので、`-sorthits`（必要なら `-sorthsps`）を移植し、NCBI BLAST+ 2.17.0 の outfmt 0 とバイト一致で確かめる（推奨。CLI に無い値を ABI だけで出すと、バイトで確かめる手段が無い）。移植しない program・値は明示的な拒否にする。
3. **ABI v2 の追加**（`docs/web/abi_v2.md` と `web/app/src/ports/engine.ts` を同じコミットで変える。互換を壊さない追加にし、壊すなら版を上げる）：
   - query と subject の組ごとに、NCBI が書く文字列のままの Max score、Total score、Query cover、E value、Per. ident（と HSP の数）。stream 1 のレコードの種類を足すか、新しい stream にするかは、アダプタの都合で決めて記録する。
   - HSP レコードに鎖の欄（BLASTN の HSP の query の frame、+1 / −1）。BLASTN の 1 文字の HSP は start = end で、向きが座標から分からない（`docs/evidence/losat_web_e2c/AUTHORITY.md` §P、`abi_v2.md` §8）。既存の `query_frame` の意味（BLASTN は null）を変えないため、別の欄を推奨する。
4. **fixture**：
   - 集約の値：4 program × `-sorthits` 0〜4 の outfmt 0（1 subject に複数の HSP、minus 鎖、翻訳の frame、被覆が 100% を超える組（99 に切る）、bit score の同点、`-max_hsps`、範囲の指定）。
   - 鎖：AUTHORITY §P の例（query `TAGGACGG`、subject `YCAYAANTNCRGYACT`、`-task blastn -word_size 4`）と、plus 鎖の 1 文字の HSP。
   - **BLASTN の既定の outfmt 7**：`LOSAT/tests/outfmt0_manifest.tsv` に、既定（megablast）と `-task blastn` の outfmt 7 の fixture を足す。検証バッジの表（`web/app/build/verification.ts`）は今、既定の BLASTN の option の組を outfmt 0 と 6 でしか NCBI と比べた記録から読めず、既定の BLASTN の Run を「outside certified profile」と示す（W4 のゲート記録の残件）。fixture を足せば、生成する表が自動で直る。
5. **ゲート**：固定した fixture で NCBI とバイト一致。各 program の既存のゲート、Gate A、TLOSAN の Stage G、S07〜S11・SF の fixture と sweep、v1 の検査に退行なし。変えた program の V-ABI（新しい欄を含む）。V-PERF の非退行（集約は出力の時だけ）。棚卸しの表を基準にした独立監査。
6. **アプリへの申し送り**：結果画面は、エンジンが値を出したら列定義表の「採用・エンジン待ち」の列（Max score、Total score、Query cover、E value（表形式）、Per. ident）を NCBI の表形式の並びで出し、鎖の欄で BLASTN の 1 文字の HSP の向きを示す（`src/domain/result-index.ts` の `orientation`）。その作業は、この段階が merge された後の最初のアプリ側のセッションの最初に行う。

## 終了・引き継ぎ

README の規則 8 に従う。次はエンジン側の S09+（R2）。
