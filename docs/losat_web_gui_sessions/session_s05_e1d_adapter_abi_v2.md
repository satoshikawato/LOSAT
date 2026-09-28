# Session S05 — E1d：アダプタと ABI v2

## INSTRUCTION PROMPT

LOSAT Web の段階 E1d を実行する。先に [セッション README](README.md) の共通規則を読む。新しい crate `web/adapter/` は `web/AGENTS.md` に従い、そこから使うエンジン側の変更が要るときはルートの `AGENTS.md` に従う。完了条件の正本は、総合計画書 §7 の S05 の行である。設計は計画 §4.6〜§4.8 と TD-1・TD-6、ABI の下書きは `docs/web/abi_v2.md` にある。

1. `web/adapter/` に Rust の cdylib crate を作る。`LOSAT = { path = "../../LOSAT", default-features = false }` とし、threaded 版では `parallel` と `wasm-threads` を明示して有効にする。reactor の CRT の設定は `LOSAT/build.rs` をパスで共有する（`build = "../../LOSAT/build.rs"` が使えるか確かめ、使えなければ最小の写しにする）。
2. ビルドの同一性（TD-6）：`[profile.release]`（`LOSAT/Cargo.toml`）と Wasm の rustflags（`LOSAT/.cargo/config.toml` の `+simd128` など）を adapter 側に写し、共通の依存の版・profile・rustflags・link の引数が LOSAT と一致することを検査するスクリプトを作る。結果を成果物の JSON に記録する（`LOSAT/tests/build_wasi_artifacts.py` が記録している項目に合わせる）。
3. `docs/web/abi_v2.md` の export を実装し、文書を実装に合わせて確定する（状態を draft から外す）。v1 の `losat_web_*` は同じモジュールに残し、状態は共有しない（TD-1）。`losat_host.emit` への 1 MiB ごとの書き出し、観測者の出来事からの `out6` / `out0` の範囲の記録、`describe`（clap の定義から生成し、program が対応する出力形式を含める）、`validate`、`register`（その program の解析器で解析し、レコードの ID と長さを返す）、`scan_*`（解析器の種類を受け取ってレコード表を作る。TD-8）を作る。`run` は program が対応する形式だけを出す。HSP レコード（`hits`）を求めると、エンジンはアラインメントを描画する（S02 の BLASTP の実装）。描画できない HSP があると、outfmt 0 を求めた CLI と同じく失敗する。Web は常に outfmt 0 を求めるので、この振る舞いは CLI の outfmt 0 と一致する。表形式だけを求める使い方を ABI で許すなら、そのときは `hits` を求めないこと（S02 の独立監査の指摘）。観測者が報告する outfmt 0 の節は、スコアの行とアラインメントだけで、各 subject の最初の HSP の前に NCBI が書く subject の見出し（`x_DisplayAlnvecInfo` から呼ばれる `x_ShowAlnvecInfo` の中、`showalign.cpp:3613-3632`）を含まない。詳細画面が見出しを原文のまま表示できるように、見出しの範囲を得る方法（例：観測者に subject の見出しの開始と終了を足す）をここで決めて、`abi_v2.md` に書く。スレッド版の共有メモリの最大値は、認証済みの threaded ビルドと同じにする（TD-7）。
4. `scan_*` の性質試験：LF / CRLF、行長の揃わないレコード、空行、非 ASCII のヘッダー、巨大なレコードの境界で、`bio::io::fasta` 型の `scan` の表が `bio::io::fasta` の結果と一致すること。NCBI 型（BLASTX の `LOSAT/src/algorithm/blastx/input.rs`）の性質試験は、BLASTX を扱う SX で足す。`register` での照合（食い違えば止める）の試験も置く。
5. 2 つの reactor（`losat-web-serial.wasm`、`losat-web-threads.wasm`）を、出力先を worktree の外にしてビルドする。
6. V-ABI：Node で 2 つの reactor を動かし（threaded は `LOSAT/tests/wasi_thread_host.js` と同じ方式）、`docs/web/verification_cells.tsv` の BLASTP・TBLASTN・BLASTN・TBLASTX の、その時点で対応する全形式の升目を、スレッド 1/2/4、Subject を保持した連続実行（query や条件を変える）で実行して、期待値と SHA-256 で比べる。TBLASTX の `wasm-threads` 専用の出力箇所を通る場合を含める。スクリプトは `web/adapter/tests/` に置き、CI（`.github/workflows/web.yml`。対象のパスには `LOSAT/` の変更がすでに含まれている）に足す。
7. `web/app/src/ports/engine.ts` の型が ABI と一致していることを確かめ、違えば同じコミットで直す。

完了条件は計画 §7 の S05 の行による。記録は `docs/evidence/losat_web_e1d/README.md`。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S06 — BLASTN outfmt 0：権威と fixture](session_s06_e2a1_blastn_outfmt0_authority.md)。確定した ABI と、V-ABI の実行方法を、S09（ブラウザでの実行基盤）の指示書に書き足す。
