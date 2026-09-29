# LOSAT Web E2a-1（Session S06）ゲート記録

- 段階：E2a-1 BLASTN outfmt 0：権威と fixture（[総合計画書](../../losat_web_gui_plan.md) §7 の S06、[指示書](../../losat_web_gui_sessions/session_s06_e2a1_blastn_outfmt0_authority.md)）
- ブランチ：`feature/losat-web-gui`。変更前は S05 の記録の後の `ccdc471e0`。LOSAT のコードは変えていない（エンジンは S05 の `edf825d36` のまま）。この記録と fixture は同じコミットに入る
- 権威の記録：[`AUTHORITY.md`](AUTHORITY.md)（英語。NCBI の経路の表、LOSAT の部品の対応、見つかった不具合、S07 のゲート）
- 実行記録：[`run-20260929T000046Z/`](run-20260929T000046Z/)。ファイルのハッシュは [`evidence.sha256`](evidence.sha256)（リポジトリのルートで `sha256sum --check docs/evidence/losat_web_e2a/evidence.sha256` を実行する）
- 判定：**完了条件を満たした**。ただし、fixture を確かめる途中で、outfmt 0 とは別の既存の不具合を 5 つ見つけた（下の「見つかった既存の不具合」）。どれも推奨案で扱いを決め（S07 と、新しい S07+）、TLOSAN v0.2.0 の認証済みのコードに関わる 1 つを保守者に報告する

## 成果物

| ファイル | 内容 |
|---|---|
| [`AUTHORITY.md`](AUTHORITY.md) | NCBI の `blastn -outfmt 0` の呼出し経路（`blastn_app.cpp` → `CBlastFormat` → `showdefline.cpp` / `showalign.cpp`）をファイルと行の番号と逐語の断片で表にしたもの（§A）、BLASTP の報告との違い（§B）、`LOSAT/src/report/pairwise.rs` の部品の使い道（§C）、見つかった不具合と扱い（§D）、oracle の出力の例（§E）、S07 のゲートと拒否を続けるオプションと移植の順（§F） |
| [`LOSAT/tests/outfmt0_manifest.tsv`](../../../LOSAT/tests/outfmt0_manifest.tsv) | fixture の manifest。34 件（BLASTN 32、BLASTP 1、TBLASTN 1）。列は ID、program、task、query、subject、追加の引数、stdout の SHA-256 と大きさ、stderr の SHA-256、確かめる内容。すべて `LOSAT/` から実行する。`program` の列があるので、S08 の TBLASTX の fixture も同じ manifest に足せる |
| [`LOSAT/tests/fixtures/outfmt0/`](../../../LOSAT/tests/fixtures/outfmt0/) | 固定した NCBI の出力（`<ID>.out`、stderr があるものは `<ID>.err`）。合計 1.8 MB |
| [`LOSAT/tests/fasta/outfmt0/`](../../../LOSAT/tests/fasta/outfmt0/) | 作った入力 19 ファイル（`make_inputs.py` が固定の乱数の種からバイト単位で同じものを作る） |
| [`make_inputs.py`](make_inputs.py) | 入力を作る |
| [`run_oracle.py`](run_oracle.py) | manifest の全件を NCBI で実行する。既定は manifest のハッシュとの照合、`--freeze` は manifest のハッシュの列を書く（S06 だけ）。NCBI の出力を変える環境変数と `.ncbirc` があると実行しない。`search_argv` は S07 が LOSAT に同じ引数を渡すのに使う |
| [`precheck_hits.py`](precheck_hits.py) | 各 fixture を outfmt 6 で LOSAT と NCBI の両方で実行して比べる（outfmt 0 で一致するには、同じ HSP が見つかっている必要がある） |
| [`scoring_sweep.py`](scoring_sweep.py) | 既定以外の reward / penalty / gap の組合せを outfmt 6 で比べる（§D.5） |

## 完了条件と結果

| 完了条件（計画 §7 の S06） | 結果 | 証拠 |
|---|---|---|
| 経路の対応表 | 作成した。`AUTHORITY.md` §A の表の NCBI の断片は、引用した行の範囲にそのまま含まれることを機械的に確かめた（出力の文字列やオプションの名前を除く）。引用した NCBI のファイルは、固定 commit `598d8ae6` と改行以外に差が無い | `AUTHORITY.md` §A |
| 固定した fixture と SHA-256 | 34 件を NCBI BLAST+ 2.17.0 で固定した。2 回目の実行で全件が同じバイトになった。最長の実行は 0.43 秒 | `run-20260929T000046Z/oracle-freeze.log`、`oracle-check.log`、manifest |
| （指示書 5）S07 のゲートと、明示的に拒否するオプションの一覧 | `AUTHORITY.md` §F に書いた。一覧のオプションはすべて LOSAT の解析器がすでに拒否することを `edf825d36` で確かめた | `AUTHORITY.md` §F |

### fixture が確かめる範囲

指示書 3 の各項目と fixture の対応（詳細は manifest の `covers` の列）：

| 項目 | fixture |
|---|---|
| plus / minus 鎖 | `strand.*`、`multi.*`、`LC738874_LC738870.megablast`、`AP027152_AP027202.blastn` |
| 1 つの subject に複数の HSP | `strand.*`、`multi.*`（`msD` は plus と minus を 1 つずつ）、2 つのゲノムの組 |
| 複数の query と subject | `compact.*`（2 × 2）、`multi.*`（3 × 6。ヒットの無い query を挟む）、`many.*`（1 × 260） |
| ヒット無し | `compact.nohit.*`、`multi.*` の `mq2`、`edge.*` |
| `-lcase_masking`、`-dust` | `mask.*`（DUST、`-dust no`、query と subject の小文字、その組合せ）、`edge.alllower_lcase.blastn` |
| ギャップ | `strand.*`（行頭のギャップ）、`LC738874_LC738870.megablast` |
| 長い defline | `longdef.megablast`、`strand.*`、`LC738874_LC738870.megablast` |
| `-max_target_seqs` の境界 | `multi.mts1/mts3`（5 未満の警告）、`many.blastn`（説明 260・アラインメント 250）、`many.mts255`、`many.mts500`（省略と 500 の区別） |
| task megablast / blastn | ほぼすべての入力を両方の task で |
| そのほか | 座標の桁の境界（`width.*`）、無効な query（`edge.allN.blastn`）、HSP の選別のオプション（`multi.maxhsps1/perc97/besthit`）、既定以外の得点と 0 のギャップの末尾の式（`multi.r2p3g00.megablast`） |

入力は `LOSAT/tests/blastn_parity_manifest.tsv` の既存の入力（compact と 2 つのゲノムの組）をまず使い、既存の入力に無い場合（鎖、マスク、境界など）だけ小さな入力を作った。

除外した候補：`LC738873_LC738871.blastn`（7.2 MB、11,668 HSP）、`LC738874_LC738870.blastn`（1.3 MB、2,729 HSP）、`LC738874_LC738870.nodust.megablast`（1.1 MB、3,454 HSP）。どれも同じ書式の分岐を小さな fixture が確かめているので、リポジトリの大きさを優先した。outfmt 6 では後の 2 つも LOSAT と NCBI が一致した。

## 見つかった既存の不具合（outfmt 0 の writer の外。詳細と NCBI の根拠は `AUTHORITY.md` §D）

扱いは、セッションの規則（保守者に尋ねず推奨案で進め、記録する）に従って決めた。

| | 内容 | 扱い |
|---|---|---|
| D.1 | BLASTP と TBLASTN の outfmt 0 の座標の桁数が、NCBI（0 始まりの最大値）と違い 1 始まりの最大値から求められている。最大の座標がちょうど 10 の累乗のとき、空白が 1 つ多い（`width.blastp`、`width.tblastn`）。BLASTX は正しい | S07 で、NCBI の規則の 1 つの関数にまとめ、BLASTP・TBLASTN・BLASTN（S08 で TBLASTX）から使う。BLASTX は変えない（DW-10）。変わる凍結ハッシュはすべて列挙し、それぞれ NCBI の出力に一致することを示す |
| D.2 | BLASTN の `-task blastn -word_size 7`（または 8）で、ギャップ付き伸長の開始点を探す処理が subject の先頭を越えて panic する（`gapped.rs:559`）。NCBI は subject の両端の番兵のバイトで止まる | S07 の最初の作業（エンジンの変更、回帰試験は `multi.ws7_e1000.blastn`）。同じ修正を試しに入れたビルドは、この入力の word size 4〜8 と 34 件の fixture すべてで outfmt 6 が NCBI とバイト一致した（`probe-gapped-start-bound.diff`、`precheck-probe.tsv`。worktree のファイルは元に戻した） |
| D.3 | BLASTN は NCBI の警告（`-max_target_seqs` が 5 未満、無効な query）を stderr に出さない（どの形式でも） | S07 で足し、stderr を記録した 3 件で stderr もゲートにする |
| D.4 | BLASTN の `-max_target_seqs` は既定値 500 を持つので、省略と 500 を区別できない（outfmt 0 の表示数に必要） | S07 |
| D.5 | 認証済みの既定の得点（blastn 2/−3・5/2、megablast 1/−2・0/0）以外で：NCBI が拒否する組合せ 34 を LOSAT は受け付けて実行し、NCBI が受け付ける組合せのうち 17 は結果が違う（違うものはすべて reward 1 で gap が 0 でない） | 新しいセッション **S07+**（段階 E2c、S07 の直後）で、NCBI の検査（同じ拒否と文言）を移植し、違いの原因を調べる。終わっても一致しない組合せは明示的に拒否し、アプリが認証されていない BLASTN の得点を出さないようにする。S07 の fixture は認証済みの範囲に留める |

**保守者への報告（判断は不要）**：D.1 の TBLASTN の桁数は、TLOSAN v0.2.0 で認証した TBLASTN の outfmt 0 と同じコードにある。認証の fixture はこの境界に達しないので、認証の主張の範囲は変わらない。このブランチで S07 に直す。

## 実行の方法（再現）

```bash
# 入力を作る（既存のファイルと同じバイトになる）
python3 docs/evidence/losat_web_e2a/make_inputs.py
# NCBI の出力を manifest と照合する（固定するときだけ --freeze と --out LOSAT/tests/fixtures/outfmt0）
python3 docs/evidence/losat_web_e2a/run_oracle.py --bin-dir /home/kawato/micromamba/bin --out <dir>
# LOSAT と NCBI の outfmt 6 の比較、得点の組合せの比較
python3 docs/evidence/losat_web_e2a/precheck_hits.py --bin-dir /home/kawato/micromamba/bin --losat <LOSAT>
python3 docs/evidence/losat_web_e2a/scoring_sweep.py --bin-dir /home/kawato/micromamba/bin --losat <LOSAT>
```

`<LOSAT>` は `edf825d36` から作ったネイティブの release の実行ファイル（`cargo +1.92.0 build --release --locked`）。

## 引き継ぎ（S07 へ）

- S07 の指示書に、`AUTHORITY.md` §F の移植の順（D.2 の修正と D.1 の関数から始める）とゲートを書き足した。
- 新しいセッション S07+ の指示書を作り、計画 §7 と README の表に行を足した。
- S08 は、TBLASTX の fixture をこの manifest（`program` の列）に足し、`run_oracle.py` で固定する。
