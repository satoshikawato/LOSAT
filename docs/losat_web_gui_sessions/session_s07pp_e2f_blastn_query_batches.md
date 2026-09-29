# Session S07++ — E2f：BLASTN の query の batch

## INSTRUCTION PROMPT

LOSAT の段階 E2f を実行する。BLASTN（LOSATN）が複数の query を、NCBI と同じ query の batch に分けて検索するようにするセッションである（TD-14）。先に [セッション README](README.md) の共通規則を読み、特に規則 4（`AGENTS.md`、`verify-ncbi-parity-and-speed`）に従う。完了条件の正本は、総合計画書 §7 の S07++ の行である。

現状（S07+ の第 7 回の独立監査、`docs/evidence/losat_web_e2c/AUTHORITY.md` §K）：NCBI の `blastn` は query を batch に分け、batch ごとに検索する（`c++/src/app/blast/blastn_app.cpp:261-300`）。最初の batch は約 5000 残基（`CBatchSizeMixer`、`blast_app_util.hpp:57-95`、`blast_app_util.cpp:65-81`）で、その後の batch の大きさは、前の batch の成功した ungapped 伸長の数（`CLocalBlast::GetNumExtensions`、`local_blast.cpp:301-307` の `good_init_extends`）から決まる。LOSAT は全 query を 1 つの塊として検索する。そのため、batch に依存する次の振る舞いが NCBI と違う：

- 近似の ungapped 伸長は、context ではなく query の塊の端で止まる（`na_ungapped.c:164,286,317`）。batch の境の query の端の seed で、LOSAT は次の query まで伸長し、既定のオプションでも HSP が増える（例：MelaMJNV[186013:191352] と `ACGTACGTAC` の 2 つの query、LC738874[178822:189290]、`-task blastn -outfmt 6` で NCBI の 62 行に 63 行。NCBI に `BATCH_SIZE=100000000` を与えると LOSAT と同じ）。
- lookup table の種類（`BlastChooseNaLookupTable` の query の長さ）、対角線の表か hash か（8000 残基）、表の大きさは、batch の query の塊から決まる。
- S07+ で batch に依存するために明示的に拒否しているもの（`scoring.rs` の `gap_x_dropoffs` の最初の batch を超える場合、`run.rs` の `unsearched_queries` が決まらない場合）。

S07+ の終わりに読んだ NCBI の経路（下の 1 で確かめ、記録する）：

- `CBatchSizeMixer`（`blast_app_util.hpp:57-95`、`blast_app_util.cpp:65-81`）：目標の hit 数は `max(1000000, subject の総文字数 / 3000)`、最初の batch は目標の 1/200（普通は 5000）。次の batch は `ratio = (hits + 1) / 前の batch の大きさ`（前の `ratio` があれば 0.3 と 0.7 で混ぜる）、`(Int4)(目標 / ratio)` を [100, chunk − 1000] に丸め、丸めたときは `ratio` を −1 に戻す。chunk は blastn 1000000、megablast 5000000（`local_blast.cpp:54-73`）。hits が 0 か 1 のとき `(Int4)` の変換は `int` を超える（x86 では `INT_MIN` になり、次の batch は 100 になるはず）。oracle の出力で確かめる。
- `good_init_extends`：`BlastNaWordFinder` の最後の `Blast_UngappedStatsUpdate(…, init_hitlist->total)`（`na_ungapped.c:1688`、`blast_diagnostics.c:102-115`。lookup の hit が 0 の subject の chunk は数えない）を、subject の chunk ごと・スレッドごとに足す（`blast_diagnostics.c:118-135`）。索引の megablast の経路（`na_ungapped.c:2144-2149`）は LOSAT に無い。
- query の分割（`CQuerySplitter`、`split_query_cxx.cpp:50-62`、`split_query_aux_priv.cpp:51-145`、`prelim_stage.cpp:230-260`）：batch の query の総文字数が 2 × (chunk − 100) 以上なら、NCBI は batch を重なり 100 の塊に分けて予備の段階を行い、HSP を合わせる。LOSAT に分割は無い。S07+ の終わりの確認（EDL933 の 2.1 Mb を query に、分割の境をまたぐ subject、`-task blastn`・megablast）では outfmt 6 が NCBI とバイト一致した。分割の境の近くの HSP で差が出るかを sweep で調べ、差があれば移植するか明示的に拒否する。

1. NCBI の batch の経路をソースで追い、`docs/evidence/losat_web_e2f/AUTHORITY.md` に記録する：`CBlastInput::GetNextSeqBatch` の詰め方、`CBatchSizeMixer::GetBatchSize` の式（`Int4` の変換を含む）、`good_init_extends` を数える箇所（`na_ungapped.c`）、`SplitQuery_GetChunkSize` の上限と query の分割、環境変数 `BATCH_SIZE`。
2. **batch ごとの検索を移植する。** `search` が NCBI と同じ batch に query を分け、batch ごとに query の塊・lookup table・対角線の表を作り、同じ subject を検索する。各 batch の `good_init_extends` を NCBI と同じ箇所で数え、次の batch の大きさを NCBI の式で求める。直す箇所の直上に NCBI のファイル・行と断片を書く。
3. S07+ の batch に依存する明示的な拒否を、移植した batch で置き換える（`gap_x_dropoffs`、`unsearched_queries`、outfmt 7 の件数、Karlin-Altschul の表の失敗の繰り返しの回数）。
4. 数の一致を確かめる。NCBI の batch の境は出力に直接は出ないので、境で結果が変わる入力（上の例と、境の前後に query の端がある入力の sweep）を作り、`BATCH_SIZE` を与えた NCBI（境を固定できる）と比べる。
5. 試験：複数の query の sweep（`docs/evidence/losat_web_e2c/slice_sweep.py` を複数の query に広げる）が NCBI とバイト一致。S07+ のすべての検査（`check_inputs.py`、`scoring_sweep.py`、`word_size_sweep.py`、`slice_sweep.py`）、既存の BLASTN のゲート、S07 の fixture に退行なし。V-PERF（NCBI は batch ごとに subject を検索し直す。`blastn-many` の変化を記録し、+5% を超えるなら原因と判断を記録する）。BLASTN の全升目の V-ABI。独立監査。

記録は `docs/evidence/losat_web_e2f/`（`README.md`、`evidence.sha256`）。変更前の成果物は S07+ の成果物である。

## 終了・引き継ぎ

README の規則 8 に従う。次は [S08 — TBLASTX outfmt 0/7](session_s08_e2b_tblastx_outfmt0_7.md)。
