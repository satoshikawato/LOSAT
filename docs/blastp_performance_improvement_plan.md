# BLASTP Performance Improvement Plan

## 目的

LOSATP の BLASTP 経路を NCBI BLASTP に近い速度へ改善する。ただし、LOSAT の最優先要件は
NCBI BLAST+ との bit-perfect output parity であり、速度改善は常に NCBI C/C++ 実装の確認後に
行う。

この計画では、まず遅さの原因を測定で切り分け、その後に parity リスクの低い順で改善する。

## 絶対制約

- NCBI BLAST C/C++ 実装を唯一の挙動リファレンスにする。
- LOSAT runtime から NCBI BLAST+ バイナリ、ライブラリ、FFI、subprocess へ委譲しない。
- 実装変更時は、該当 Rust コードの直上に NCBI ファイルパスと行番号つきの参照コメントを置く。
- 出力順、score、bit score、E-value、座標、HSP pruning の結果を変えない。
- 並列化しても merge 順序を scheduler/channel 受信順に依存させない。
- 既知の parity divergence がある状態で、速度だけを目的にアルゴリズムを簡略化しない。

## 現時点の仮説

LOSATP が NCBI BLASTP より遅い主因は、Rust だからではなく、BLASTP 経路がまだ
parity-first の移植段階で、NCBI BLASTP の成熟した実行モデルと低レベル最適化に追いついて
いないためと考える。

優先度順の仮説:

1. BLASTP subject-level preliminary search がまだ十分に並列化されていない。
2. FASTA subject を毎回読み込み、encode し、preformatted BLAST database 相当の恩恵を受けていない。
3. amino-acid lookup, ungapped extension, diagonal table, gapped preliminary extension の hot loop が
   NCBI の cache-aware 実装ほど最適化されていない。
4. composition-based statistics / Kappa redo alignment が single-thread かつ allocation-heavy になっている。
5. pairwise output や alignment view 構築が検索後の時間を押し上げている可能性がある。

## 参照すべき NCBI ソース

実装前に、必ずローカル NCBI ソースで該当箇所の行番号を再確認する。

- `ncbi-blast/c++/src/algo/blast/core/blast_engine.c`
  - subject loop and HSP stream write: around `1409-1498`
  - per-thread preliminary setup: around `1633-1688`
- `ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c`
  - BLASTP one-hit/two-hit ungapped extension
- `ncbi-blast/c++/src/algo/blast/core/BLAST_aalookup.c`
  - amino-acid lookup table finalize, PV array, compressed lookup behavior
- `ncbi-blast/c++/src/algo/blast/core/blast_extend.c`
  - diagonal table lifecycle and init HSP sorting
- `ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c`
  - restricted/exact preliminary gapped alignment and traceback
- `ncbi-blast/c++/src/algo/blast/core/blast_kappa.c`
  - composition-based statistics, score normalization, redo alignment
- `ncbi-blast/c++/src/algo/blast/core/blast_traceback.c`
  - traceback batching and multithread behavior

## フェーズ 0: 測定基盤を固定する

### 目標

どの段階が遅いかを、同一入力・同一オプション・同一出力形式で測る。

### 作業

- `cargo build --release` で release binary を使う。
- LOSAT と NCBI BLASTP の比較コマンドを固定する。
- まず `outfmt 6` で formatting の影響を減らす。
- `-num_threads 1` と複数 thread の両方を測る。
- `/mnt/c/...` 上の I/O と WSL native filesystem 上の I/O を分けて測る。
- `LOSAT_TIMING=1` が不足している場合は、NCBI 処理段階に対応する timing counter を追加する。

### 測定軸

- query FASTA 読み込み時間
- subject FASTA 読み込み時間
- query encoding / SEG masking
- subject encoding
- lookup table construction
- subject scan
- ungapped extension
- preliminary gapped extension
- composition/Kappa redo
- final traceback
- output formatting/write

### 成果物

- `docs/blastp_performance_baseline.md` に、入力、コマンド、環境、elapsed/user/sys、stage timing を保存する。
- 速度改善前の output diff を保存し、以後の変更で parity regression がないことを確認する。

## フェーズ 1: BLASTP subject-level preliminary search の並列化

### 優先度

最優先。期待速度改善が最も大きく、既存の `docs/multithreading_opportunities.md` でも BLASTP の
最大の未実装改善点として扱われている。

### 方針

`LOSAT/src/algorithm/blastp/blast_engine.rs` の subject loop 本体を、global hitlist を直接更新しない
subject-local 関数へ切り出す。

推奨形:

1. worker は `(subject_index, Vec<(query_index, BlastpHspList)>)` を返す。
2. worker-local に `offset_pairs`, `diag_array`, `interval_tree`, `init_hsps`, `GapAlignScratch`,
   composition/gap scratch を持たせる。
3. main thread が `subject_index` 昇順に sort してから、既存の `BlastpHitList::update` を呼ぶ。
4. Kappa redo、final pruning、output は既存順序を維持する。

### NCBI parity チェックポイント

- subject 内の scan order は変えない。
- diagonal table の `diag_offset` 更新タイミングを NCBI と一致させる。
- two-hit state と restricted alignment retry の順序を subject-local 内で維持する。
- hitlist merge は NCBI subject order に固定する。
- thread 数を変えても output が byte-identical であることを確認する。

### 完了条件

- `-num_threads 1` と `-num_threads N` で output が一致する。
- NCBI BLASTP との比較で既存 parity を悪化させない。
- multi-subject dataset で preliminary search の elapsed time が明確に改善する。

## フェーズ 2: subject preprocessing / encoding の並列化

### 優先度

中。検索本体の並列化後に実施する。

### 方針

subject FASTA records は互いに独立しているため、protein encoding と defline/id 抽出を並列化できる。
ただし、最終的な `subjects` vector の順序は FASTA 入力順のまま固定する。

### 注意点

- NCBI BLAST database の preformatted input と完全に同じ性能特性にはならない。
- LOSAT runtime に NCBI database reader を依存として入れない。
- encoded subject cache を導入する場合は、NCBI にない user-visible behavior を作らない。内部最適化に限定する。

### 完了条件

- subject order が変わらない。
- `outfmt 6` と pairwise output が既存結果と一致する。
- large FASTA subject で preprocessing time が改善する。

## フェーズ 3: amino-acid lookup / scan hot loop の最適化

### 優先度

高。ただし parity リスクがあるため、フェーズ 1 後に profile を見てから実施する。

### 方針

NCBI の `BLAST_aalookup.c` と現在の Rust lookup table 実装を対応づけ、以下を確認する。

- protein lookup table type selection
- compressed amino-acid lookup table behavior
- PV array test のタイミング
- offset list memory layout
- query offset chain traversal
- word size 3 と blastp-fast word size 5 の分岐

改善候補:

- BLOSUM62 専用 score lookup を hot path で完全に direct table 化する。
- bounds check が残っている inner loop を source-verified な形で削減する。
- empty lookup cell の fast reject を NCBI と同じ順序で強化する。
- `offset_pairs` の reuse と capacity sizing を NCBI の longest-chain logic に合わせる。

### 完了条件

- profile 上の scan/lookup 時間が低下する。
- seed count、init HSP count、preliminary HSP count が変更前と一致する。
- NCBI との parity diff が悪化しない。

## フェーズ 4: ungapped / preliminary gapped extension の allocation 削減

### 優先度

中から高。profile で DP scratch や Vec allocation が目立つ場合に実施する。

### 方針

`aa_ungapped.c`, `blast_extend.c`, `blast_gapalign.c` の scratch lifecycle と現在の Rust 実装を比較する。

改善候補:

- per-subject `GapAlignScratch` の reuse を徹底する。
- init HSP vector の clear/reuse を subject-local scratch に寄せる。
- interval tree の rebuild 回数とタイミングを NCBI と照合する。
- restricted gapped alignment の inconclusive retry で不要な再計算を避ける。ただし NCBI と同じ retry
  条件に限定する。

### 完了条件

- allocation profile が改善する。
- preliminary gapped extension の score-only result が変更前と一致する。
- restricted/exact retry が NCBI と同じタイミングで発生する。

## フェーズ 5: composition-based statistics / Kappa redo の並列化または軽量化

### 優先度

中。ただし parity リスクが高いため、最後に近い段階で扱う。

### 方針

NCBI の `blast_kappa.c` と `blast_traceback.c` を確認し、redo alignment の batch 単位、heap 更新順序、
score normalization 順序を合わせる。

候補:

- query 単位で独立な redo alignment を並列化し、merge は query order に戻す。
- `BlastCompositionWorkspace` を thread-local にする。
- adjusted matrix construction の reuse 可能性を測る。

### 注意点

floating-point と heap admission の順序依存があり得るため、単純な parallel iterator 化は禁止する。

### 完了条件

- default `comp_based_stats=2` で output が byte-identical。
- `comp_based_stats=0`, `1`, `2`, `3` の代表ケースで regression がない。
- Kappa redo が profile 上で短縮する。

## フェーズ 6: output formatting の最適化

### 優先度

低。検索本体の parity と速度が安定してから行う。

### 方針

`outfmt 6/7` と pairwise `outfmt 0` を分けて測る。formatting が支配的なケースのみ、rendering の
precompute を検討する。

### 注意点

- final write order は NCBI order に固定する。
- pairwise alignment line wrapping、spacing、floating-point format を変えない。

### 完了条件

- output byte diff がない。
- output-heavy case でのみ measurable improvement がある。

## 検証マトリクス

| Case | 目的 |
|---|---|
| small query vs small subject | startup/overhead regression の検出 |
| many short subjects | subject-level parallelism の効果確認 |
| few long subjects | subject-local hot loop と DP cost の確認 |
| low-complexity proteins | SEG/lookup/extension explosion の確認 |
| `task=blastp` | default path |
| `task=blastp-fast` | compressed lookup/chaining path |
| `comp_based_stats=0` | Kappa なしの検索本体測定 |
| `comp_based_stats=2` | default composition adjustment |
| `outfmt=6` | formatting 影響を抑えた測定 |
| `outfmt=0` | pairwise formatting/traceback 表示測定 |

## 推奨作業順

1. baseline timing を保存する。
2. BLASTP subject loop を subject-local result 関数へ切り出す。
3. `-num_threads 1` の output が切り出し前と同じことを確認する。
4. Rayon subject-level parallel path を追加する。
5. deterministic merge を実装し、thread count independent output を確認する。
6. profile を取り直し、次の支配的 bottleneck を選ぶ。
7. subject preprocessing, lookup hot loop, gap scratch, Kappa redo の順で改善する。

## すぐに避けること

- NCBI BLASTP を LOSAT runtime の fallback として呼ぶこと。
- output が変わる高速化。
- subject 内 scan order の変更。
- channel receive order に依存した hitlist update。
- Kappa redo の単純な並列化。
- NCBI ソース行番号コメントなしの Rust 実装変更。

## 成功基準

- 既存 parity が悪化しない。
- `-num_threads N` が BLASTP 検索本体に実効的に効く。
- representative multi-subject benchmark で LOSATP elapsed time が明確に改善する。
- speedup の根拠が stage timing と profile で説明できる。
- 変更ごとに NCBI 参照、比較結果、残リスクが文書化されている。

## 実装メモ: 2026-05-11

### 実装済み

- `LOSAT/src/algorithm/blastp/blast_engine.rs` に `BlastpPreparedSubject` と
  `prepare_blastp_subjects_preserving_order` を追加した。
- subject preprocessing で subject ID、title、NCBISTDAA protein encoding をまとめて準備する。
- `parallel` feature が有効で、`-num_threads` の実効値が 2 以上、かつ subject が複数ある場合だけ
  Rayon thread pool で subject records を並列処理する。
- `par_iter().map(...).collect()` の結果を FASTA 入力順の `Vec` として受け取り、その順に
  `subject_ids`、`subject_titles`、`subjects` へ分解する。これにより NCBI subject oid order と
  downstream hitlist/output order を変えない。
- serial path は従来通り `subject_records.iter()` で同じ準備処理を行う。

### NCBI 参照コメント

- 追加コード直上に以下の NCBI 参照コメントを置いた。
  - `ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1426`
    - subject sequence iteration と `BlastSeqSrcGetSequence`。
  - `ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1659-1668`
    - multithread search で thread-local setup/diagnostics を使う考え方。

### 検証

- `rustfmt src/algorithm/blastp/blast_engine.rs --edition 2021`
- `cargo build --release`
- LOSAT-only thread determinism check:
  - `target/release/LOSAT blastp -q tests/fasta/CoBV.faa -s tests/fasta/PajaWSV.faa -n 1 --outfmt 6`
  - `target/release/LOSAT blastp -q tests/fasta/CoBV.faa -s tests/fasta/PajaWSV.faa -n 2 --outfmt 6`
  - `cmp -s` で byte-identical、両方 624 行。

### 未実施 / 残リスク

- NCBI BLASTP との比較 oracle は今回実行していない。
- full `cargo fmt` は repository 全体の `rustfmt` が長時間走り続けたため中断し、変更ファイルのみ
  `rustfmt` した。
- subject preprocessing の効果は subject FASTA の record 数と長さに依存する。次回 timing では
  `LOSAT_TIMING=1` の `blastp_subject_encoding` を、`-num_threads 1` と複数 thread の両方で比較する。

