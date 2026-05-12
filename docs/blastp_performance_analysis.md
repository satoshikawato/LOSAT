# LOSATP が BLASTP より遅い原因と改善策

## 結論

`LOSAT/tests/plots/execution_time_comparison_all.png` の BLASTP パネルでは、
WSSV/PajaWSV などの小さいウイルス同士では LOSATP が BLASTP と同等かやや速い一方、
AP027xxx/NZ_CP006932 系の大きいタンパク質セットでは LOSATP が数秒から十数秒遅い。

この差は起動や出力だけでは説明しにくく、主因は BLASTP の探索本体、特に
候補 HSP が増えた後の gapped alignment と composition-based statistics/Kappa redo の
データ移動、再アラインメント、HSP リスト処理にある。

比較スクリプト `LOSAT/tests/run_comparison.sh` では LOSATP と BLASTP のどちらも
`-num_threads 1` 相当で実行されているため、今回の図に限れば「BLAST+ だけが
マルチスレッドだから速い」という説明ではない。ただし LOSAT 側の BLASTP は
`-n`/`--num_threads` を受け取るものの、BLASTP の subject-level 並列化はまだ
探索本体に入っていないため、今後大きい入力で差を縮める最優先候補になる。

## 入力条件

- 参照したプロット: `LOSAT/tests/plots/execution_time_comparison_all.png`
- 計測コマンド: `LOSAT/tests/run_comparison.sh`
- LOSATP: `../target/release/LOSAT blastp ... -n 1 --outfmt 6`
- BLASTP: `blastp ... -num_threads 1 -outfmt 6`
- 実装入口: `LOSAT/src/algorithm/blastp/blast_engine.rs`
- NCBI 参照入口:
  `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1498`

## 遅い原因

### 1. BLASTP の subject ループがまだ並列化されていない

`docs/multithreading_opportunities.md` でも整理されている通り、LOSAT の BLASTN と
TBLASTX には subject-level 並列化が入っているが、BLASTP の subject-level
preliminary search は未実装のまま残っている。

NCBI の探索エンジンは subject を順に取り出し、subject ごとに search core を呼んで
HSP stream に流す。

NCBI reference:

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1498`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1633-1688`

今回のプロットは `-n 1`/`-num_threads 1` なので直接の差分要因ではないが、
AP027xxx のような subject 数と候補数が多いケースでは、LOSAT 側だけ
`-n 8` にしても BLASTP 本体の重い部分が並列化されない可能性が高い。
これは実運用の速度差を広げる構造的な問題。

### 2. Gapped alignment と traceback の Rust 実装が C のホットループより重い

LOSAT の BLASTP は NCBI の `blast_gapalign.c` に合わせて
score-only preliminary gapped alignment と traceback を Rust で実装している。
AP027xxx 系では候補 HSP が多く、ここが時間を支配しやすい。

NCBI reference:

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4209-4319`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4549-4715`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:797-931`

LOSAT 側で重くなりやすい点:

- DP/traceback scratch が `Vec`/`Vec<Vec<u8>>` ベースで、NCBI の再利用メモリ構造より
  allocator と bounds/branch の影響を受けやすい。
- score 計算がセルごとに `ScoringMatrix` 分岐や関数呼び出しを通る箇所があり、
  NCBI の `matrix_row[*b_ptr]` 型の内側ループより命令数が増えやすい。
- traceback の edit script を `Vec<GapEditOp>` として保持し、後段で clone/copy する
  場面がある。

出力を維持するには DP の範囲、x-drop、tie-break、edit script の意味は変えられない。
改善は「同じ演算をより NCBI に近いメモリ配置と行列アクセスで実行する」方向に限定する。

### 3. Composition-based statistics/Kappa redo のデータ移動が大きい

BLASTP のデフォルトは composition-based statistics mode 2
(`CompositionMatrixAdjust`) で、LOSAT も NCBI に合わせて Kappa redo を実行している。
これはパリティ上は正しいが、候補 HSP が増える大きいケースでは、redo alignment の
範囲構築、subject SEG、adjusted matrix、再 traceback が大きな負荷になる。

NCBI reference:

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2389-2394`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3675-3716`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1510-1704`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c`

LOSAT 側では `BlastCompoAlignment` の boxed linked list、`BlastCompoSequenceData` の範囲
コピー、`GapEditOp` の clone、subject range ごとの SEG 処理が見える。
NCBI に合わせる必要があるため CBS/Kappa を省略してはいけないが、scratch と range
buffer を再利用し、最終 HSP に残すまで edit script の clone を遅延させる余地がある。

### 4. HSP list pruning/common endpoint purge に O(n^2) になりやすい処理がある

`LOSAT/src/algorithm/blastp/hsp.rs` の common endpoint purge は NCBI の
`Blast_HSPListPurgeHSPsWithCommonEndpoints` を移植しているが、Rust 側では
`Vec::remove` を使う箇所がある。これは 1 件削除ごとに後続要素を詰めるため、
同じ query/subject endpoint に候補が集中するケースでは O(n^2) になりやすい。

NCBI reference:

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2455-2535`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2392-2452`

NCBI 側は pointer array を sort し、削除対象を解放または compaction する設計なので、
Rust 側も「削除マーク後に一括 compaction」または「NCBI の tail-swap compaction 順序を
明示的に再現」に寄せるべき。単純な `retain` は順序を変えない点では安全だが、
NCBI の tie-break と一時的な trimmed HSP 追加順序を壊さない確認が必要。

### 5. スコア行列アクセスが NCBI の内側ループほど特殊化されていない

NCBI の gapped alignment は query residue ごとに matrix row を取り、subject 側の
residue で直接参照する。

NCBI reference:

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:737-957`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3409-3438`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:804-883`

LOSAT の BLASTP は BLOSUM62 fast path を持っているが、ungapped extension、DP、
identity/positive 再計算などの複数箇所で matrix 判定を通る。デフォルト BLASTP は
BLOSUM62/11/1 が大半なので、BLOSUM62 専用の row pointer 風 API を hot loop に渡すと
出力を変えずに命令数を減らせる。

### 6. `blastp-fast` 用の compressed lookup は今後の注意点

`LOSAT/src/algorithm/blastp/args.rs` は NCBI と同様、word size が 5 以上なら
`CompressedAaLookupTable` を選ぶよう解決している。

NCBI reference:

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_prot_options.cpp:55-83`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:275-285`

ただし現在の実測プロットは default `blastp`、word size 3 のため、この点は今回の
遅さの主因ではない。将来 `blastp-fast` を計測する場合は、args だけでなく探索本体が
NCBI の compressed AA lookup と同じ経路を使っているかを先に確認する。

## 改善策

### 汎用性能方針: sparse duplicate を速くする

BLASTP の大きい入力では、HSP 数そのものは増えるが common endpoint duplicate が常に
高密度になるとは限らない。したがって、特定の重複集中ケースを速くする実装でも、
重複が少ない多数 query/subject ケースで固定費を増やすなら汎用性能としては失敗になる。

今後の最優先方針:

- duplicate がない、または少ないケースの固定費を増やさない。
- allocation は duplicate 検出後に遅延する。
- compaction は削除対象が実際に存在する場合だけ行う。
- 大規模 multi-FASTA では「1 subject あたりの HSP 分布」を見る。全体 HSP 数だけで
  最適化を決めない。
- 高 duplicate ケースの改善は、sparse duplicate ケースの regression gate を通ったものだけ
  残す。

### 優先度 1: BLASTP 専用 timing/counter を汎用性能ゲートとして整備

実装を変える前に、次を `LOSAT_TIMING=1` などで出す。

- FASTA read/encoding
- lookup build
- subject scan hits
- ungapped extension 件数
- preliminary gapped alignment 件数と時間
- Kappa redo 件数と時間
- HSP purge/pruning 時間
- final formatting 時間
- common endpoint purge の入力 HSP 数、削除数、trimmed 数、duplicate group 数
- duplicate-free subject 数、duplicate-sparse subject 数、duplicate-dense subject 数
- purge で allocation/compaction が発生した回数

BLAST+ 実行は比較 oracle としてのみ使い、LOSAT の runtime path からは呼ばない。
この計測で「gapped DP が支配的か」「Kappa redo が支配的か」「HSP purge が詰まっているか」
だけでなく、「purge 最適化の固定費が duplicate-free ケースで増えていないか」を見る。

### 優先度 2: HSP purge を lazy/adaptive compaction に直す

`Vec::remove` の逐次削除は高 duplicate ケースで O(n^2) になりうるため避ける。ただし、
現在の `removed` bitmap + 一括 compaction は、duplicate がないケースでも bitmap 確保と
active prefix scan の固定費を増やす。このままだと query/subject が増えて duplicate が
疎な入力で不利になりうる。

次に行うべき実装:

1. sorted HSP を走査し、duplicate を見つけるまでは `removed` bitmap も `removed_indices`
   も確保しない。
2. 最初の duplicate を見つけた時点で `removed_indices: Vec<usize>` を作り、以後は削除対象
   index だけを記録する。
3. `removed_indices.is_empty()` なら compaction せず、ソート済み `hits` をそのまま返す。
4. 削除対象がある場合だけ一括 compaction する。
5. duplicate が極端に多い場合だけ bitmap に切り替える閾値を検討する。目安は
   `removed_indices.len() > hits.len() / 8` だが、採用前に計測で決める。

期待する計算量:

- duplicate-free: sort + single scan、追加 allocation なし、compaction なし。
- duplicate-sparse: sort + scan + `removed_indices` 小容量 allocation + compaction 1 回。
- duplicate-dense: sort + scan + bitmap または indices + compaction 1 回。

注意点:

- comparator は `s_QueryOffsetCompareHSPs` と `s_QueryEndCompareHSPs` の順序を維持する。
- `purge=false` 時に切り出した trimmed HSP の追加順序を変えない。
- HSP 数だけでなく、bit score、E-value、座標、出力順の diff を必ず確認する。
- no-duplicate fixture を追加し、purge 後の vector identity/order と allocation-free path を
  確認する。

### 優先度 3: Gapped alignment scratch と scoring row の特殊化

出力を変えずにできる改善:

- BLOSUM62/11/1 の hot path では `matrix == Blosum62` 判定を外側で済ませ、row slice を
  DP 内側ループに渡す。
- `GapAlignScratch` の trace row を arena/flat buffer 化し、`Vec<Vec<u8>>` の row ごとの
  clear/resize を減らす。
- score-only preliminary alignment と traceback alignment で同じ scratch policy を共有する。
- x-drop、subject range adjustment、fence sentry、terminal gap pruning の順序は NCBI
  reference のまま維持する。
- scratch reuse は subject/query が増えても allocator pressure を下げるため、汎用性能への
  悪影響が出にくい。ただし巨大 buffer の保持でメモリ常駐量が増えないよう上限を設ける。

### 優先度 4: Kappa redo の clone/copy 削減

CBS/Kappa は削れないので、同じ処理を軽くする。

- `BlastCompoSequenceData` の range copy buffer を query/subject ごとに再利用する。
- `BlastCompoAlignment` の linked list allocation を arena か reusable node pool に寄せる。
- `GapEditOp` は最終的に残る HSP だけ clone し、redo 中は借用または arena 所有にする。
- subject SEG は NCBI の適用条件を維持しつつ、同じ subject range に対する重複マスクを避ける。
- clone/copy 削減は high-hit ケースに効くが、低 hit ケースで arena 初期化や大きい scratch
  保持が固定費にならないよう lazy 初期化する。

### 優先度 5: BLASTP subject-level 並列化

探索本体の subject loop を関数化し、worker が subject ごとの provisional HSPList を返す
形にする。

安全な形:

1. query lookup、contexts、cutoffs、search space、encoded subjects は read-only 共有。
2. worker ごとに `offset_pairs`、diag table、interval tree、gapped scratch、
   preliminary HSP list、composition scratch を持つ。
3. worker の戻り値を `(subject_index, Vec<(query_index, BlastpHspList)>)` にする。
4. main thread で `subject_index` 昇順に sort してから既存の `BlastpHitList::update` を呼ぶ。

この順序なら出力順序は現在の単一スレッド実装を維持しやすく、NCBI の subject stream
順にも近い。ただし subject 数が少ない入力や小さい subject では thread pool/scratch 作成が
固定費になるため、`num_threads > 1` でも subject 数・総 subject residues・推定 seed hit 数の
閾値を満たす場合だけ有効化する。

#### 実装状況: final HSP common-endpoint purge の一括 compaction

`LOSAT/src/algorithm/blastp/hsp.rs` の final HSP common-endpoint purge について、
`Vec::remove` による逐次削除をやめ、削除対象を `removed` bitmap にマークした後で
残存 HSP を一括 compaction する実装に置き換えた。

汎用性能上の評価:

- dense duplicate では逐次 `Vec::remove` の O(n^2) shift を避けられる。
- duplicate-free または duplicate-sparse では、`removed` bitmap allocation と compaction scan が
  固定費になる。
- query/subject 数が増え、各 subject の duplicate 密度が低い入力では、この固定費がメリットを
  上回る可能性がある。
- したがって現在の実装は最終形ではなく、次の作業で lazy/adaptive compaction に置き換える。

実装箇所:

- `compact_hits_after_endpoint_purge`
- `purge_sorted_hits_with_common_endpoints`
- `purge_hits_for_subject_ex`

NCBI reference:

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478-2502`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2504-2529`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2493-2500`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2520-2527`

保持した挙動:

- `s_QueryOffsetCompareHSPs` 相当の開始端点ソート後に、同一 context/query start/subject
  start の後続 HSP を削除対象にする。
- `s_QueryEndCompareHSPs` 相当の終了端点ソート後に、同一 context/query end/subject end
  の後続 HSP を削除対象にする。
- `purge=false` では NCBI の `s_CutOffGapEditScript` 相当処理を通して trimmed HSP を
  `extra_start` 以降へ追加する。
- BLASTP の subject は plus-only なので、NCBI の `subject.frame` 比較は `Hit` 側では
  canonical value として扱う。

検証:

- `rustfmt --edition 2021 src/algorithm/blastp/hsp.rs`: pass
- `cargo build --release`: pass
- `cargo build`: pass
- `cargo test purge_hits_with_common_endpoints --lib`: pass

#### 2026-05-12 追記: final HSP common-endpoint purge の実コード反映

上記方針を `LOSAT/src/algorithm/blastp/hsp.rs` に反映し、final HSP
common-endpoint purge の開始端点パス、終了端点パスから `Vec::remove(j)` による逐次削除を
除去した。実装は、NCBI の active prefix と tail への pointer compaction を Rust 側では
`removed` bitmap と一括 prefix compaction で再現する。

ただし汎用性能を重視するなら、この版はまだ固定費が大きい。duplicate を検出する前から
bitmap を確保し、削除対象がない場合でも compaction helper を通るため、duplicate-free な
大規模 multi-query/multi-subject では regress しうる。次の変更では、duplicate が見つかるまで
allocation しない lazy path にする。

追加・変更した実装箇所:

- `compact_hits_after_endpoint_purge`
- `purge_sorted_hits_with_common_endpoints`
- `purge_hits_for_subject_ex`

NCBI reference:

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478-2502`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2504-2527`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2493-2500`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2520-2525`

保持した挙動:

- `s_QueryOffsetCompareHSPs` 相当の開始端点ソート後、同一 context/query start/subject
  start の後続 HSP を削除対象にする。
- `s_QueryEndCompareHSPs` 相当の終了端点ソート後、同一 context/query end/subject end の
  後続 HSP を削除対象にする。
- `purge=false` の trimmed HSP は、既存 LOSAT API と同じく active HSP prefix の後ろに
  追加する。
- BLASTP runtime では NCBI の `purge |= (program != eBlastTypeBlastn)` により
  final purge は削除モード相当なので、trimmed HSP path は主に shared helper の互換性用。

検証:

- `rustfmt --edition 2021 src/algorithm/blastp/hsp.rs`: pass
- `cargo test purge_hits_with_common_endpoints --lib`: pass
- `rg -n "hits\\.remove\\(|remove\\(j\\)" LOSAT/src/algorithm/blastp/hsp.rs`: no matches

次に直すべき点（2026-05-12 lazy/adaptive 化前の未解決メモ。下記追記で
`purge_sorted_hits_with_common_endpoints` の lazy allocation と 3 系統の unit test は対応済み）:

- `purge_sorted_hits_with_common_endpoints` を lazy allocation にする。
- duplicate-free の場合は compaction せず、ソート済み `hits` をそのまま返す。
- sparse duplicate は `removed_indices`、dense duplicate は bitmap、という adaptive path を
  計測結果に基づいて採用する。
- unit test に no-duplicate / sparse-duplicate / dense-duplicate の 3 ケースを追加する。
- `LOSAT_TIMING=1` で duplicate-free subject 数と compaction 発生数を出す。

未解決:

- `LOSAT/tests/run_blastp_comparison.sh` は `LOSAT/tests` から実行すると
  `AP027078.AP027131.blastp` で 1 hit 差分が出る。
- 差分は `BDU67827.1 -> BDV02650.1` の high-scoring hit が現在出力から落ちる形。
- この public final purge 関数は現行の tabular streaming path では使われていないため、
  差分は既存 BLASTP divergence と見ているが、受け入れ条件上は未達。

#### 2026-05-12 追記: final HSP common-endpoint purge の lazy/adaptive 化

`LOSAT/src/algorithm/blastp/hsp.rs` の final HSP common-endpoint purge を、
duplicate 検出まで allocation しない lazy path に置き換えた。前回の bitmap 一括
compaction 版は duplicate-free でも `removed` bitmap を確保していたが、今回の実装では
duplicate-free なら removal tracking も compaction も行わず、ソート済み `hits` をそのまま
次の pass へ渡す。

実装内容:

- `EndpointRemovalMarks` を追加し、最初の duplicate を見つけた時点でのみ
  `Indices(Vec<usize>)` を確保する。
- 削除対象が dense 化した場合だけ `Bitmap(Vec<bool>)` に切り替える。
  閾値は `hits.len() >= 256 && removed_count > hits.len() / 8`。
- compaction は削除対象がある場合だけ 1 回実行する。
- compaction は新規 active `Vec` を作らず、既存 `Vec<Hit>` に対する `retain` 1 pass で
  NCBI と同じ active-prefix order を残す。
- `EndpointPurgePassStats` / `EndpointPurgeStats` を unit test 用に追加し、
  duplicate-free / sparse / dense の各 path を明示的に検証できるようにした。

2026-05-12 regression fix:

- 実測で微小な slowdown が見えたため、小さい HSP list では bitmap に切り替えないよう
  `ENDPOINT_REMOVAL_BITMAP_MIN_HITS = 256` を追加した。これにより通常の少数 HSP subject では
  `Vec<bool>` allocation と初期化を避け、sparse `Indices` path のまま処理する。
- `swap` + `truncate` の compaction は `Vec::retain` に変更した。NCBI の削除後 active-prefix
  order は維持しつつ、`Hit` 全体の swap を繰り返す固定費を削るため。
- `EndpointPurgePassStats` / `EndpointPurgeStats` の更新は test/instrumented path のみで行う。
  通常の `purge_hits_for_subject_ex` は `purge_hits_for_subject_ex_impl::<false>` を通し、
  stats 用の count/flag 更新を const generic で落とす。

#### 2026-05-12 追記: BLASTP gapped alignment / preliminary purge hot path の追加高速化

`LOSAT/src/algorithm/blastp/gapalign.rs` と
`LOSAT/src/algorithm/blastp/blast_engine.rs` に、`-n 1` でも効く hot path の固定費削減を追加した。

実装内容:

- traceback DP の trace row は、以前は row capacity 分を `resize(..., 0)` で先に埋めていた。
  NCBI の `state_struct->used += ...` と同じく、実際に書いた traceback cell だけを
  `gap_trace_row_write` で push し、必要な末尾 sentinel 分だけ `SCRIPT_GAP_IN_A` で補う。
- DP 内側ループの query/subject residue 取得は、NCBI の `matrix_row[*b_ptr]` に近づけるため、
  明示的な境界確認後に `get_unchecked` を使う。subject 側の out-of-window は従来通り sentinel
  `0` を返す。
- `prelim_edit_ops_to_gap_edit_script` は forward traceback ops 追加前に capacity を確保する。
- terminal gap pruning は `Vec::remove(0)` の反復をやめ、leading gap ops を数えて 1 回だけ
  `drain` する。
- preliminary HSP common-endpoint purge は `BlastpPreliminaryHsp: Copy` を使い、`swap` ではなく
  overwrite compaction に変更した。NCBI の left-shift compaction と同じ active-prefix order を
  残しつつ、削除対象を tail に退避する余分な書き込みを避ける。

NCBI reference:

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:560-665`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:689-726`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4681-4707`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2493-2500`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2519-2525`

検証:

- `rustfmt --edition 2021 src/algorithm/blastp/gapalign.rs src/algorithm/blastp/blast_engine.rs`: pass
- `cargo test gapalign --lib`: pass
- `cargo test prepare_preliminary_hits_for_kappa --lib`: pass

NCBI reference:

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478-2502`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2504-2527`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2488-2500`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2514-2525`

追加した regression tests:

- `test_purge_hits_with_common_endpoints_no_duplicate_allocates_no_removal_marks`
- `test_purge_hits_with_common_endpoints_sparse_duplicate_uses_index_marks`
- `test_purge_hits_with_common_endpoints_dense_duplicate_switches_to_bitmap`

検証:

- `rustfmt --edition 2021 src/algorithm/blastp/hsp.rs`: pass
- `cargo test purge_hits_with_common_endpoints --lib`: pass

残る注意点:

- この final purge helper は現行 tabular streaming path ではまだ main runtime hot path ではない。
- 次は `LOSAT_TIMING=1` の BLASTP counters に final/preliminary purge の input HSP 数、
  duplicate-free/sparse/dense 件数、allocation/compaction 発生数を出し、実測で支配的な path を
  判断する。

#### 2026-05-12 追記: Kappa redo subject SEG range cache と callback hot path 削減

ユーザー指摘の通り、final purge と gapped alignment まわりの前回変更は AP027xxx 系の
BLASTP 図に出るほど効かなかった。`LOSAT_TIMING=1` で AP027078 vs AP027131 を測ると、
支配的なのは `blastp_kappa_redo` だった。

変更前の同一条件:

- `blastp_total: 31.797s`
- `blastp_kappa_redo: 23.281s`
- `blastp_kappa_subject_seg_range: 37,395,261 bytes (calls=80,611, masked_residues=1,245,183)`

今回の実装内容:

- `LOSAT/src/algorithm/blastp/kappa.rs` に `BlastpKappaSubjectRangeCache` を追加した。
- Kappa redo の `get_range` callback で、同じ `(subject index, subject_range.begin,
  subject_range.end)` に対する subject SEG マスク結果を再利用する。
- `should_test_identical` が必要な near-identical 判定は NCBI と同じ順序で実施し、
  SEG 適用が決まった後だけキャッシュを参照する。near-identical により SEG をスキップする
  挙動は変えない。
- `postprocess_preliminary_hits` へ subject range cache を渡し、`-n 1` の sequential path では
  query/subject pair をまたいで同じ worker cache を再利用する。
- Kappa callback context の `GapAlignScratch` と subject range cache は、NCBI の
  `gapping_params->context` と同じ同期 callback lifetime に限定し、`RefCell` から
  `NonNull` raw pointer に変更した。これにより `get_range` と `redo_one_alignment` の
  80k 回級 hot path から動的 borrow check を外した。

NCBI reference:

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1614-1646`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:1188-1208`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1910-1916`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2389-2394`

実測結果（AP027078 vs AP027131, `-n 1 --outfmt 6`, release build）:

- subject SEG cache 後: `blastp_total: 26.375s`, `blastp_kappa_redo: 18.989s`
- raw pointer callback 化後: `blastp_total: 24.957s`, `blastp_kappa_redo: 17.891s`
- `blastp_kappa_subject_seg_range: 201,757 bytes (calls=595, masked_residues=6,627)`
- `blastp_kappa_subject_seg_cache: hits=80,016 misses=595`
- `/tmp/losat_ap027078_ap027131.out` との `diff -q` は差分なし

残る支配要因:

- subject SEG の重複実行はほぼ消えたが、Kappa redo 自体はまだ 17.891s と最大要素。
- 残りは 80k 件級の redo traceback、adjusted matrix 適用、linked-list alignment 処理が中心。
- 次は `blast_redo_one_match_with_workspace` 内の redo alignment 呼び出し件数、
  `blast_adjust_scores` 件数、`with_distinct_ends` の containment scan 件数を timing/counter 化し、
  実測で最も大きい path だけを削る。

## 受け入れ条件

- `LOSAT/tests/run_blastp_comparison.sh` の outfmt 6 diff がゼロ。
- `LOSAT/tests/run_blastp_pairwise_comparison.sh` も対象にする場合は pairwise 出力 diff がゼロ。
- NCBI BLAST+ は比較 oracle としてのみ実行し、LOSAT runtime/build/fallback には使わない。
- hit count 改善だけを合格条件にしない。bit score、E-value、座標、出力順まで確認する。
- 速度改善の PR でも、該当 Rust コードには NCBI C/C++ reference コメントを追加する。
- 汎用性能改善では、small viral pair、AP027xxx/NZ_CP006932、synthetic duplicate-free、
  synthetic duplicate-dense の少なくとも 4 系統を測る。
- duplicate-free/sparse ケースで HSP purge の wall time と allocation count が増える変更は
  reject する。dense duplicate の改善だけでは受け入れない。
- `-n 1` と `-n 8` の出力完全一致を維持し、thread overhead が小入力を遅くする場合は自動で
  single-thread path を選ぶ。

## 推奨作業順

1. BLASTP timing/counter を拡張し、purge 入力 HSP 数、duplicate 数、duplicate density、
   allocation/compaction 発生数を出す。
2. 現在の final HSP common-endpoint purge を lazy/adaptive compaction に改める。
   duplicate-free では allocation も compaction もしない。
3. no-duplicate / sparse-duplicate / dense-duplicate の unit test と micro fixture を追加し、
   汎用性能 regression gate を作る。
4. AP027131.NZ_CP006932 と WSSV.PajaWSV の `LOSAT_TIMING=1` 内訳を取り、HSP purge が本当に
   支配的なケースだけ purge 最適化を続ける。
5. Gapped alignment scratch と BLOSUM62 row hot path は、allocator 固定費を増やさない形で進める。
6. Kappa redo の buffer reuse と edit script clone 遅延を入れる。arena/scratch は lazy 初期化し、
   小入力で常駐メモリを増やさない。
7. subject-level parallelism は最後に入れる。subject 数・総 residues・推定 seed hit 数の閾値で
   小入力を single-thread に戻し、`-n 1` と `-n 8` の出力完全一致を確認する。

重要な判断:

- 以前の「2 と 3 を先に入れる」は dense duplicate や特定 hot loop を前提にしすぎていた。
- 汎用性能では、まず duplicate-free/sparse ケースの固定費を消す。
- その後、計測で支配的な hot path だけを順に削る。
