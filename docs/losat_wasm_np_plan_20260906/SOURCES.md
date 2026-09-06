# 参照元と検証範囲

作成時確認日: 2026-09-06。コードは明示したcommitに固定して参照する。外部runtime文書は随時更新されるため、実行時の固定依存版・runtime版で確認し直す。

総合計画と指示書のR番号は本ファイルの項目に対応する。URLは再取得用であり、参照リポジトリ全体、実行バイナリ、原始出力、NCBIソースをこのZIPへ同梱したものではない。静的読解で確認したこと、過去の保存済み実測値、今後の改善仮説を区別する。記載する関数名は探索の入口であり、変更時はcall pathと現在の行範囲まで確認する。

## R01

**LOSATのrepositoryと基準commit**

2026-09-06のGitHub branch読取でmainのSHAを確認。計画の静的読解基準であり、実行担当のcheckoutや将来のmainを固定するものではない。

```text
https://github.com/satoshikawato/LOSAT
https://github.com/satoshikawato/LOSAT/commit/7db9bb0060e4e057f9f50807bf9edc2362f20133
```

## R02

**AGENTS.md — 互換性・Native Gate A/B・開発規約**

NCBIを意味仕様の権威とすること、runtime依存禁止、変更のsourceコメント、PD-NCBI-PLATFORM-VARIANCE、Native/Wasm境界、検証方針の根拠。作業開始時は手元の現行版と局所指示を優先する。

```text
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/AGENTS.md
```

## R03

**NCBI parity/speed用のrepository skill**

SKILLとexceptions/commandsの読取を基に検証の順序と境界を整理。evidence.mdは実行時に利用する既存記録の入口。認定済みraw goldenの所在は、この計画には同梱していない。

```text
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/.agents/skills/verify-ncbi-parity-and-speed/SKILL.md
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/.agents/skills/verify-ncbi-parity-and-speed/references/exceptions.md
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/.agents/skills/verify-ncbi-parity-and-speed/references/commands.md
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/.agents/skills/verify-ncbi-parity-and-speed/references/evidence.md
```

## R04

**Cargo設定**

Cargo.tomlとconfigの読解により既存release最適化とWasm SIMD設定を確認。Cargo.lockは実行時の依存版確認の入口。設定記述だけで配布artifactへの適用を証明したとは扱わない。

```text
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/LOSAT/Cargo.toml
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/LOSAT/.cargo/config.toml
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/LOSAT/Cargo.lock
```

## R05

**v0.1.0の保存済みbenchmark**

総合計画§3.1の数値・試行数・環境・effective thread分類の出典。測定対象コードはaf3e2ea837afdb8a00cf19920f68be4f0bf3bfb5。計画作成時に新規測定をしたものではない。

```text
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/benchmarks/v0.1.0/README.md
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/benchmarks/v0.1.0/plot_data.json
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/benchmarks/v0.1.0/execution_times.tsv
```

## R06

**LOSATN greedy memory/traceback**

GreedyNonAffineMem::alloc_traceback_row、NonAffineGreedyRow/persist等に関する静的読解の根拠。copy削減の価値とproduction到達は未計測。S02で対応する実行時call pathを確定する。

```text
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/LOSAT/src/algorithm/blastn/alignment/greedy.rs
```

## R07

**LOSATN ungapped extension**

extend_hit_ungapped_approx_ncbi/exact_ncbi、4塩基pack、score_table、X-drop/exact再計算の根拠。新たな近似アルゴリズムを導入する提案ではない。

```text
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/LOSAT/src/algorithm/blastn/extension.rs
```

## R08

**LOSATN sequence comparison**

比較dispatchとfind_first_mismatch_exの読解の根拠。source上のscalar記述から生成コードやproduction実行費用を断定しない。

```text
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/LOSAT/src/algorithm/blastn/sequence_compare.rs
```

## R09

**LOSATN run coordinatorと既存並列化**

subject/chunk、MAX_DBSEQ_LEN、Wasm並列化判定、pool、入力準備とreductionの根拠。新規partitionの正当性はこのファイルの存在だけでは証明されない。

```text
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/LOSAT/src/algorithm/blastn/blast_engine/run.rs
```

## R10

**LOSATP engine・計測・preliminary/Kappa並列経路**

blastp_timing_env_enabled、get_or_build_blastp_thread_pool、prepare_independent_subject、prelimのfor_each_init、single-query match-redo、ordered heap replayの根拠。実行経路と性能比率はS01で測る。

```text
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/LOSAT/src/algorithm/blastp/blast_engine.rs
```

## R11

**LOSATP gapped alignmentとmatrix参照**

BLOSUM62特殊化、BlastpScoreMatrix、adjusted matrixのrow取得、GapAlignScratch関連の根拠。既存最適化を再実装しないための参照。

```text
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/LOSAT/src/algorithm/blastp/gapalign.rs
```

## R12

**旧改善計画・性能分析（履歴資料）**

開発上の検討事項を探す資料。未実装と記載された内容が現行コードでは実装済みの場合がある。本文の現状認定、採否、現在のタイミング値の権威にしない。

```text
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/docs/wasm_multithreading_improvement_plan.md
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/docs/blastp_wasm_current_thread_pool_plan.md
https://github.com/satoshikawato/LOSAT/blob/7db9bb0060e4e057f9f50807bf9edc2362f20133/docs/blastp_performance_analysis.md
```

## R13

**gbdraw側の既存LOSAT host（特定SHAの参考）**

module共有、thread workerの準備・起動・終了、outer/inner thread予算を確認する入口。LOSATのSHAとは独立。S10で現在のgbdraw SHA、実asset、全entrypointを再確認する。

```text
https://github.com/satoshikawato/gbdraw/blob/f9c4654bd00c4fe940974d90b85407c1d3a2d9a0/gbdraw/web/js/services/losat.js
https://github.com/satoshikawato/gbdraw/blob/f9c4654bd00c4fe940974d90b85407c1d3a2d9a0/gbdraw/web/js/workers/losat-threaded-worker.js
```

## R14

**Rayon ParallelIteratorの初期化契約**

公式API文書。job局所の初期化をOS workerごとの一回初期化と見なさないための参照。固定Cargo.lockが別版なら対応版の文書・実装で確認する。

```text
https://docs.rs/rayon/1.11.0/rayon/iter/trait.ParallelIterator.html#method.for_each_init
https://docs.rs/rayon/1.11.0/rayon/iter/trait.ParallelIterator.html#method.map_init
```

## R15

**Rayon ThreadPoolBuilderの寿命**

2026-09-06の公式文書読取ではlatestは1.12.0。local poolにuse_current_threadを使う場合のregistry注意を確認。対象repositoryが同版であるとの主張ではなく、実際のlock版でも照合が必要。

```text
https://docs.rs/rayon/latest/rayon/struct.ThreadPoolBuilder.html#method.use_current_thread
```

## R16

**Cargo configurationとRust flags**

公式Cargo文書。環境変数とtarget別rustflags等の適用元を確認する入口。性能実験では設定ファイルだけでなく実commandとartifactを記録する。

```text
https://doc.rust-lang.org/cargo/reference/config.html#buildrustflags
```

## R17

**V8 Wasm compilation pipeline**

V8の一次資料。baseline/optimizing tierとWasmでのOSRに関する説明を、初回実行と再利用時の区別に用いる。別engineへ一般化しない。

```text
https://v8.dev/docs/wasm-compilation-pipeline
```

## R18

**SharedArrayBufferとbrowser capability**

Mozillaの実装者文書。secure context、cross-origin isolation、shared memoryの確認に用いる。配信headerがあるだけで実threadingや性能が保証されるわけではない。

```text
https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/SharedArrayBuffer
https://developer.mozilla.org/en-US/docs/WebAssembly/Reference/JavaScript_interface/Memory
```

## R19

**Rust wasm32-wasip1-threads target**

Rust公式target資料。threads targetの前提を確認する入口。browserにこのtarget用の適切なhostが不要になるという意味ではない。

```text
https://doc.rust-lang.org/rustc/platform-support/wasm32-wasip1-threads.html
```

## 引用・記録の運用

実行reportでは対象コードの実在commitと関数、NCBI側の版・file・行範囲、検証command、raw evidenceを併記する。本文の説明だけで実行を証明しない。既存sourceのcommentに書かれた個人PCのabsolute pathや古い行番号を、確認せずに新しい根拠として再利用しない。

本パッケージには、実装済み改善、実行済み新規benchmark、原本へのbyte一致証明は含まれない。それらは各sessionの成果物として生成・検証する。
