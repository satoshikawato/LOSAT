# Session F — 並列・Wasm

## INSTRUCTION PROMPT

`feature/tlosan-tblastn-v0.2.0` で TLOSAN v0.2.0 の段階 F を実行する。
開始時にブランチと作業ツリー、[総合計画書](../tlosan_tblastn_v0.2.0_plan.md)、
`AGENTS.md`、`verify-ncbi-parity-and-speed`、段階 A〜E の証拠を確認する。
段階 E の serial 出力契約が未確定なら並列化を始めない。

1. NCBI の入力順、フレーム内・フレーム間 HSP 統合、linking、
   Kappa redo、hitlist の比較・削除順を保持して独立した query/subject 仕事を
   分配し、NCBI 順に集約する。候補集合や演算順が変わる近似を入れない。
   Rust の移植箇所の直上に NCBI ファイル・行番号と C/C++ 断片を記す。
2. 複数の独立仕事がある固定 fixture で、Native
   `-num_threads 1,2,4,8` の 0/6/7 出力 SHA-256 を比較し、
   同点と複数 HSP のケースは繰り返して決定性を検証する。
   並列経路が実際に動いた計測記録を残す。
3. plain `wasm32-wasip1` は serial として扱う。
   配布対象に含める場合、実スレッドは `wasm32-wasip1-threads` と
   `wasm-threads` の組合せで検証し、同じ fixture の Native 出力と比較する。
   使用不能な target は理由を記録して認証範囲から除く。
4. NCBI ローカル `-subject` の `-num_threads 4` を
   NCBI の並列性能証拠として扱わない。コード 1 の NCBI 生バイト契約と
   非標準コードの限定契約を全スレッドで維持する。

成果物は target・thread ごとのコマンドと SHA-256、実際の並列稼働証拠、
繰り返し実行の結果、未認証 target の明示である。段階 G に固定 fixture と
出力 checksum を渡す。
