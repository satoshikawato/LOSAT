# ファイル分割計画

## 1. tblastx/utils.rs (3645行) → 6ファイルに分割

### 分割後の構造:
```
tblastx/utils/
├── mod.rs          # モジュール定義と公開API
├── types.rs        # 構造体定義 (DiagStruct, OffsetPair, InitHSP, WorkerState) + SIMDヘルパー
├── tracing.rs      # トレース/デバッグ関連 (trace_hsp_target, trace_match_target等)
├── scan.rs         # スキャン関数 (s_blast_aa_scan_subject, s_blast_aa_scan_subject_lazy)
├── hsp_processing.rs # HSP処理 (adjust_initial_hsp_offsets, get_ungapped_hsp_list, reevaluate_ungapped_hsp_list, purge_hsps_with_common_endpoints)
├── run.rs          # メイン run() 関数 (標準モード)
└── run_neighbor_map.rs # run_with_neighbor_map() 関数
```

### 推定行数:
- types.rs: ~330行 (構造体定義 + SIMDヘルパー)
- tracing.rs: ~150行
- scan.rs: ~500行
- hsp_processing.rs: ~600行
- run.rs: ~1100行
- run_neighbor_map.rs: ~1100行

## 2. tblastx/sum_stats_linking.rs (2034行) → 3ファイルに分割

### 分割後の構造:
```
tblastx/sum_stats_linking/
├── mod.rs          # モジュール定義と公開API
├── types.rs        # 構造体 (LinkingParams, LhHelper, HspLink, BufferPools)
├── cutoffs.rs      # カットオフ計算 (calculate_link_hsp_cutoffs_ncbi等)
└── linking.rs      # リンキングロジック (apply_sum_stats_even_gap_linking, link_hsp_group_ncbi)
```

### 推定行数:
- types.rs: ~150行
- cutoffs.rs: ~200行
- linking.rs: ~1600行

## 3. blastn/utils.rs (1802行) → 3ファイルに分割

### 分割後の構造:
```
blastn/utils/
├── mod.rs          # モジュール定義と公開API
├── runner.rs       # メインrun()関数
├── filtering.rs    # filter_hsps関数
└── search.rs       # 検索ロジック (ワーカー処理部分)
```

### 推定行数:
- runner.rs: ~800行
- filtering.rs: ~200行
- search.rs: ~800行

## 4. tblastx/lookup.rs (1120行) → 2ファイルに分割

### 分割後の構造:
```
tblastx/lookup/
├── mod.rs          # モジュール定義と公開API
├── table.rs        # BlastAaLookupTable, BackboneCell, build_ncbi_lookup
└── neighbor.rs     # NeighborMap, NeighborLookup, build_direct_lookup
```

### 推定行数:
- table.rs: ~500行
- neighbor.rs: ~600行

## 実装のポイント

1. **mod.rsの活用**: 各サブモジュールに`mod.rs`を作成し、`pub use`で既存APIを維持
2. **依存関係の整理**: `types.rs`を先に分離し、他のモジュールから参照
3. **テストの移動**: インラインテストがある場合は対応するファイルに移動

## 実装順序

1. tblastx/utils.rs (最大のファイル)
2. tblastx/sum_stats_linking.rs
3. blastn/utils.rs
4. tblastx/lookup.rs

