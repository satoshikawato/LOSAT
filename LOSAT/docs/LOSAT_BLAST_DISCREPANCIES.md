## LOSAT BLAST Discrepancies (NCBI Parity Audit)

### 目的

`TBLASTX_NCBI_PARITY_DEVLOG.md` の **直近5セッション**で行われた修正について、
**NCBI BLAST 実装（ground truth）**: `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/`
と LOSAT 実装を突き合わせて検証し、**非準拠（出力差分につながる/仕様解釈が違う/移植が未完了/誤植や誤読）** を網羅的に列挙する。

### 運用ルール

- **NCBI 実装が唯一の正解**。ここに書くのは「NCBI と一致していない点」のみ。
- 推測は避け、可能な限り
  - **NCBI 側の該当ファイル/関数/行**
  - **LOSAT 側の該当ファイル/関数**
  - **差分の根拠（条件分岐/定数/データ構造差による論理差）**
  を明記する。
- 不確実なものは **Needs verification** として別枠に隔離し、確証が取れた時点で確定セクションへ移す。

### 対象セッション（devlogの直近5セッション）

devlog 末尾から（ファイル順）:

- `2026-01-02 追記: sum_stats_linking NCBI完全準拠修正`
- `2026-01-01 セッション: 侵入型リンクリスト実装と NCBI パリティ修正`
- `2026-01-02 追記: best 選択バグの修正と can_skip 最適化の安全化`
- `2026-01-02 セッション終了: sum_stats_linking の根本的なバグ発見`
- `2026-01-02 追記: 重大バグ発見と修正 - extension で masked sequence を使用していた`

---

## Discrepancies（確定）

### 1) devlogの主張がNCBI実装と矛盾: 「extensionはunmasked(sequence_nomask)を使う」

- **devlog**: `2026-01-02 追記: 重大バグ発見と修正 - extension で masked sequence を使用していた` にて
  - 「LOSAT は extension 時に SEG masked された `aa_seq` を使用していた。NCBI は `sequence_nomask` (unmasked) を使う」
  - という主張＋ `aa_seq_nomask` を使う修正案が書かれている。
- **NCBI ground truth**:
  - マスク設定で `query_blk->sequence_start_nomask` / `query_blk->sequence_nomask` を作るのは事実だが、これは **非マスク配列の保持**であり（例: identity 計算や出力用途）、
    extension で参照される配列は **`query->sequence`（マスク済み）**。
  - `aa_ungapped.c` の `s_BlastAaExtendTwoHit()` は `Uint1 *q = query->sequence;` を参照し、`matrix[q[q_right_off + i]]...` を評価する（`query->sequence_nomask` は参照しない）。
- **NCBI refs**:
  - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_filter.c:1350-1406`
    - `BlastSetUp_MaskQuery()` が `sequence_start_nomask` を作りつつ、マスクは `query_blk->sequence`（working buffer）に適用している。
  - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:1089-1158`
    - `Uint1 *q = query->sequence;` を用いて ungapped extension を実施。
- **LOSAT現状**:
  - `LOSAT/src/algorithm/tblastx/utils.rs` の normal/neighbor-map の two-hit extension は `aa_seq`（working buffer）を渡している。
  - これは **NCBI と整合**する可能性が高い。
- **結論**:
  - devlog の「NCBI は sequence_nomask を extension に使う」は **誤り**。
  - もし実装が devlog 記述どおり `aa_seq_nomask` を extension に使うようになっている箇所が残っていれば、それは **NCBI 非準拠**（出力差分要因）となる。

### 2) `sum_stats_linking` の link cutoffs 計算が NCBI `CalculateLinkHSPCutoffs()` と一致していない

- **対象**:
  - LOSAT: `LOSAT/src/algorithm/tblastx/sum_stats_linking.rs:calculate_link_hsp_cutoffs()`
  - NCBI: `blast_parameters.c:CalculateLinkHSPCutoffs()`（コメントでも参照されている）
- **不一致点（確定）**:
  - **`cutoff_small_gap` の下限**:
    - NCBI: `cutoff_small_gap = MAX(word_params->cutoff_score_min, floor(log(x)/Lambda)+1)`
    - LOSAT: `cutoff_small = floor(log(x)/lambda)+1`（`word_params->cutoff_score_min` 相当が反映されていない）
  - **`gap_prob` の扱い（small search space 時）**:
    - NCBI: `search_sp <= 8*window^2` の場合 `link_hsp_params->gap_prob = 0; cutoff_small_gap=0;`
    - LOSAT: `ignore_small_gaps` は判定しているが、`gap_prob` は常に定数 `GAP_PROB_UNGAPPED` を使用しており、
      `prob[1]` の「(1-gap_prob) で割る」補正が **NCBI と一致しない**（NCBIは gap_prob=0 → 割らない）。
  - **`y_variable` の計算分岐**:
    - NCBI: database search の場合 `db_length > subject_length` で `log(db_length/subject_length)` を使う分岐がある。
    - LOSAT: `y_variable = log((subject_length + expected_length)/subject_length) * K / gap_decay_rate` のみ（db_length 分岐なし）。
  - **expected length の丸め**:
    - NCBI: `expected_length = BLAST_Nint(log(K*q*s)/H)`（整数丸め）
    - LOSAT: 浮動小数のまま差し引き（`max(1.0)` でクランプ）
- **NCBI refs**:
  - `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:998-1082`

（注）この差分は **E-value 分布/最終 hit 数**に直結するため、Strict Algorithmic Fidelity Protocol 上「最優先で是正」対象。

---

## Needs verification（未確定）

TBD

---

## 修正計画（本ドキュメントに基づく網羅計画）

TBD


