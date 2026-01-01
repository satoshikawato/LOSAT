# NEXT SESSION PROMPT — 2025-12-31

# Instruction: Strict Algorithmic Fidelity Protocol

**Role:** You are an elite Systems Engineer tasked with correcting the LOSAT codebase to match the NCBI BLAST implementation perfectly.

**The Golden Rule:**
The NCBI C/C++ implementation is the **GROUND TRUTH**. Your goal is to make the Rust code behave **identically** to the NCBI code in terms of logic, side effects, and performance characteristics.

**Operational Directives (Non-Negotiable):**
1. **Zero Tolerance for Simplification:**
   * Never simplify the logic because it is "complex" or "hard to read." Complexity in BLAST is often a deliberate optimization.
   * If the NCBI code uses a convoluted state machine or weird buffer management, **you must implement equivalent logic in Rust.**
   * **Prohibited Excuses:** "It requires major refactoring," "The current code mostly works," "It's too complex." (Violation results in immediate termination of the session).

2. **Allowed Deviations (Performance Only):**
   * You are authorized to change the structure of the *current* Rust code (LOSAT) drastically if it differs from NCBI.
   * You may adapt the implementation to Rust's memory model *only* if a direct C-port would severely degrade performance (e.g., fighting the borrow checker causing excessive cloning).
   * **Condition:** The output and algorithmic complexity (Big O) must remain exactly the same.

3. **Methodology:**
   * Do not guess. Do not try to be "original."
   * Compare the NCBI source and the LOSAT source.
   * If a discrepancy is found, **fix it immediately**. Tear down existing structs or functions if necessary.
   * Prioritize correctness and speed over code aesthetics.

**Action:**
Review the provided NCBI source (Target) and the current LOSAT source (Draft). Rewrite the LOSAT code to strictly match the NCBI logic. If major structural changes are needed to achieve parity, **DO IT.**


## TBLASTX: “LOSATのヒットが微妙に少ない/多すぎる”を NCBI パリティで安定化する

### いまの状況（確定事実）
- `NZ_CP006932 self` で、BLAST+ 62,053 hits に対して LOSAT baseline は 39,522 hits と少ない。
- 不足が最大の帯域は **pident 40–45% × length 16–32 aa**。
- 代表欠損HSPの少なくとも1つは **`--include-stop-seeds` をONにすると復活**する（stop絡みseed除外が寄与）。
- ただし `--max-hits-per-kmer` を上げると **LOSATがBLAST+より“多すぎる”**方向にも崩れる（例: `--max-hits-per-kmer 5000` で 73,204 行）。
- overlong tail（異常に長い高一致HSP）の再発は、現時点の実験出力で **安全柵チェックOK**。

関連サマリ:
- `LOSAT/docs/SESSION_SUMMARY_2025-12-31_TBLASTX_HIT_PARITY_FLAGS.md`

---

### ゴール（次セッション）
1. **high-frequency word suppression を “NCBIパリティ” に収束**させる  
   - “少ない”→“多すぎる”の両端に振れるのを止め、安定して BLAST+ に近づける
2. `NZ_CP006932 self` の worst bin（40–45% × 16–32 aa）を改善しつつ、**overlong tail が再発しない**ことを必ず維持する

---

### 次にやること（作業順）
#### 1) まず差分を再確認（baseline / parity / hf）
比較対象:
- BLAST+: `LOSAT/tests/blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out`
- LOSAT baseline: `LOSAT/tests/losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.newbaseline.out`
- stop seeds: `LOSAT/tests/losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.stopseeds.out`
- stop seeds + stopstop ncbi: `LOSAT/tests/losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.stopseeds.stopstopncbi.out`
- hf5000: `LOSAT/tests/losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.stopseeds.stopstopncbi.hf5000.out`

コマンド例（repo rootから）:
```bash
cd LOSAT/tests
python3 analyze_tblastx_hit_gap.py --blast ./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out --losat ./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.newbaseline.out
python3 analyze_tblastx_hit_gap.py --blast ./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out --losat ./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.stopseeds.stopstopncbi.hf5000.out
python3 check_tblastx_overlong_tail.py ./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.stopseeds.stopstopncbi.hf5000.out
```

#### 2) “多すぎる”の内訳を把握（hf5000の膨張源を特定）
- hf5000 で増えたのがどの帯域か（pident/length/bitscore）を確認  
  - 特に `bitscore 16–32` が膨れていないか（閾値周りの疑い）
- 可能なら追加する分析:
  - 「BLAST+にはないがLOSATに大量にある」HSPを抽出するスクリプト（逆差分）

#### 3) NCBIの high-frequency suppression 実装を照合して仕様を決める
参照:
- NCBI: `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c`
観点:
- 抑制が **どの段階**（exact match grouping / neighboring 生成 / finalize）で効いているか
- `BLASTAA_SIZE=28` / compressed alphabet など、LOSATが簡略化している部分の影響

実装案（次セッションで決める）:
- (a) `--max-hits-per-kmer` の意味を NCBI と同じ場所に寄せる（適用タイミング修正）
- (b) NCBI と同じ “圧縮alphabet + scansub” を実装（手間大だが必須）

#### 4) 代表欠損HSPでの再現を維持しながら改善を評価
- `LOSAT/tests/find_tblastx_missing_hsp.py` で選んだ代表欠損HSPが、\n  suppression 調整で「復活する/しない」を確認

---

### 実験に使う主要フラグ（現状）
- `--include-stop-seeds`（default OFF）
- `--ncbi-stop-stop-score`（default OFF）
- `--max-hits-per-kmer N`（default 1000）
- `LOSAT_DIAGNOSTICS=1`（段階カウンタ）

---

### 注意（安全柵）
- パリティ寄せのフラグをONにするたびに:
  - `python3 LOSAT/tests/check_tblastx_overlong_tail.py <out>` を必ず実行
  - `max(length)` と “高一致×長大HSP” を確認してから次に進む


