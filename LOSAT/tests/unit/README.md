# Unit Tests for LOSAT

## Overview

This directory contains unit tests for LOSAT modules, organized by component.

## Important Discovery: Implementation Already Ported

**✅ LOSAT already has NCBI BLAST's length adjustment implementation ported**

- **File**: `src/stats/length_adjustment.rs`
- **Function**: `compute_length_adjustment_ncbi()`
- **Reference**: NCBI BLAST `blast_stat.c`'s `BLAST_ComputeLengthAdjustment`

This means we can aim for **exact match** (or <0.1% difference) rather than just "tolerance-based comparison".

See `helpers/SUMMARY.md` for details on the verification process and tools.

## Test Structure

- `common/` - Tests for common modules (diagnostics, evalue, chaining)
- `blastn/` - Tests for BLASTN modules
- `tblastx/` - Tests for TBLASTX modules
- `helpers/` - Test utilities and helpers

## NCBI BLAST Integration

### Testing Strategy (Multi-Layered Approach)

単純にλとKを逆算するだけでは不十分です。以下の4つの問題に対処します：

#### 1. **Effective Search Space Verification** ⚠️ 最重要
- NCBI BLASTのヘッダーから実効探索空間を抽出
- LOSATの計算結果と比較
- これが最も重要な検証項目（length adjustmentのバグを検出）

#### 2. **Composition-based Statistics Avoidance**
- TBLASTX/BLASTPでは `-comp_based_stats 0` で実行
- 純粋な統計値（固定のλ, K）と比較

#### 3. **Multiple Parameter Set Testing**
- 複数のスコアリングパラメータセットでテスト
- デフォルト値だけでなく、様々な組み合わせを検証

#### 4. **Length Variation Testing**
- 短い配列（100bp）、中程度（1kbp）、長い配列（1Mbp）でテスト
- 長さ調整が正しく機能するか確認

### Current Status

**Integration Tests (Manual)**: 
- `tests/run_comparison.sh` runs both LOSAT and NCBI BLAST and compares outputs
- `tests/plot_comparison.py` visualizes differences between LOSAT and NCBI BLAST

**Unit Tests (Automated)**:
- Basic functionality tests are implemented
- **NCBI BLAST reference data comparison framework is ready**
- Reference data extraction scripts are available

### Adding NCBI BLAST Reference Data

To add NCBI BLAST reference data for unit tests:

1. **Run NCBI BLAST with proper options**:
   ```bash
   # BLASTN
   blastn -query test.fasta -subject db.fasta -outfmt 7 > ncbi_output.txt
   
   # TBLASTX (⚠️ 重要: -comp_based_stats 0 を指定)
   tblastx -query test.fasta -subject db.fasta \
           -comp_based_stats 0 \
           -outfmt 7 > ncbi_output.txt
   ```

2. **Extract test cases**:
   ```bash
   python tests/unit/helpers/extract_ncbi_cases_enhanced.py \
       ncbi_output.txt \
       <query_length> <db_length> <db_num_seqs> \
       [blastn|tblastx]
   ```

3. **Populate reference data**:
   - Copy generated Rust code to `helpers/ncbi_reference.rs`
   - Add to `get_ncbi_blastn_evalue_cases()` or `get_ncbi_tblastx_evalue_cases()`
   - **Important**: Include `expected_effective_space` from NCBI BLAST header

4. **Enable NCBI comparison tests**:
   - Remove `#[ignore]` from tests in `common/evalue_ncbi.rs`
   - Run: `cargo test --test unit_tests evalue_ncbi`

### Test Role Division

**Unit Tests** (this module):
- Purpose: Verify calculation logic correctness
- Verifies:
  - Bit score calculation accuracy
  - E-value calculation accuracy
  - **Effective search space calculation** (most important)
  - Parameter set consistency
  - Behavior with various lengths

**Integration Tests** (`tests/run_comparison.sh`):
- Purpose: Verify end-to-end NCBI BLAST compatibility
- Verifies:
  - Hit counts on real sequences
  - E-value distribution
  - Bit score distribution
  - Alignment coordinates

## Running Tests

```bash
# Run all unit tests
cargo test --test unit_tests

# Run specific test module
cargo test --test unit_tests common::diagnostics

# Run with output
cargo test --test unit_tests -- --nocapture
```

## Test Coverage Goals

- Common modules: > 80% coverage
- BLASTN modules: > 70% coverage for critical paths
- TBLASTX modules: > 70% coverage for critical paths
- NCBI BLAST compatibility: All critical calculations verified against reference data

