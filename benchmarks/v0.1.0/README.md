# LOSAT v0.1.0 Benchmark Suite & Evaluation Snapshot

This directory contains benchmark data, performance metrics, and evaluation figures comparing LOSAT against official NCBI BLAST+ (v2.17.0+) across representative genomic and proteomic datasets.

---

## Overview & Key Findings

The benchmark evaluates three execution environments across wall-clock time, multi-threading scaling, and numerical output parity:
1. **NCBI BLAST+ (v2.17.0+)**: Reference baseline (single-thread `-num_threads 1` and multithreaded `-num_threads 8`).
2. **LOSAT Native (Rust)**: Standalone native binary (`-num_threads 1` and Rayon-backed `-num_threads 8`).
3. **LOSAT WASI WebAssembly**: Single-threaded serial (`wasm32-wasip1`) and multithreaded WASI (`wasm32-wasip1-threads` n8).

### Key Takeaways
- **Bit-Perfect Parity**: Across 1,371,516 verified alignment rows, LOSAT produces identical coordinates, bit scores, and E-values matching NCBI BLAST+.
- **Native Efficiency**: LOSAT matches or outperforms NCBI BLAST+ on pairwise FASTA comparisons, eliminating database indexing overhead and sequence extraction latency.
- **Near-Native WebAssembly**: The multithreaded WASI build delivers competitive throughput, enabling client-side BLAST analyses in web browsers and sandboxed environments.

---

## Representative Benchmark Cases

The benchmark suite tests six representative biological cases spanning viral and bacterial genomes:

| Case ID | Program / Task | Query vs Subject | Biological Context |
|:---|:---|:---|:---|
| `Sakai.MG1655.megablast` | `blastn` (megablast) | *E. coli* Sakai vs *E. coli* MG1655 | Large bacterial chromosome alignment (~5.5 Mb) |
| `PesePMNV.MjPMNV.task_blastn` | `blastn` (blastn) | Penaeid shrimp viral isolates | Divergent nucleotide sequence alignment |
| `pairwise_default_serial` | `blastp` | Baculovirus proteome | Pairwise protein homology search |
| `p11_avclpv_psclpv` | `tblastx` | Clarireovirus isolates | Translated 6-frame viral genome synteny |
| `p03_mela_pemojnva` | `tblastx` | Nudivirus genomes | 6-frame translated comparative genomics |
| `d06_ap027131_ap027133_db4` | `tblastx` | Bacterial genomes (gencode 4) | Translated alignment with Mycoplasma genetic code |

---

## Timing Protocol & Environment

- **Hardware**: Intel Core i9-14900HX (32 logical threads), Ubuntu 24.04.3 LTS (WSL2).
- **Software**: Rust 1.92.0, Node.js v24 LTS (WASI runtime), NCBI BLAST+ 2.17.0+.
- **Protocol**: 1 warmup run followed by 5 timed repetitions per case and execution mode. Bars represent median wall-clock time, with whiskers indicating the full min–max range across samples.
- **Modes Evaluated**: NCBI BLAST+ (n1/n8), LOSAT Native (n1/n8), LOSAT WASI Serial (n1), and LOSAT WASI Threaded (n8).

---

## Data Completeness & Integrity

| Dataset | Status | Scope & Interpretation |
|:---|:---|:---|
| **Alignment Output** | `AVAILABLE_EXACT` | 1,371,516 normalized rows from 43 retained native contracts; exact hashes recorded in `metadata.json`. Matches Linux x86_64 runtime and NCBI BLAST+ 2.17.0+. |
| **Current Timing Samples** | `AVAILABLE_CURRENT` | 180 timed wall-clock samples across 6 cases, 6 execution modes, and 5 repetitions. |
| **Historical Comparison** | `AVAILABLE_HISTORICAL` | Preserved historical timing data points for longitudinal performance tracking. |

---

## Directory Contents

- `metadata.json`: Dataset provenance, commit SHAs, sequence checksums, compiler flags, and execution commands.
- `alignment_results.tsv.gz`: Normalized alignment rows (`program`, `case_id`, `qseqid`, `sseqid`, `pident`, `length`, etc.).
- `execution_times.tsv`: Individual timed repetitions with platform provenance and effective-thread records.
- `plot_data.json`: Deterministic summary statistics (median, mean, min, max) and histogram bins for CI verification.
- `hit_distribution.png`: Visual alignment distribution across percent identity and match length.
- `execution_time.png`: 2×2 faceted execution time comparison (TBLASTX, MegaBLAST, BLASTN, BLASTP).
- `render_manifest.json`: Checksums for reproducible figure rendering.

---

## Reproducing Figures

To re-render the benchmark figures using the deterministic data snapshot (requires Python 3, NumPy >= 2.2.5, and Matplotlib >= 3.10.3):

```bash
python scripts/render_benchmark_plots.py \
  --snapshot benchmarks/v0.1.0 \
  --output benchmarks/v0.1.0
```

The rendering script operates strictly on the local frozen snapshot data without requiring network access, binary builds, or re-running sequence searches.
