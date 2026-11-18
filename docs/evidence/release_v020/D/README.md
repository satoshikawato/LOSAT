# v0.2.0 release test set D: the `-subject` rows re-run against NCBI BLAST+ 2.17.0

Set D of [the release test plan](../../../release/v0.2.0_certification_plan.md) (section 2) on the
candidate commit `C` = `703caed171944b0763c2efa4b3ef58805099fc47`, 2026-11-19.

- Where: the Linux x86_64 workstation (Intel Xeon Gold 6338, Linux 6.8), CPUs booked untimed on NUMA
  node 0. NCBI BLAST+ 2.17.0 is the official Linux x64 archive; every binary has the SHA-256 of the
  oracle on the release machine (`ncbi_bin.sha256`, compared file by file on 2026-11-19).
- LOSAT: `cargo build --release --locked` of `C` with Rust 1.92.0 (`toolchain.txt`, `losat_bin.sha256`).
- Rows: `LOSAT/tests/outfmt0_manifest.tsv`, 2,602 rows (blastn 635, blastp 591, blastx 529, tblastn 404,
  tblastx 348, psiblast 95).

| Step | Command (from the repository root) | Result |
| --- | --- | --- |
| NCBI re-run | `TMPDIR=/tmp python3 docs/evidence/losat_web_e2a/run_oracle.py --bin-dir <NCBI bin> --out <tmp>` | `# fixtures=2602 differing=[]`: every NCBI output has the frozen hash (1,171 s) |
| Frozen files untouched | `git status --porcelain` | empty |
| LOSAT, 1 thread | `python3 docs/evidence/losat_web_e2a/check_losat.py --losat <losat> --threads 1` | `# fixtures=2602 differing=[]` |
| LOSAT, 2 threads | same with `--threads 2` | `# fixtures=2602 differing=[]` |
| LOSAT, 4 threads | same with `--threads 4` | `# fixtures=2602 differing=[]` |

Files: `run_oracle.log.gz` and `check_losat_t{1,2,4}.log.gz` (one line per row), `git_status_porcelain.txt`
(empty), `ncbi_bin.sha256`, `losat_bin.sha256`, `toolchain.txt`. `MANIFEST.tsv` lists every file with its
size and SHA-256; its SHA-256 is `local_evidence.D.evidence_sha256` in `docs/release/v0.2.0_rc_contract.json`.
