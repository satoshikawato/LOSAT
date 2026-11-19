# v0.2.0 release test set E: the `-db` rows re-run against NCBI BLAST+ 2.17.0

Set E of [the release test plan](../../../release/v0.2.0_certification_plan.md) (section 2) on the
candidate commit `C` = `93ca2f652617acac09d416c28c293cc8a1a501d1`, 2026-11-19.

- Where: the Linux x86_64 workstation (Intel Xeon Gold 6338, Linux 6.8), CPUs booked untimed on NUMA
  node 0. NCBI BLAST+ 2.17.0 is the official Linux x64 archive; every binary has the SHA-256 of the
  oracle on the release machine (`ncbi_bin.sha256`, compared file by file on 2026-11-19).
- LOSAT: `cargo build --release --locked` of `C` with Rust 1.92.0 (`toolchain.txt`, `losat_bin.sha256`).
- Rows: `LOSAT/tests/db_search_manifest/` (eight programs: blastn, blastp, blastx, tblastn, tblastx,
  psiblast, rpsblast, rpstblastn), NCBI at `-num_threads 1`.

| Step | Command | Result |
| --- | --- | --- |
| NCBI re-run into a scratch checkout of `C` | `python3 LOSAT/tests/db_search_fixtures.py freeze --ncbi-bin <NCBI bin>` | exit 0 (379 s) |
| Frozen files unchanged | `git status --porcelain` in the scratch checkout | empty: NCBI writes every frozen file again byte for byte |
| LOSAT against the frozen bytes | `python3 LOSAT/tests/db_search_fixtures.py check --losat <losat>` | `# rows=6690 differing=0 threads=1` |

Files: `freeze.log.gz` and `check.log.gz` (one line per row), `git_status_porcelain.txt` (empty),
`ncbi_bin.sha256`, `losat_bin.sha256`, `toolchain.txt`. `MANIFEST.tsv` lists every file with its size and
SHA-256; its SHA-256 is `local_evidence.E.evidence_sha256` in `docs/release/v0.2.0_rc_contract.json`.
