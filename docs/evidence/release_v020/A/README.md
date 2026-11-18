# v0.2.0 release test set A: Gate A v0.1.1, Linux x64 native and serial Wasm

Set A of [the release test plan](../../../release/v0.2.0_certification_plan.md) (section 2) on the
candidate commit `C` = `703caed171944b0763c2efa4b3ef58805099fc47`, 2026-11-19, 04:45-06:03 JST (4,656 s).

- Where: the Linux x86_64 workstation (Intel Xeon Gold 6338, Linux 6.8), two CPUs booked untimed on NUMA
  node 0. NCBI BLAST+ 2.17.0 is the official Linux x64 archive; every binary has the SHA-256 of the
  oracle on the release machine (`ncbi_bin.sha256`; the script checks blastn, blastp and tblastx again,
  `oracle_identities.json`).
- LOSAT: `cargo build --release --locked` and the serial command (`cargo build --release --bin losat
  --target wasm32-wasip1 --target-dir target/serial-command --no-default-features --locked`) of `C` with
  Rust 1.92.0; Node 24.21.0 runs the Wasm command (`cert_toolchain.json`, `losat_bin.sha256`).
- Command, from the repository root at `C`:
  `python3 LOSAT/tests/certify_integrated_runtime_v010.py --expected-sha C --output-dir /tmp/<empty>
  --native-bin LOSAT/target/release/losat --wasm-bin LOSAT/target/serial-command/wasm32-wasip1/release/losat.wasm
  --node <node 24.21.0> --blastn-oracle ... --blastp-oracle ... --tblastx-oracle ... --source-archive
  ncbi-blast-2.17.0+-src.tar.gz` (Python 3.14).

Result `INTEGRATED_RUNTIME_CERTIFIED` (`summary.json`):

| Part | Cases | Result |
| --- | --- | --- |
| Native against NCBI (`native_contracts.tsv`) | 43: blastn 14, blastp 9, tblastx 20 | blastn 14 and blastp 9 the same text as NCBI; tblastx p01-p14 the same text, d01-d06 the approved difference of PD-TBLASTX-LOCAL-SUBJECT-GENCODE (`HSP_SET_DIFF`) |
| Serial Wasm against native (`native_wasm_equalities.tsv`) | 41: blastn 14, blastp 7, tblastx 20 | all the same bytes |
| blastp at 1 and 4 threads (`blastp_native_thread_equalities.tsv`) | 2 | the same bytes |
| Repeats (`repeatability.tsv`) | 6 classes x 2 targets x 3 runs | all the same |
| Runtime boundary (`zero_delegation.stdout`) | | no temporary production finding, the 2 reviewed dependency signals |

Files: the script's summary files above, `certify.log.gz` (its log), `output_manifest.tsv` (path, bytes and
SHA-256 of each of the 618 files the run wrote, 83 MB in all). The whole output is kept as
`A_output.tar.gz` (40,527,413 bytes, SHA-256 `7411441a24b203f4f71aa214bb065a30ac735667ca018339f77950cda94c1b99`).
`MANIFEST.tsv` lists the files of this directory; its SHA-256 is `local_evidence.A.evidence_sha256` in
`docs/release/v0.2.0_rc_contract.json`.
