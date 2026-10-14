# Authority v0.1.1: Linux arm64 characterization

`LOSAT/tests/ncbi_platform_variance_v011.json` (`ncbi-platform-variance-v0.1.1`,
`PD-NCBI-PLATFORM-VARIANCE` 1.2) is authority v0.1.0 plus one platform,
`linux-aarch64` (GitHub `ubuntu-24.04-arm`, NCBI BLAST+ 2.17.0
`ncbi-blast-2.17.0+-aarch64-linux.tar.gz`, published MD5
`5f9d49b999346e0ddb6637b71af4b22a`).

## Characterization run

- Run `38563517687` of `ncbi-only-platform-characterization.yml`, dispatched at
  candidate `d81f0209a5d91af9a820832d8ef56e7a1cbe621e` (branch
  `characterize/linux-aarch64`, not merged, as the v0.1.0 candidate
  `be672531` was not). The candidate is v0.1.0's analyzer commit `cfd213f5`
  with the platform set of `characterize_ncbi_platform_variance_v010.py`
  replaced by `linux-aarch64`, the driver's `linux-aarch64` PlatformSpec and an
  ELF reader in `binary_architecture`, and nothing else.
- Same fixtures (PR 6 V3 Git blobs, checked against the input identities
  authority v0.1.0 binds), same retained Linux reference (exported, not
  rerun), same commands; 6 authoritative + 6 diagnostic official NCBI
  searches, no LOSAT execution.
- `classify-added-platform` wrote `linux-aarch64-authority-record.json` (this
  folder, byte-identical to the run's artifact
  `ncbi-characterization-added-platform-record`; SHA-256
  `6600253b8cbb2ab4d984b79d571b54b070f341cd895420a7c4e34fbba5772d39`, bound as
  `derived_record_sha256`). Its classes follow the rule that reproduces every
  v0.1.0 diagnostic: identical bytes, then identical after CRLF/CR to LF, then
  the same row keys in another order, then the same keys with other fields.

## Result

| representative | Linux arm64 raw SHA-256 | against the retained Linux x64 reference |
| --- | --- | --- |
| blastn PesePMNV.MjPMNV.task_blastn | `fcd9236c…` | identical |
| blastn Sakai.MG1655.megablast | `ac817776…` | identical |
| blastp pairwise_default_serial | `fd4b0108…` | identical |
| tblastx p03_mela_pemojnva | `2533575e…` | row order (114 rows), new fingerprint |
| tblastx p14_ap027131_ap027133_query4 | `2a07a323…` | row order (586 rows), new fingerprint |
| tblastx d06_ap027131_ap027133_db4 | `7d3086ea…` | row order (422 rows), = macOS arm64's |

Executables (ELF, arm64): blastn `9f1c3376…`, blastp `28cb6cc9…`, tblastx
`f278c296…`. A probe of the full 61-step plan on the same runner image
(run `38549611962`, branch `bench/aarch64-cert`) gave the same six
fingerprints and executable hashes.

## Rebuild

```sh
python3 docs/evidence/ncbi_platform_variance_v011/build_authority_v011.py \
  --v010 LOSAT/tests/ncbi_platform_variance_v010.json \
  --record docs/evidence/ncbi_platform_variance_v011/linux-aarch64-authority-record.json \
  --run-id 38563517687 --output /tmp/v011.json
cmp /tmp/v011.json LOSAT/tests/ncbi_platform_variance_v011.json
```

`test_certify_platform_native_v010.py` checks that v0.1.1 carries v0.1.0's
representatives, three platform records, 18 diagnostics and characterization
provenance unchanged, that the Linux arm64 parts equal this record, and that
the v0.1.0 file is unedited.
