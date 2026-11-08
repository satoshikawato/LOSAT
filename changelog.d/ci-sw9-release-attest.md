- Release: the v0.2.0 release-readiness handoff adds a CycloneDX SBOM of the dependencies
  (`LOSAT-0.2.0.cdx.json`, one for every target and feature) and two GitHub artifact attestations:
  the build provenance of the eight release files and the SBOM, and the SBOM attestation of the
  eight files. A downloaded file can be checked with
  `gh attestation verify <file> --repo satoshikawato/LOSAT`. The workflow is dispatched on
  `release/0.2`, whose head must be the candidate. Binaries are still not code-signed or notarized.
  No engine change.
