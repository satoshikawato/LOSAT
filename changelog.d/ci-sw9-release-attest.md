- Release: the v0.2.0 release-readiness handoff adds a CycloneDX SBOM of the dependencies
  (`LOSAT-0.2.0.cdx.json`, one for every target and feature) and signs GitHub artifact attestations
  for each of the eight release files: build provenance and the SBOM. A downloaded file can be checked
  with `gh attestation verify <file> --repo satoshikawato/LOSAT`. Binaries are still not code-signed
  or notarized. No engine change.
