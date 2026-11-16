- Unused code that no build reaches was removed, or moved under `#[cfg(test)]` where tests use it as a helper:
  the old BLASTN gapped-extension twin with its scratch wrapper, the SIMD sequence-compare helpers, the
  unused BLASTN HSP filter and packed-diagonal key, the common and TBLASTX HSP chaining helpers with their
  tests, unused TBLASTX cutoff and constant definitions, unused command-line value parsers, a
  query-batch helper and a few sequence-database accessors. No user-visible change; outputs do not change.
