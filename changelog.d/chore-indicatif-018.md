- Dependencies: `indicatif` 0.17 -> 0.18 (0.18.6, with `console` 0.16), which drops the unmaintained
  `number_prefix` (RUSTSEC-2025-0119) from both lock files; the Nightly `cargo audit` job no longer
  ignores it. No output change.
