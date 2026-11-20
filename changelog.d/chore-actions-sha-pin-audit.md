- CI: every third-party action in the 11 workflows is pinned to a commit SHA with its release tag as a
  comment (`uses: actions/checkout@<sha> # v4.4.0`), and Nightly gains a `cargo audit` job (RustSec,
  `--deny warnings`, both lock files; the four unmaintained crates that come through `bio` and
  `indicatif` 0.17 are named in its `--ignore` list). Pull-request checks are unchanged.
