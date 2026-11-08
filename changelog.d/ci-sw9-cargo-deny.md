- CI: `deny.toml` (cargo-deny 0.20.2) for the lock files of LOSAT, web/adapter and python. Every pull
  request checks their licenses (an allow list), sources (crates.io only) and bans (no `*`
  requirement; a second version of a crate is reported) in `rust lint`; Nightly checks their RustSec
  advisories with `cargo deny` instead of `cargo audit`, now for python's lock file too. No engine
  change.
