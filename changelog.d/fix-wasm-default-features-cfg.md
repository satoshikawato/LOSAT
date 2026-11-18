- Build: `cargo build --target wasm32-wasip1` with the default features compiles again (the
  `-db` stream's index-length call was left out of that build by an earlier narrowing; shipped
  builds were not affected).
