- CI: pull requests into `dev` run the `rust` tests on the reference paths (`LOSAT_X_DEFAULT=0`) again, about 2 minutes
  per shard since the tests run in two cargo-nextest shards at opt-level 3; a test that failed only there turned `dev`
  red after its merge. The reference fast regressions stay on the pull request into `main` and on pushes. The
  `cargo test (k/2)` jobs save the dependency cache also when their tests fail.
