- CI: the runs with every performance switch off (`LOSAT_X_DEFAULT=0`) are gone: the `rust` cargo-test step
  and the fast regressions in `ci.yml`, and the `reference` mode of the platform checks (outputs and cargo
  test); `LOSAT/tests/reference_step_targets.sh` is deleted. Platform checks compare the default outputs at 1
  and 4 threads and across platforms; the jobs are named `outputs (<platform>)` and `cargo test (<platform>)`.
