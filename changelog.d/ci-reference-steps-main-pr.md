- CI: pull requests into `dev` no longer run the `LOSAT_X_DEFAULT=0` steps of `ci.yml` (cargo test on the reference
  paths, about 7-12 minutes of the `rust` job, and the reference fast regressions, about 2 minutes); the pull request
  into `main`, pushes that start `ci.yml` and the platform checks still run them. The default run's reach check
  (`LOSAT_X_REACH`) requires that it takes the reference side of every conditional switch that has a fixture for it.
