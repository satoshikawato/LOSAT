- `LOSAT_X_DEFAULT` is removed: a run ignores it, and each switch of the strict set is off only with its own
  `LOSAT_X_<NAME>=0` (as before, a diagnostic aid). Entries above that name `LOSAT_X_DEFAULT=0` describe the
  state when they merged. Outputs do not change. The tests compare each switch with its own `=0`; arms that
  only tested `LOSAT_X_DEFAULT` went.
