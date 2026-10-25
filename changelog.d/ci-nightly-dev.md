- CI: the scheduled nightly run also starts a nightly run of `dev` (until now it tested only `main`)
  and fails when the last completed nightly run of `dev` did not pass.
