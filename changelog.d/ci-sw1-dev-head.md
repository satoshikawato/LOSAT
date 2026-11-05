- CI: a new workflow, Dev head (`dev-head.yml`), starts `ci.yml` on the head of `dev` every four hours
  when that commit has no run that passed or is still running (a merge made by Auto-merge starts no
  push run), and fails while the last run of `ci.yml` on `dev` that ended did not pass. `ci.yml` gains
  `workflow_dispatch` for it. `ci.yml` and `web.yml` declare a read-only token
  (`permissions: contents: read`, the repository's default today); clippy, nextest and the doc tests
  run with `--locked`; a test shard stops after 80 minutes. The checks themselves are unchanged.
