- CI: rust caches are saved only by runs on dev and main, and pull requests restore dev's copy;
  their own copies, which no other pull request can read, had filled the 10 GiB of Actions cache.
  The semver report job runs without a cache. python.yml can be started by hand, and its runs on
  dev and main are no longer cancelled by a newer push.
