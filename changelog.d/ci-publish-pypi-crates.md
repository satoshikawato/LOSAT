- Publication on PyPI and crates.io, by hand on the release's commit (the tag is checked, not checked out):
  the `Python wheels` workflow also builds the sdist of pylosat (and a wheel from it, which runs the frozen
  search) and uploads the five wheels and the sdist by PyPI Trusted Publishing (no stored token) from the
  environment `pypi`, or to TestPyPI as a rehearsal; the new workflow `Publish crate` packages the crate
  `losat`, holds it to the 10 MiB limit of crates.io, compares it byte for byte with the GitHub release's
  `.crate` file, runs `cargo publish --dry-run`, and uploads from the environment `crates-io`. Every pull
  request now checks the crate's size against that limit (`rust lint`).
