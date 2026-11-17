- CI: the nightly advisory check (`cargo deny check advisories`) ignores RUSTSEC-2026-0543 (chrono unmaintained, no
  vulnerability) by its id. chrono reaches only `python/Cargo.lock`, through arrow-array (Apache Arrow, pylosat's Arrow
  tables), which requires it in every release; the entry goes when arrow-array drops chrono. No product change.
