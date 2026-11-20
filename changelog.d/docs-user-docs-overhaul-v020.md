- Documentation: README.md is short (what LOSAT is, installation, a quick start, the programs, benchmarks,
  libraries) and links to the user manual for each program's options and differences. The manual gains
  pages for common workflows, environment variables, and the Rust, Python and WebAssembly libraries.
  `benchmarks/v0.2.0` holds the v0.2.0 measurements as TSV files with four figures (standard searches,
  scaling, published tools, WebAssembly) and the scripts that build them. `docs/release/v0.2.0.md` has the
  release's "What's new" section. `python/pyproject.toml` lists the manual, changelog and issue tracker as
  project links. The manual's check (`scripts/check_manual_links.py`) also checks the links and wording of
  the READMEs and the release record.
  The test that every option LOSAT refuses is documented now reads the manual's "What LOSAT refuses"
  page instead of README.md.
