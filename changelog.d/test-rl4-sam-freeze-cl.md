- Test tooling: `db_search_fixtures.py freeze` keeps a frozen SAM (`-outfmt 17`) file whose NCBI bytes differ only in
  the `@PG` CL tag (NCBI's argv[0], the path of the NCBI build that ran), the rule `check` compares with. A freeze
  with the oracle no longer rewrites the 32 SAM rows frozen from other NCBI 2.17.0 builds. No product change.
