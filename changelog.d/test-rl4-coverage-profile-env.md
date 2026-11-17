- Test tooling: four tests that failed only in the nightly coverage jobs pass there: `ncbi_toolkit_dryrun` keeps
  `LLVM_PROFILE_FILE` for the command it runs with a cleared environment (the profile no longer lands in the case's
  directory), and `tax_report_streams` does not compare the profile runtime's "LLVM Profile Error" line that the file
  size limit causes. No product change.
