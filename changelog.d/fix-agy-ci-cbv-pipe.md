- Test-only: `cleanup_blastdb_volumes_fixtures::unlink_failure_flushes_stdout_before_stderr_and_matches_closed_pipe_exit_status`
  creates its closed stdout pipe with `std::io::pipe()` and drops the read end before spawning
  `losat` directly instead of racing a `bash` subshell (`exec 0<&-`), eliminating the load-dependent
  `nvol=10` exit-code flake (44 of 100 loaded runs failed before; 200 of 200 pass after). No
  program output changes.
