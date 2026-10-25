- The trace, debug, dump and profile variables that no test, script or CI step sets (the TBLASTX
  and BLASTP HSP traces, the `LOSAT_DEBUG_*` chaining, extension and scan families, the stage dump, the DP
  capture, the Newton and composition-adjustment counters, `LOSAT_STARTUP_TRACE`; list: docs/agents/diagnostics.md, "Diagnostic build") now need a build with
  `--features diagnostics`; a default build prints one stderr line for each one that is set and ignores it.
  `LOSAT_DEBUG_CUTOFFS`, `LOSAT_DEBUG_BLASTN` and one `LOSAT_X_STATS` read are read once per process. No
  output change.
