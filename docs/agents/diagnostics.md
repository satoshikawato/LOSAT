# LOSAT debug and diagnostics environment variables

Moved from `AGENTS.md` on 2026-10-08 so that it is read only when needed.

- `LOSAT_TRACE_HSP="qstart,qend,sstart,send"` trace a specific TBLASTX HSP.
- `LOSAT_TRACE_HSP_MASKS=1` print mask coverage for the traced TBLASTX HSP.
- `LOSAT_TRACE_CHAIN_HSP="qstart,qend,sstart,send"` trace TBLASTX chain
  selection for a specific HSP.
- `LOSAT_TRACE_LINK_SELECTIONS=1` print TBLASTX link-selection details.
- `LOSAT_DUMP_TBLASTX_STAGE=<dir>` append TBLASTX stage snapshots as TSV files.
- `LOSAT_TRACE_BLASTN_HSP="qstart,qend,sstart,send"` trace a specific BLASTN
  HSP.
- `LOSAT_TRACE_BLASTN_SEED="q,s"` trace a specific BLASTN seed.
- `LOSAT_TRACE_BLASTN_CONTEXT=<context_idx>` restrict BLASTN tracing by context.
- `LOSAT_TRACE_BLASTN_SUBJECT=<subject_id_or_index>` restrict BLASTN tracing by
  subject.
- `LOSAT_TRACE_BLASTN_STAGE=<seed|ungapped|prelim|traceback|purge|hitlist|all>`
  restrict BLASTN tracing by stage.
- `LOSAT_DEBUG_CUTOFFS=1` cutoff calculations (tblastx + blastn).
- `LOSAT_DEBUG_CUTOFFS_ALL=1` verbose TBLASTX cutoff diagnostics.
- `LOSAT_DEBUG_CHAINING=1` chaining debug (legacy; tblastx).
- `LOSAT_DEBUG_EXTENSION=1` tblastx extension debug.
- `LOSAT_DEBUG_HSP_SAVING=1` TBLASTX HSP-save diagnostics.
- `LOSAT_DEBUG_OUTPUT_FILTER=1` TBLASTX output filter diagnostics.
- `LOSAT_DEBUG_BLASTN=1` blastn hit loss diagnostics.
- `LOSAT_DEBUG_COORDS=1` blastn coordinate transforms.
- `LOSAT_DEBUG_COORDS_START=<int>` narrow selected BLASTN coordinate diagnostics.
- `LOSAT_DEBUG_SCAN_SOFF=<int>` tblastx scan debug center subject offset.
- `LOSAT_DEBUG_SCAN_WINDOW=<int>` tblastx scan debug window size.
- `LOSAT_TIMING=1` timing breakdown.
- `LOSAT_DIAGNOSTICS=1` general diagnostics counters.
- `LOSAT_STARTUP_TRACE=1` startup trace.
- `LOSAT_WASI_THREADS_DEBUG=1` threaded-WASI scheduling diagnostics.
- `LOSAT_TBLASTX_PARALLEL_CHUNKS=1` force TBLASTX subject-chunk parallel path
  for diagnostics.
- `LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS=1` diagnostic-only sequential TBLASTX
  scan-interior chunking; this does not enable parallel scan work.
