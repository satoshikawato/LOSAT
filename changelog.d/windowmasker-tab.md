- LOSAT Web: added a WindowMasker main tab after Results for running `windowmasker` in the browser
  across all three modes: Make counts (`-mk_counts`), Mask (`-ustat`), and Convert counts
  (`-convert`), driven by the engine's `-xmlhelp` option descriptions (`ToolEngine.describeTool`).
- LOSAT Web: WindowMasker runs share the single run queue with BLAST searches, stream outputs and
  inflated compressed inputs through session OPFS without holding full files in page memory, honor
  `-out` / `-logfile` file sharing, and report outputs, previews (512 KiB cap; `Binary, <size>` for
  binary outputs), downloads, warnings, and Run details in Results.
- LOSAT Web: counts files made in the session can be sent straight to BLASTN (`Use in BLASTN`,
  attaching `-window_masker_db`) or to Mask (`Mask with these counts`), and BLASTN's repeats row
  shows a `Make counts…` link when `Counts file` is chosen.
