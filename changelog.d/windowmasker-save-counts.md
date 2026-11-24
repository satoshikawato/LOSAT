- LOSAT Web: WindowMasker counts written to an `-out` file by `-mk_counts` or `-convert` can now
  be saved for a taxonomy ID (`Save for taxonomy ID`) in the browser's OPFS counts store
  (`windowmasker/<taxid>/wmasker.obinary`, streaming in 8 MiB slices under the ID's Web Lock),
  with a native confirmation dialog before replacing an existing copy, in-place refusal when
  another tab is downloading that ID, cross-tab synchronization, `Made in this browser` status in
  BLASTN's taxonomy ID row, `Use in BLASTN` switching to `-window_masker_taxid`, and saved counts
  listed in WindowMasker's `-ustat` and `-in` pickers.
- LOSAT Web: WindowMasker tool runs now display a `Supported` / `Not yet supported` status badge
  generated at build time from the NCBI-frozen exit-0 rows of `LOSAT/tests/windowmasker/cases.tsv`.
- LOSAT Web: fixed the BLASTN `Counts file` hint wrapping so the help sentence and `?` stay on a
  full-width line below the file controls, matched tool output download button labels to the tab
  names (`Download output` and `Download warnings`), reference-counted inflated compressed tool
  session files across the form and queued/running runs, cleaned up Data-side staged tool runs on
  cancel during preparing or finalizing, and re-read `ToolDraft` state after awaiting validation
  in `submit()`.
