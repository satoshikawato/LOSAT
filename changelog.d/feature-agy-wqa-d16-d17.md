- LOSAT Web:
  - Record exclusions in the query and subject selectors now persist across search mode changes
    (including switching between nucleotide and protein readers); only records newly parsed under
    the new reader enter as included, while clearing, replacing, or removing inputs resets
    exclusions cleanly.
  - Combined inputs and record selections with mixed CR, CRLF and LF line endings now run, read as
    the command line reads the same file; Run details and the Results summary say when that reading
    differs from the record tables of the files read alone.
