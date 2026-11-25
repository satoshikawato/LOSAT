- LOSAT Web:
  - Record exclusions in the query and subject selectors now persist across search mode changes
    (including switching between nucleotide and protein readers); only records newly parsed under
    the new reader enter as included, while clearing, replacing, or removing inputs resets
    exclusions cleanly.
  - Combined inputs and record selections with mixed CR/CRLF/LF line endings no longer fail with a
    record-table mismatch against the engine: pure-CR joins preserve CR line endings, and assembled
    run inputs scan their exact bytes with the engine scanner so the run and exported input file
    match the native CLI.
