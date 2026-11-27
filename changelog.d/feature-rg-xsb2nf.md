- BLASTX's refusal of a `-threshold` whose words fill NCBI BLAST+'s overflow banks of the compressed lookup table
  (`-word_size` 5 and more) in a later query batch now has frozen rows among the tests (tests only; no output or
  behavior changes): two query batches whose second fills the banks at `-threshold 1`, in `-outfmt` 0, 6 and 7, where
  `blastx` writes the first batch's report as NCBI BLAST+ writes it before it crashes and then stops with the refusal
  (exit status 1), and the same two batches at `-threshold 2`, which NCBI BLAST+ and `blastx` run to the end.
