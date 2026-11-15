- Engine: the composition-based-statistics redo of a BLASTX query batch's matches (the
  read of the matches, the redo ahead of time on the search pool, the heaps and their pour-out)
  moves out of `batch_end` into its own module. No output changes.
