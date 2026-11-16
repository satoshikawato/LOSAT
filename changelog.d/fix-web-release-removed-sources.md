- LOSAT Web: a source that the search form no longer shows (removed, cleared, replaced by newer
  pasted text, replaced by an Example, or a superseded or failed read) is dropped from the Data
  worker with its record tables, and the form's cached checks of it go too; before, the worker kept
  every source until the tab closed, so adding and removing a 32 MB FASTA of 100,000 records ten
  times left 306 MB of heap behind (Chromium), and now 0.5 MB. A run keeps its inputs in its own
  snapshot, as before, so its results, outputs and Run details stay. Searches and outputs are
  unchanged.
