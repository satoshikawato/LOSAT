- LOSAT Web: Cancel ends a run at once while it is "Preparing", and the next run in the queue starts. Before, Cancel
  did nothing until the part of LOSAT that keeps the results answered: seconds later while it read a large input, and
  never if it did not answer, so the run stayed "Preparing" for good. Nothing of the cancelled run is kept.
