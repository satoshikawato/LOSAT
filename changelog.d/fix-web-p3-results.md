- LOSAT Web: Cancel stops a run that is still preparing its outputs; before, it waited until the
  search had ended. When the browser refuses to keep the screen on, the refusal stays shown for the
  rest of the runs instead of being asked again at each step of a run. The E value and bit score
  filters compare the values that the lists show, so "E value ≤ 0" keeps the HSPs written as 0.0
  (NCBI writes 0.0 below 1e-180). The Dot Plot's HSP popup opens right of the plot, or under it where
  there is no room, so that it no longer covers other HSPs' lines. A query without HSPs says why in
  the Dot Plot tab too ("No significant similarity found for this query."). "Download All" lists
  outfmt 0, 6 and 7 and saves the one chosen, and saved outputs are named with the task and the
  time of the run (`losat-run3-blastn-megablast-20261104T111530Z.outfmt6.txt`).
