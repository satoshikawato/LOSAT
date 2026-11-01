- LOSAT Web: the HSP records of a search share the alignment rows with the reports instead of copying them
  (`ReportOutputs::hits` takes `Arc<Vec<PairwiseHit>>`), and BLASTN's rows are their exact size. The guest memory peak of a
  megablast search of 100 reads against 28,000 16S sequences falls from 472 to 259 MiB; the linear memory does not shrink
  between searches, so this stays lower for the whole page.
