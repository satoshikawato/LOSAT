- Web tests: the Alignments tests no longer fail now and then. They scroll a Range's block (it stays)
  instead of its section (replaced when its text arrives: "Element is not attached to the DOM", 10 CI jobs
  in 7 days), and they take the sections as settled only when every read asked for is answered and shown,
  waiting for the block they scrolled to first (4 CI jobs). The application is unchanged.
