- Web tests: the storage tests start their search only once the form stops moving (the describe that
  comes after the page opens moved "Run LOSAT", and a click across the move started no run: 7 CI jobs in
  7 days), and the engine smoke test gets 120 s instead of 30 s (WebKit on CI took 10-23 s to open the
  page alone: 3 CI jobs). The application is unchanged.
