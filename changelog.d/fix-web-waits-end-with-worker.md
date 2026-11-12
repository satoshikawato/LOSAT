- Web: in Safari (WebKit), cancelling a search that runs on several threads no longer uses up
  the page's shared memory. A cancelled search's thread whose wait had no timeout stayed alive in
  WebKit, and with it the 4 GiB memory reservation of its search; a page holds about 15, so after
  some cancels a search of several threads ran on one thread ("could not reserve ... shared
  memory"). The shared-memory guard of the threaded module (`LOSAT/tests/wasi_shared_memory.js`,
  version 4) now makes every wait wait in slices of at most 100 ms with a call of the module's
  `sched_yield` import between them, where the browser ends the thread, and the Engine worker's own
  wait for a thread to start waits in the same slices. In Playwright's WebKit, 20 cancels in one page
  kept every search at 4 threads and freed every memory; before, 10 to 15 memories stayed after 12 to
  16 cancels in 3 of 4 pages, and one page ran its 16th search on one thread. Outputs do not change;
  Chromium and Firefox ended such threads already.
