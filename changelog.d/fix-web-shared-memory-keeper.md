- Web: in WebKit, searches no longer stop running threaded after some threaded runtimes in one page ("the
  threaded engine could not start: the browser could not reserve 1 GiB of shared memory (Out of memory)"; the
  search then ran on one thread). The Engine worker had sent the main thread the shared memory of each threaded
  runtime, and the main thread kept every memory it had been sent until it collected garbage; a new memory did
  not make it collect. One page of Playwright's WebKit 26.5 on Linux held at most 15 shared memories of 4 GiB, so
  each renewed, cancelled or re-made runtime (for example after a search on one thread) used up one more: in a
  test page, the 23rd to 32nd runtimes fell back to one thread. The memory now stays in the workers of its
  runtime, which end with it. A small worker of the runtime holds it, and grows it to its maximum when the main
  thread ends a runtime during a search, so that WebKit does not crash the page (as before).
