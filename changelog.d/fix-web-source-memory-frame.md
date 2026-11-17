- Web tests: the test that adds and removes a 32 MB source ten times (`source-memory.spec.ts`) reads the
  JavaScript heap once the page has drawn two frames. The click on Remove starts a CSS transition of the
  button, and Chromium keeps a removed element's transition, and with it the removed card and its source's
  record table (24.5 MB), in the page's animation timeline until the next frame; on CI that frame sometimes
  came after the heap was read (2 of 41 CI jobs failed on 2026-11-15/16). The application and the test's
  limit are unchanged.
