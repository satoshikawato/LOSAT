- Tests: the Web E2E tests use Playwright 1.65.0 (Chromium 158, Firefox 159, WebKit 27.2 build 2388) instead
  of 1.61.1. Its WebKit did not crash once in 15,000 page starts on CI's runners, where WebKit build 2311 crashed
  44 times in the same runs: twice in JIT code on the page's main thread (seen on CI as a run stuck "Preparing" or
  a source stuck "indexing"), the others as processes ended after their test. The OPFS BlockStore contract now runs
  in a profile on disk in Firefox and WebKit, so WebKit runs it for the first time (WebKit refuses OPFS in
  Playwright's incognito contexts, where the application keeps results in memory and says so).
