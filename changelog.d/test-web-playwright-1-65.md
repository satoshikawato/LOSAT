- Tests: the Web E2E tests use Playwright 1.64.1 (Chromium 156, Firefox 157, WebKit 27.2 build 2370) instead of
  1.61.1, whose WebKit build 2311 crashed in a diagnostic on CI's runners (15,000 page starts with pastes): 44 core
  dumps, 2 of them crashes in JIT code on the page's main thread, which CI saw as a page that stopped answering. The
  OPFS BlockStore contract runs in a profile on disk, where WebKit passes it for the first time (WebKit refuses OPFS in
  Playwright's incognito contexts; that failure had reverted 1.64.1 in #847). Not 1.65.0: its Firefox 159 crashes
  during the axe check of `announce.spec.ts:449` (30 of 30 runs; Firefox 157, 0 of 30).
