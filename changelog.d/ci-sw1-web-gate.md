- CI: only a pull request's own run of web.yml reports the required check `web result`; a push run
  reports `web result (push)`. A push run of a branch that an open pull request covers skips its
  jobs, and its passing gate on the same commit let Auto-merge merge pull requests before their Web
  jobs ended (#843, #844). `@playwright/test` is back at 1.61.1: with 1.64.1 WebKit fails the OPFS
  BlockStore test in a dedicated worker. Dependabot proposes no `@playwright/test` 1.62 or later
  until a newer release, tried by hand, passes that test.
